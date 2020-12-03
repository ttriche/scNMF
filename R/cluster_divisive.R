#' Fast divisive clustering of columns in a large sparse matrix
#'
#' @param input A dgCMatrix matrix or a Seurat object with counts for clustering in the default assay counts slot
#' @param min.samples The minimum number of samples permitted in a cluster. Real valued between 2 and ncol(input), defaults to 10
#' @param subsample Number of samples to subsample for SVD and delta distance calculations, usually only useful for large matrices (> 10000 samples)
#' @param min.dist Minimum delta distance, which is a proportion of how much more similar samples are to their assigned cluster than their opponent cluster. NULL by default. Values which approximate Newman-Girvan modularity of 0 are usually between 0.01 and 0.05.
#' @param seed Random seed for subsampling (NULL by default)
#' @param idents.prefix Prefix for ident names, useful to prevent later runs of cluster.divisive on the same object from overwriting previous idents
#' @param verbose 0 = no output, 1 = output for each generation, 2 = progress bar for each generation, 3 = details for each division
#' @returns If a Seurat object was input, a Seurat object with idents set to cluster IDs. If a sparse matrix was input, a list of cluster assignments for each sample.
cluster.divisive <- function(input, min.samples = 10, subsample = 10000, seed = NULL, verbose = 2, min.dist = NULL, idents.prefix = "") {
  require(Matrix)
  require(wordspace)
  require(irlba)
  require(qlcMatrix)

  if (subsample < 5000) warning("Subsampling too few samples can result in non-representative or low-confidence clusterings. Consider using at least 10000 samples at least for calculating clusters, and more if your dataset is very large or heterogenous. There is very little to no speed gain for subsampling less than 10000. Subsampling is not necessary unless your dataset is larger than at least 100,000 samples and highly heterogeneous.")

  if (!is.null(seed)) set.seed(seed)

  if (class(input)[1] == "Seurat") {
    require(Seurat)

    # load data from the "counts" slot of the default assay of input
    d <- GetAssayData(object = input, slot = "counts")

    # if d is not a sparse matrix, coerce it to one
    if (class(d)[1] != "dgCMatrix") d <- Matrix(d, sparse = TRUE)
  } else if (class(input)[1] == "dgCMatrix") {
    d <- input
  } else {
    stop("input is not a supported class, either dgCMatrix or Seurat")
  }

  prev_iter_nodes <- c("1")
  nodes <- list("1" = colnames(d))
  while (length(prev_iter_nodes > 0)) {
    if (verbose > 1) message("\n")
    if (verbose > 0 && length(prev_iter_nodes) > 1) message(paste0("Running divisive clustering on ", length(prev_iter_nodes), " clusters from generation ", nchar(prev_iter_nodes[1])))
    if (verbose > 0 && length(prev_iter_nodes) == 1) message(paste0("Running divisive clustering on 1 cluster from generation ", nchar(prev_iter_nodes[1])))
    if (verbose == 2) pb <- txtProgressBar(char = "=", style = 3, max = length(prev_iter_nodes), width = 50)
    iter_nodes <- list()
    for (i in 1:length(prev_iter_nodes)) {
      node.name <- prev_iter_nodes[i]
      node.samples <- nodes[[node.name]]

      # no attempt is made to split a cluster if it contains less than twice the minimum number of samples
      if (length(node.samples) < 2 * min.samples) {
        if (verbose > 2) message(paste0("   will not attempt to split cluster ", paste0(node.name), " (", length(node.samples), " samples < ", min.samples * 2, " min.samples * 2)"))
      } else {
        # spectral.cluster is the worker function for cluster.divisive
        # parallelizing with doParallel or future resulted in greater multithreading load than benefit due to multithreading
        curr.node <- spectral.cluster(node.name = node.name, node.samples = node.samples, d = d[, node.samples], min.dist = min.dist, min.samples = min.samples, subsample = subsample)
        if (verbose > 2) message(curr.node$message)
        if (!is.null(curr.node$child_nodes)) iter_nodes <- c(iter_nodes, curr.node$child_nodes)
      }
      if (verbose == 2) setTxtProgressBar(pb = pb, value = i)
    }
    nodes <- c(nodes, iter_nodes)
    prev_iter_nodes <- names(iter_nodes)
  }

  if (class(input)[1] == "Seurat") {
    # assign samples to Idents in the input input
    if (verbose > 1) message("\n\nUpdating Seurat Idents with cluster assignments")

    # find number of generations (length of longest node.name)
    num_generations <- max(sapply(names(nodes), nchar))
    if (verbose > 1) pb <- txtProgressBar(char = "=", style = 3, max = num_generations, width = 50)

    # assign idents by generation, permitting samples with no membership in current generation to maintain their leaf node membership
    for (i in 1:num_generations) {
      # get all nodes of length i (after replacing zeros)
      node_names <- unlist(lapply(names(nodes), function(x) if (nchar(x) == i) x))
      for (j in 1:length(node_names))
        Idents(input, cells = nodes[[node_names[j]]]) <- node_names[j]
      input[[paste0(idents.prefix, "gen", i, ".idents")]] <- Idents(input)
      if (verbose > 1) setTxtProgressBar(pb = pb, value = i)
    }
    input[[paste0(idents.prefix, "leaf.idents")]] <- Idents(input)
    return(input)
  } else {
    return(nodes)
  }
}


spectral.cluster <- function(node.name, node.samples, d, min.dist, min.samples, subsample) {

  # apply subsampling if d is 1.5x larger than subsample
  # note that there is computational overhead to subsampling (non-sampled samples must be matched to clusters)
  # thus, we only apply subsampling if 1.5x larger than subsample
  ifelse(subsample > 0 && length(node.samples) > subsample * 1.5,
    d_node <- d[, sample(1:length(node.samples), subsample)],
    d_node <- d)

  # fast truncated SVD with sparse matrix support using irlba (almost 2x faster than sparsesvd)
  svd_vec <- irlba(d_node, nv = 2)$v[, 2]
  names(svd_vec) <- colnames(d_node)

  # assign cluster idents based on the sign of the second right svd vector
  ident1 <- which(svd_vec > 0)
  ident2 <- which(svd_vec < 0)

  node.centers <- NULL
  item.message <- "continue"

  # sometimes svd produces singlets. Pick these out early and save the hassle of further calculations if this is the case
  if (length(ident1) < 2 || length(ident2) < 2) {
    item.message <- paste0("   could not split cluster ", node.name, " (clustering on subsampling population resulted in a singlet)")
  } else {

    # if d_node was subsetted, calculate cluster centers and match all samples to a center based on cosine similarity
    # in this case, wordspace::dist.matrix is faster than qlcMatrix::cosSparse
    if (ncol(d_node) != length(node.samples)) {
      node.centers <- cbind(Matrix::rowMeans(d_node[, ident1]), Matrix::rowMeans(d_node[, ident2]))
      all.idents <- apply(wordspace::dist.matrix(node.centers, d, byrow = FALSE), 2, which.max)
      ident1 <- which(all.idents == 1)
      ident2 <- which(all.idents == 2)
    }

    # if fewer than min.samples occur in either of the clusters, this cluster cannot be split
    if (length(ident1) < min.samples || length(ident2) < min.samples) {
      item.message <- paste0("   could not split cluster ", node.name, " (", length(ident1), " or ", length(ident2), " samples is < ", min.samples, " min.samples)")
    }
  }

  # measure the mean similarity of samples to their assigned cluster vs the other cluster
  # when clusters become very similar (i.e. low modularity of the bipartition) this value approaches zero
  if (item.message == "continue" && !is.null(min.dist)) {
    if (is.null(node.centers)) {
      node.centers <- cbind(Matrix::rowMeans(d_node[, which(svd_vec > 0)]), Matrix::rowMeans(d_node[, which(svd_vec < 0)]))
    }
    # in this case, qlcMatrix::cosSparse is faster than wordspace::dist.matrix
    delta.dist <- mean(c(
            apply(qlcMatrix::cosSparse(node.centers, d_node[, which(svd_vec > 0)]), 2, function(x)(x[1] - x[2]) / x[1]),
            apply(qlcMatrix::cosSparse(node.centers, d_node[, which(svd_vec < 0)]), 2, function(x)(x[2] - x[1]) / x[2])))
    if (delta.dist < min.dist) {
      item.message <- paste0("   could not split cluster ", node.name, " (", round(delta.dist, 5), " dist is < ", min.dist, " min.dist)")
    }
  }

  if (item.message == "continue") {
    child_nodes <- list()
    child_nodes[[paste0(node.name, "1")]] <- names(ident1)
    child_nodes[[paste0(node.name, "2")]] <- names(ident2)
    txt <- paste0("   split cluster ", node.name)
    if (!is.null(min.dist)) txt <- paste0(txt, ", dist: ", round(delta.dist, 7))
    item.message <- paste0(txt, ", num.samples: ", length(ident1), " ", length(ident2))
    return(list("child_nodes" = child_nodes, "message" = item.message))
  } else {
    return(list("child_nodes" = NULL, "message" = item.message))
  }
}