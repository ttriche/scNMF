#' Run recursive bipartite clustering by spectral decomposition with sparse matrix methods
#'
#' @param data Seurat object, the default assay counts slot will be used
#' @param stop.modularity Optional, stopping criteria for modularity (default = NULL, modularity not used as stopping criteria). Most logical stopping criteria (if used) is 0. Newman-Girvan modularity ranges from -1 to 1.
#' @param stop.delta.dist Optional, stopping criteria based on distance. Faster calculation than modularity, measures the mean fold change in similarity for all cells between a child cluster center and a parent cluster center using cosine similarity.
#' @param stop.min.cells Required, the minimum number of cells permitted in a cluster. Real valued between 1 and ncol(A). Default 50.
#' @param subsample.size The maximum number of cells to sample for determining cluster centers by ssvd (default 20000)
#' @param Q.sample.size The maximum number of cells in each clustering for which to compute a cosine similarity graph for calculation of Q (only applicable if stop.modularity is set)
#' @param seed Random seed for subsampling, set to NULL by default
#' @returns A Seurat object where `object@ident` has been updated with new cluster info. Individual ident classes are created and stashed for each generation (i.e. seuratObj@meta.data$gen1.idents). Clusters are given names as a sequence of 1 and 2 corresponding to bipartite splits, and leaf clusters names are terminated with a 0
FindClusters.Divisive <- function(seuratObj, stop.modularity = NULL, stop.delta.dist = NULL, stop.min.cells = 50, subsample.size = 20000, seed = NULL, verbose = 1, shift = FALSE, Q.sample.size = 1000){
  library(Matrix)
  library(qlcMatrix)
  library(Seurat)
  library(sparsesvd)

  d <- Matrix(GetAssayData(object = seuratObj, slot = "counts"), sparse = TRUE) # load data from the "counts" slot of the default assay of the given Seurat object as a sparse Matrix

  # run bipartite.cluster on the initial dataset
  dist <- Q <- FALSE
  if(!is.null(stop.delta.dist)) dist <- TRUE
  if(!is.null(stop.modularity)) Q <- TRUE
  if(verbose > 0) message("\nRunning divisive clustering on input dataset")
  first.split <- cluster.divisive(d, subsample.size = subsample.size, seed = seed, Q = Q, dist = dist, Q.sample.size = Q.sample.size, shift = shift)
  nodes <- list("1" = names(first.split$clusters[first.split$clusters == 1]), "2" = names(first.split$clusters[first.split$clusters == 2]))
  prev.iter.nodes <- c()
  if(!is.null(stop.modularity) && first.split$Q < stop.modularity){
    warning("Modularity of first two clusters was less than stop.modularity.\n")
  } else if(!is.null(stop.delta.dist) && first.split$delta.dist < stop.delta.dist){
    warning("Change in distance between cells to clusters vs parent was less than stop.delta.dist.\n")
  }else if(ncol(d) < stop.min.cells*2) {
    warning("Number of cells in input seuratObj are not more than twice the value of stop.min.cells, thus the object cannot be further divided.")
  } else{
    prev.iter.nodes <- c("1","2")
  }

  while(length(prev.iter.nodes > 0)){
      iter.nodes <- c()
      if(verbose > 0) message(paste0("\nRunning divisive clustering on ", length(prev.iter.nodes)," clusters from generation ", nchar(prev.iter.nodes[1])))
      if(verbose == 1) pb <- txtProgressBar(char = '=', style = 3, max = length(prev.iter.nodes))
      for(i in 1:length(prev.iter.nodes)){
        # if the node has fewer than stop.min.cells, just make it a leaf node right away
        # if either of the child nodes meet stop criteria, make the node a leaf node
        stop <- FALSE
        node.name <- prev.iter.nodes[i]
        node.cells <- nodes[[node.name]]
        if(stop == FALSE){
            node.children <- cluster.divisive(d[,nodes[[node.name]]], subsample.size = subsample.size, seed = seed, Q = Q, dist = dist, Q.sample.size = Q.sample.size, shift = shift)
            node1.cells <- names(node.children$clusters[node.children$clusters == 1])
            node2.cells <- names(node.children$clusters[node.children$clusters == 2])
            if(!is.null(stop.modularity) && node.children$Q < stop.modularity) stop <- TRUE
            if(!is.null(stop.delta.dist) && node.children$delta.dist < stop.delta.dist) stop <- TRUE
            if(min(length(node1.cells),length(node2.cells)) < stop.min.cells) stop <- TRUE
        }
        if(length(node.cells) < (2 * stop.min.cells)) stop <- TRUE
        if(stop == FALSE){
          if(length(node1.cells) > (2*stop.min.cells)){
            iter.nodes <- c(iter.nodes, paste0(node.name,"1"))
            nodes[[paste0(node.name,"1")]] <- node1.cells
          } else{
            nodes[[paste0(node.name,"10")]] <- node1.cells
          }
          if(length(node2.cells) > (2*stop.min.cells)){
            iter.nodes <- c(iter.nodes, paste0(node.name,"2"))
            nodes[[paste0(node.name,"2")]] <- node2.cells
          } else{
            nodes[[paste0(node.name,"20")]] <- node2.cells
          }
          if(verbose > 1) {
            txt <- paste0("   ...split cluster", node.name)
            if(!is.null(stop.modularity)) txt <- paste0(txt, ", Q: ", trunc(node.children$Q,5))
            if(!is.null(stop.delta.dist)) txt <- paste0(txt, ", delta.dist: ", trunc(node.children$delta.dist,5))
            txt <- paste0(txt, ", num.cells: ", length(node1.cells), " " , length(node2.cells))
            txt <- paste0(txt,"\n")
            cat(txt)
          }
        }
        if(stop == TRUE){
          # add a 0 to the node name to designate it as a leaf node, add the cluster center to the centers matrix
          if(verbose > 1) cat("   ...could not divide cluster",node.name,"\n")
          nodes[[paste0(node.name,"0")]] <- nodes[[node.name]]
          nodes[[node.name]] <- NULL
        }
        if(verbose == 1) setTxtProgressBar(pb = pb, value = i)
      }
      prev.iter.nodes <- iter.nodes
    }
    if(verbose > 0) message("\nUpdating Seurat object with new identities")

  # assign cells to Idents in the input SeuratObj
  num.generations <- max(sapply(names(nodes), nchar))-1   # find number of generations (length of longest node.name)
  if(verbose == 1) pb <- txtProgressBar(char = '=', style = 3, max = num.generations)
  for(i in 1:num.generations){  # assign idents by generation, permitting cells with no membership in current generation to maintain their leaf node membership
    nodenames <- unlist(lapply(names(nodes), function(x) if(nchar(sub("0","",x)) == i) x)) # get all nodes of length i (after replacing zeros)
    for(j in 1:length(nodenames)) Idents(seuratObj, cells = nodes[[nodenames[j]]]) <- nodenames[j]
    seuratObj[[paste0("gen",i,".idents")]] <- Idents(seuratObj)
    if(verbose == 1) setTxtProgressBar(pb = pb, value = i)
  }
  seuratObj[["leaf.idents"]] <- Idents(seuratObj)
  return(seuratObj)
}
