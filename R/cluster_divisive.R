#' A single instance of divisive clustering by sparse SVD
#'
#' @param d A sparse dgCMatrix with genes in rows and cells in columns
#' @param subsample.size Number of cells to subsample for determination of cluster centers by SVD
#' @param seed Random seed for subsampling
#' @param shift Boolean, shift the clustering vector to achieve optimal bipartitioning about 0 (shift determined by the mean of both kmeans centers on the bimodal density distribution of the clustering vector)
#' @param Q Boolean, calculate Newman-Girvan modularity
#' @param Q.sample.size How many cells at most to use for calculating the modularity. A square cosine similarity matrix of this size is computed. Default 1000.
#' @param dist Whether to calculate the change in distance (delta.dist) between cells and their respective child cluster center relative to the parent center
#' @returns A list of cells in either cluster and the dip test result for bimodality between cluster centers
cluster.divisive <- function(d, subsample.size = 5000, seed = NULL, shift = TRUE, dist = TRUE, Q = FALSE, Q.sample.size = 1000){
  require(qlcMatrix)
  require(sparsesvd)
  res <- list()

  if(Q.sample.size > 5000) warning("Setting a large Q.sample.size (i.e. > 2000) can result in large memory load and slow execution time. Consider using only 1000-2000 cells at most to construct a similarity matrix for calculation of Q.")
  if(Q.sample.size < 250) warning("Setting too small of a Q.sample.size (i.e. < 1000) can result in inaccurate calculations of Q due to too much subsetting. A small value should only be used if properly set in proportion to your dataset size.")
  if(subsample.size < 1000) warning("Subsampling fewer than 1000 cells can result in non-representative or low-confidence clusterings. Consider using between 3-5000 cells at least for calculating clusters, and more if your dataset is very large")

  # run sparse svd on the input data
  if(subsample.size > 0 && ncol(d) > subsample.size){
    if(!is.null(seed)) set.seed(seed)
    d.sub <- d[,sample(1:ncol(d), subsample.size)]
    c <- sparsesvd(d.sub, rank = 2)$v[,2]
    parent.center <- Matrix::rowMeans(d.sub)
    names(c) <- colnames(d.sub)
  } else{
      c <- sparsesvd(d, rank = 2)$v[,2]
      parent.center <- Matrix::rowMeans(d)
      d.sub <- NULL
      names(c) <- colnames(d)
  }

  res[["density"]] <- density(c)

  # since 0 isn't exactly the minima in a bimodal density, find the center and shift the clustering vector so the best partition is about 0
  if(shift == TRUE){
      shift.by <- mean(kmeans(c, centers = 2)$centers)
      res[["shift"]] <- shift.by
      c2 <- c - shift.by
      # only shift if it doesn't move all cells to the right or left of zero
      if(length(c2[c2 > 0]) > 1 && length(c2[c2 < 0]) > 1) c <- c2
  }

  # cluster centers are the rowmeans of all cells in that cluster (where clusters are defined by the sign of the (shifted) singular vector)
  c1.ncells <- length(names(c[c > 0]))
  c2.ncells <- length(names(c[c < 0]))

  if(c1.ncells > 1 && c2.ncells > 1){
    if(!is.null(d.sub)) {
      res[["centers"]] <- cbind(Matrix::rowMeans(d.sub[,names(c[c > 0])]), Matrix::rowMeans(d.sub[,names(c[c < 0])]))
    } else{
      res[["centers"]] <- cbind(Matrix::rowMeans(d[,names(c[c > 0])]), Matrix::rowMeans(d[,names(c[c < 0])]))
    }
  } else{
    if(c1.ncells == 1){
      if(!is.null(d.sub)) {
        center1 <- d.sub[,names(c[c > 0])]
      } else{
        center1 <- d[,names(c[c > 0])]
      }
    } else{
      if(!is.null(d.sub)) {
        center1 <- Matrix::rowMeans(d.sub[,names(c[c > 0])])
      } else{
        center1 <- Matrix::rowMeans(d[,names(c[c > 0])])
      }
    }
    if(c2.ncells == 1){
      if(!is.null(d.sub)) {
        center2 <- d.sub[,names(c[c < 0])]
      } else{
        center2 <- d[,names(c[c < 0])]
      }
    } else{
      if(!is.null(d.sub)) {
        center2 <- Matrix::rowMeans(d.sub[,names(c[c < 0])])
      } else{
        center2 <- Matrix::rowMeans(d[,names(c[c < 0])])
      }
    }
    res[["centers"]] <- cbind(center1, center2)        
  }


  is.singleton <- FALSE
  if(length(names(c[c>0])) == 1 || length(names(c[c < 0])) == 1) is.singleton <- TRUE
  # calculate the change in cosine distance of cells from parent cluster center to child cluster centers relative to distance to parent cluster center
  if(dist == TRUE && is.singleton == FALSE){
    if(!is.null(d.sub)){
      rel.dist1 <- mean(apply(1 - cosSparse(cbind("center" = res$centers[,1], "parent" = parent.center), d.sub[,names(c[c > 0])]), 2, function(x) (x[2] - x[1])/x[2]))
      rel.dist2 <- mean(apply(1 - cosSparse(cbind("center" = res$centers[,2], "parent" = parent.center), d.sub[,names(c[c < 0])]), 2, function(x) (x[2] - x[1])/x[2]))
    } else{
      rel.dist1 <- mean(apply(1 - cosSparse(cbind("center" = res$centers[,1], "parent" = parent.center), d[,names(c[c > 0])]), 2, function(x) (x[2] - x[1])/x[2]))
      rel.dist2 <- mean(apply(1 - cosSparse(cbind("center" = res$centers[,2], "parent" = parent.center), d[,names(c[c < 0])]), 2, function(x) (x[2] - x[1])/x[2]))
    }
    # interpretation: clusters are res$delta.dist fold more similar to their assigned clusters than to the parent cluster
    res[["delta.dist"]] <- mean(rel.dist1, rel.dist2)
  } else if(dist == TRUE && is.singleton == TRUE){
      res[["delta.dist"]] <- 0
  }

  if(!is.null(d.sub)) {
      # use cluster centers from subset of d to match all cells to their best-fit cluster based on cosine similarity
      c <- apply(cosSparse(res$centers, d), 2, which.max)
      d <- d.sub
  } else{
      # assign cluster centers from all of d based on the sign of their coefficient in the (shifted) singular vector
      c[c > 0] <- 1
      c[c < 0] <- 2
  }

  res[["clusters"]] <- c

  # calculate Newman-Girvan modularity between both clusters from at most 1000 randomly selected cells from d
  if(Q == TRUE){
    if(ncol(d) > Q.sample.size) d <- d[,sample(1:ncol(d), Q.sample.size)]
    A <- cosSparse(d)

    # remove diagonal weights as these are not true "edges" in the adjacency graph
    diag(A) <- 0
    m <- intersect(names(c[c == 1]), rownames(A))

    # L is the sum of all edges in the adjacency graph
    L <- sum(A)

    # Okk is the sum of all edges within cluster1
    Okk <- sum(A[m,m])

    # Lk is the sum of all edges for cells in cluster1, both within and between cluster1 and cluster2
    ifelse(!is.null(ncol(A[m,])), Lk <- sum(rowSums(A[m,])), Lk <- sum(A[m,]))

    # Q is the difference between the observed connectivity density among nodes in cluster1 (Okk/L) and
    # the expectation of the connectivity density when connections are randomly assigned to the pairs of nodes given the node degrees  
    Q1 <- ((Okk/L) - (Lk/L)^2)
    
    # Because Q is the same for cluster 1 and cluster 2, we can simply multiply Q * 2 to get the final Newman modularity
    res[["Q"]] <- 2*Q1
  }

  return(res)
}
