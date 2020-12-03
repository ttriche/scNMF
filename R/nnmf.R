#' Fast non-negative matrix factorization by alternating least squares with sequential coordinate descent against mean squared error loss
#' Methods are adapted from Lin and Boutros, 2018 BMC Bioinformatics with several improvements for efficiency in unregularized and non-masked use cases
#' Function is theoretically equivalent to NNLM::nnmf(..., method = "scd", loss = "mse", alpha = c(0,0,0), beta = c(0,0,0), mask = NULL, init = NULL)
#' This method offers >2x faster calculation than NNLM::nnmf on a single thread for low-rank factorization of a 10k x 200 matrix, and increasingly faster with matrix size, number of threads, and rank
#' Key differences from NNLM::nnmf are 1) an improved OpemMP multithreading loop structure, 2) a parallelized MSE error loss function, 3) fewer conditional checks, 4) no MKL error calculations, and 5) no support for regularization, masking, or non-random initialization
#'
#' @param A A matrix to be factorized. If sparse, will be converted to dense
#' @param k Decomposition rank, integer (required)
#' @param rel.tol Stop criterion, defined as the relative tolerance between two successive iterations: |e2-e1|/avg(e1,e2). (default 1e-3)
#' @param verbose boolean, give updates every trace iterations
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores).
#' @param max.iter Maximum number of alternating NNLS solutions for H and W, integer (default 1000)
#' @param trace An integer specifying a multiple of iterations at which MSE error should be calculated and checked for convergence. To check error every iteration, specify 1. To avoid checking error at all, specify trace > max.iter (default is 5, and is generally an efficient and effective value)
#' @return A list of W and H matrices
nnmf <- function(A, k = NULL, max.iter = 1000, rel.tol = 1e-3, n.threads = 0, verbose = TRUE, trace = 5) {
  if (class(A)[1] == "dgCMatrix") A <- as.matrix(A)
  if (is.null(k)) stop("Specify a positive integer value for k")
  if (rel.tol > 0.1) stop("rel.tol is greater than 0.1. Results will be unstable and unreproducible")
  if (trace < 1) trace <- 1
  if (n.threads < 0) n.threads <- 0
  inner.max.iter <- 100
  inner.rel.tol <- 1e-6

  res <- c_nnmf(
    A,
    as.integer(k),
    as.integer(max.iter),
    as.double(rel.tol),
    as.integer(n.threads),
    as.integer(verbose),
    as.integer(inner.max.iter),
    as.double(inner.rel.tol),
    as.integer(trace)
  )

  rownames(res$W) <- rownames(A)
  colnames(res$H) <- colnames(A)
  colnames(res$W) <- rownames(res$H) <- paste0("NMF_", 1:k)
  return(res)
}