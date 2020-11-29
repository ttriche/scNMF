#' Non-negative matrix factorization with alternating least squares with sequential coordinate descent against mean squared error loss
#' Methods are adapted from Lin and Boutros, 2018 BMC Bioinformatics with several modifications for speed gains
#' Function is theoretically equivalent to NNLM::nnmf(..., method = "scd", loss = "mse", alpha = c(0,0,0), beta = c(0,0,0), mask = NULL, init = NULL)
#' This method offers >2x faster calculation than NNLM::nnmf on a single thread for low-rank factorization of a 10k x 200 matrix, and increasingly faster with matrix size, number of threads, and rank
#' Key differences from NNLM::nnmf are 1) an improved OpemMP multithreading loop structure, 2) a parallelized MSE error loss function, 3) fewer conditional checks, 4) no MKL error calculations, and 5) no support for regularization, masking, or non-random initialization
#' 
#' @param A A matrix to be factorized, in dense format.
#' @param k Decomposition rank, integer.
#' @param rel.tol Stop criterion, defined as the relative tolerance between two successive iterations: |e2-e1|/avg(e1,e2). Default 1e-3.
#' @param inner.rel.tol Stop criterion for the inner sequential coordinate descent least squares loop, defined as relative tolerance passed to inner W or H during matrix updating: |e2-e1|/avg(e1,e2). Default 1e-6.
#' @param max.iter Maximum number of alternating NNLS solutions for H and W, integer. Default 1000.
#' @param inner.max.iter Maximum number of iterations passed to each inner W or H matrix updating function. Default 1e-6.
#' @param verbose 0 = no tracking, 1 = progress bar relative to max.iter, 2 = iteration info. Default 2.
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores).
#' @param trace An integer indicating how frequently the MSE error should be calculated and checked for convergence. To check error every iteration, specify 1. To avoid checking error at all, specify trace > max.iter.  Default 5.
#' @return A list of W and H matrices
nnmf <- function(A, k = NULL, max.iter = 1000, rel.tol = 1e-3, n.threads = 0, verbose = 2, inner.max.iter = 100, inner.rel.tol = 1e-6, trace = 5) {
    if(n.threads < 0) stop("Specify 0 or a positive integer for n.threads")
    if(is.null(k)) stop("Specify a positive integer value for k")
    if(rel.tol > 0.1) warning("rel.tol is greater than 0.1, results may be unstable")
    if(inner.rel.tol > 1e-4) warning("inner.rel.tol is greater than 1e-5, it may take unnecessarily long to converge. Consider a smaller inner.rel.tol.")
    if(trace < 1) stop("trace must be a positive integer")
    if(inner.max.iter < 50) warning("inner.max.iter < 50 is not recommended")
	if (is.logical(verbose)) verbose <- as.integer(verbose)
	res <- c_nnmf(A, as.integer(k), as.integer(max.iter), as.double(rel.tol), as.integer(n.threads), as.integer(verbose), as.integer(inner.max.iter), as.double(inner.rel.tol), as.integer(trace))
    rownames(res$W) <- rownames(A)
    colnames(res$H) <- colnames(A)
    colnames(res$W) <- rownames(res$H) <- paste0("NMF_",1:k)
    return(res)
}