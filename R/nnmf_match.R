#' Runs NMF on non-overlapping halves of a dataset and computes minimum possible cosine angle between both models
#'
#' @param A A matrix to be factorized. If sparse, will be converted to dense
#' @param k Decomposition rank(s), integer or array of integers. If array, match will be run at multiple ranks of k.
#' @param ind Indices in A to be assigned to one half (optional)
#' @param byrow Whether indices represent rows or columns (rows, by default)
#' @param rel.tol Stop criterion, defined as the relative tolerance between two successive iterations: |e2-e1|/avg(e1,e2). (default 1e-3)
#' @param verbose 0 = no tracking, 1 = progress bar relative to max.iter, 2 = updates every trace iterations (default 2)
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores).
#' @param seed Random seed for determining indices (if not specified)
#' @param max.iter Maximum number of alternating NNLS solutions for H and W, integer (default 1000)
#' @param trace An integer specifying a multiple of iterations at which MSE error should be calculated and checked for convergence. To check error every iteration, specify 1. To avoid checking error at all, specify trace > max.iter (default is 5, and is generally an efficient and effective value)
#' @return A list of NMF models and angle between models for each rank. If multiple ranks are specified then summary statistics will also be calculated.
nnmf.match <- function(A, ind = NULL, byrow = FALSE, k = NULL, max.iter = 1000, rel.tol = 1e-3, n.threads = 0, verbose = 2, trace = 5, seed = NULL) {
  require(qlcMatrix)
  require(RcppHungarian)
  if (class(A)[1] == "dgCMatrix") A <- as.matrix(A)
  if (n.threads < 0) stop("Specify 0 or a positive integer for n.threads")
  if (is.null(k)) stop("Specify a positive integer value or array of integer values for k")
  if (rel.tol > 0.1) warning("rel.tol is greater than 0.1, results may be unstable")
  if (trace < 1) stop("trace must be a positive integer")
  if (is.logical(verbose)) verbose <- as.integer(verbose)
  inner.rel.tol <- 1e-6
  inner.max.iter <- 100
  if (!is.null(seed)) set.seed(seed)
  if (length(ind) == nrow(A)) byrow <- TRUE
  if (length(ind) == ncol(A)) byrow <- FALSE
  if (is.null(ind) && byrow == FALSE) ind <- sample(1:ncol(A), floor(ncol(A) / 2))
  if (is.null(ind) && byrow == TRUE) ind <- sample(1:nrow(A), floor(nrow(A) / 2))
  if (byrow == TRUE) {
    A1 <- A[ind, ]
    A2 <- A[-ind, ]
  }
  if (byrow == FALSE) {
    A1 <- A[, ind]
    A2 <- A[, - ind]
  }
  summary <- list()
  for (i in 1:length(k)) {
    if (verbose > 0) cat("\n\n Running NMF for k =", k[i])

    if (verbose > 0) cat("\n\n Factorizing first subset")
    res1 <- c_nnmf(A1, as.integer(k[i]), as.integer(max.iter), as.double(rel.tol), as.integer(n.threads), as.integer(verbose) - 1, as.integer(inner.max.iter), as.double(inner.rel.tol), as.integer(trace))

    if (verbose > 0) cat("\n\n Factorizing second subset")
    res2 <- c_nnmf(A2, as.integer(k[i]), as.integer(max.iter), as.double(rel.tol), as.integer(n.threads), as.integer(verbose) - 1, as.integer(inner.max.iter), as.double(inner.rel.tol), as.integer(trace))

    if (verbose > 0) cat("\nMatching factors and calculating mean standard error")

    ifelse(byrow == FALSE,
            cosDists <- 1 - as.matrix(cosSparse(res1$W, res2$W)),
            cosDists <- 1 - as.matrix(cosSparse(t(res1$H), t(res2$H)))
        )
    matched <- HungarianSolver(cosDists)$pairs
    colnames(matched) <- c("model1", "model2")
    angles <- apply(matched, 1, function(x) cosDists[x[1], x[2]])
    summary[[i]] <- list("k" = k[i], "angles" = angles, "mean.angle" = mean(angles), "model1" = res1, "model2" = res2, "matched.factors" = matched)
  }
  if (length(summary) == 1) {
    return(summary[[1]])
  } else {
    best.k <- summary[[which.min(unlist(lapply(summary, function(x) x["mean.angle"])))]]$k
    min.mean.angle <- min(unlist(lapply(summary, function(x) x["mean.angle"])))
    mean.angles <- data.frame(cbind("k" = unlist(lapply(summary, function(x) x[["k"]])), "mean.angle" = unlist(lapply(summary, function(x) x["mean.angle"]))))
    summary[["best.k"]] <- best.k
    summary[["min.mean.angle"]] <- min.mean.angle
    summary[["mean.angles"]] <- mean.angles
    return(summary)
  }
}