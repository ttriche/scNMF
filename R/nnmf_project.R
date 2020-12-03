#' Calculate an optimal mapping weights (H) for a set of samples to a latent model (W.model)
#' This function will minimize mean((samples - W.model %*% H)^2)
#'
#' @param samples Samples to be projected onto W.model
#' @param W.model NMF factor model (genes x factors) for mapping to cells
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores)
#' @return A list of fit weights (H)
nnmf.project <- function(samples, W.model, n.threads = 0) {
  if (class(samples)[1] == "dgCMatrix") samples <- as.matrix(samples)
  if (class(W.model)[1] == "dgCMatrix") W.model <- as.matrix(W.model)
  if (n.threads < 0) n.threads <- 0
  if (nrow(W.model) != nrow(samples)) stop("Number of rows in W.model are not equal to number of rows in samples matrix")
  H <- c_project(samples, W.model, as.integer(n.threads), as.integer(500), as.double(1e-8))$H
  rownames(H) <- paste0("NMF_", 1:nrow(H))
  colnames(H) <- colnames(samples)
  return(H)
}