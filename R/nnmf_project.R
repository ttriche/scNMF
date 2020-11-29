#' Given W and a set of samples, calculate the optimal mapping H such that the objective mean((samples - WH)^2) is minimized
#'
#' @param samples Samples to be projected onto W.model. Dense matrix.
#' @param W.model NMF factor model (genes x factors) for mapping to cells
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores)
#' @return A list of fit weights (H), an imputed counts matrix (A_imputed)
nnmf.project <- function(samples, W.model, n.threads = 0){
    if(n.threads < 0) stop("Specify 0 or a positive integer for n.threads")
    H <- c_project(samples, W.model, as.integer(n.threads), as.integer(500), as.double(1e-8))$H
    rownames(H) <- paste0("NMF_",1:nrow(H))
    colnames(H) <- colnames(samples)
    return(H)
}