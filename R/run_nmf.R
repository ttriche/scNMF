#' Run NMF on cluster centers calculated from a Seurat object ident class, and project cells onto NMF coordinates
#'
#' @param object Seurat object
#' @param k Decomposition rank, integer
#' @param rel.tol Stop criterion, defined as the relative tolerance between two successive iterations: |e2-e1|/avg(e1,e2). Default 1e-3.
#' @param verbose 0 = no tracking, 1 = progress bar relative to max.iter, 2 = iteration info. Default 2.
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores).
#' @param ident Ident to average across (default is the active Idents()). Specify any valid discrete 'object@meta.data$ slot'.
#' @param assay Assay to average across (default is the active assay). Specify any valid assay name.
#' @param slot Slot to average, default is counts. Specify any valid counts slot within the assay.
#' @param return Seurat object with a dimensional reduction NMF slot
RunNMF <- function(object, k = NULL, ident = "active.ident", assay = "default", slot = "counts", verbose = 1, rel.tol = 1e-4, n.threads = 0) {

  A <- average.expression(object, ident = ident, assay = assay, slot = slot, verbose = verbose)
  if (assay == "default") assay <- DefaultAssay(object)

  if (verbose > 0) message("\n\nRunning NMF on cluster centers:")
  nmf.mod <- nnmf(A, rel.tol = rel.tol, n.threads = n.threads, verbose = verbose, k = k)

  # now project cells 1000 at a time
  data <- GetAssayData(object, assay = assay)
  H <- list()
  n.iters <- ceiling(ncol(data) / 1000)
  if (verbose > 0) message("\nProjecting cells onto NMF coordinates")
  if (verbose > 0) pb <- txtProgressBar(char = "=", style = 3, max = n.iters, width = 50)
  for (i in 1:n.iters) {
    if (i != n.iters) {
      data.iter <- data[, ((i - 1) * 1000 + 1):(i * 1000)]
    } else {
      data.iter <- data[, ((i - 1) * 1000 + 1):ncol(data)]
    }
    if (class(data.iter)[1] == "dgCMatrix") data.iter <- as.matrix(data.iter)
    H[[i]] <- nnmf.project(data.iter, W.model = nmf.mod$W, n.threads = n.threads)
    if (verbose > 0) setTxtProgressBar(pb = pb, value = i)
  }

  H <- do.call(cbind, H)

  # create a new dimensional reduction object
  nmf.dr <- CreateDimReducObject(embeddings = t(H), loadings = nmf.mod$W, key = "NMF_", assay = assay)

  # add the dimensional reduction object to Seurat
  object[["nmf"]] <- nmf.dr
  return(object)
}