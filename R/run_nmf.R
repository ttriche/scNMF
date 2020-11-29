#' Non-negative matrix factorization on the average expression of clusters in a Seurat object identity class and projection of cells onto NMF coordinates
#'
#' @param object Seurat object
#' @param k Decomposition rank, integer.
#' @param rel.tol Stop criterion, defined as the relative tolerance between two successive iterations: |e2-e1|/avg(e1,e2). Default 1e-3.
#' @param verbose 0 = no tracking, 1 = progress bar relative to max.iter, 2 = iteration info. Default 2.
#' @param n.threads Number of threads/CPUs to use. Default to 0 (all cores).
#' @param ident Ident to average across (default is the active Idents()). Specify any valid discrete 'object@meta.data$ slot'.
#' @param assay Assay to average across (default is the active assay). Specify any valid assay name.
#' @param slot Slot to average, default is counts. Specify any valid counts slot within the assay.
#' @param return Seurat object with a dimensional reduction NMF slot
RunNMF <- function(object, k = NULL, ident = "active.ident", assay = "default", slot = "counts", verbose = 1, rel.tol = 1e-4, n.threads = 0){

    A <- AverageExpression.Sparse(object, ident = ident, assay = assay, slot = slot, verbose = verbose)
    if(assay == "default") assay <- DefaultAssay(object)
    nmf.mod <- nnmf(as.matrix(A), rel.tol = rel.tol, n.threads = n.threads, verbose = verbose, k = k)

    if(verbose > 0) cat("Projecting cells onto NMF coordinates")
    # now project cells onto the learned NMF embeddings
    H <- nnmf.project(as.matrix(GetAssayData(object, assay = assay)), W.model = nmf.mod$W, n.threads = n.threads)

    # there is an issue with creating the DimReduc object

    # nmf.dr <- CreateDimReducObject(embeddings = t(H), loadings.projected = nmf.mod$W, key = "NMF", assay = assay)
    # SetAssayData(object, slot = "reductions", nmf.dr)
    # return(object)
    return(list("W" = nmf.mod$W, "H" = H, "assay" = assay, "key" = "NMF"))
}