#' Average feature expression across clustered samples in a Seurat object using fast sparse matrix methods
#'
#' @param object Seurat object
#' @param ident Ident with sample clustering information (default is the active ident)
#' @param assay Assay to average (default is the active assay)
#' @param slot Slot to average (default is counts)
#' @param verbose Boolean or integer, show progress bar (default is TRUE)
#' @return dense matrix of cluster centers
average.expression <- function(object, ident = "active.ident", assay = "default", slot = "counts", verbose = TRUE) {
  require(Matrix)

  ifelse(assay == "default",
    data <- GetAssayData(object, slot = slot),
    data <- GetAssayData(object, assay = assay, slot = slot))

  # convert data to a sparse matrix if it isn't already
  if (class(data)[1] != "dgCMatrix") data <- Matrix(data, sparse = TRUE)

  # get the idents class over which to average
  ifelse(ident == "active.ident", idents <- as.vector(Idents(object)), idents <- as.vector(object@meta.data[[ident]]))

  # loop through all idents, averaging them in data
  ident.names <- unique(idents)

  if (verbose > 0) message("Averaging expression of ", length(ident.names), " clusters")

  if (verbose > 0) pb <- txtProgressBar(char = "=", style = 3, max = length(ident.names), width = 50)
  m <- list()
  for (i in 1:length(ident.names)) {
    m[[i]] <- Matrix::rowMeans(data[, which(idents == ident.names[i])])
    if (verbose > 0) setTxtProgressBar(pb = pb, value = i)
    }
  result <- do.call(cbind, m)
  colnames(result) <- ident.names
  return(result)
}