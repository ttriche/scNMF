#' Average cluster expression from Seurat object using fast sparse matrix methods
#' 
#' @param object Seurat object
#' @param ident Ident to average across (default is the active Idents()). Specify any valid discrete 'object@meta.data$ slot'.
#' @param assay Assay to average (default is the active assay). Specify any valid assay name.
#' @param slot Slot to average, default is counts. Specify any valid counts slot within the assay
AverageExpression.Sparse <- function(object, ident = "active.ident", assay = "default", slot = "counts", verbose = 1){
  ifelse(assay == "default", data <- GetAssayData(object, slot = slot), data <- GetAssayData(object, assay = assay, slot = slot))
  data <- Matrix(data, sparse = TRUE)   # convert data to a sparse matrix
  ifelse(ident == "active.ident", idents <- as.vector(Idents(object)), idents <- as.vector(object@meta.data[[ident]])) # get the idents class to average over
  ident.names <- unique(idents)
  result <- NULL
  if(verbose == 1) pb <- txtProgressBar(char = '=', style = 3, max = length(ident.names))
  for(i in 1:length(ident.names)){   # loop through all idents, averaging them in data
    m <- Matrix::rowMeans(data[,which(idents == ident.names[i])])
    ifelse(is.null(result), result <- m, result <- cbind(result, m))
    if(verbose == 1) setTxtProgressBar(pb = pb, value = i)
  }
  colnames(result) <- ident.names
  return(Matrix(result, sparse = TRUE))
}