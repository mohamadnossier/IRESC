
# pre-processing of scRNA-seq datasets (raw counts and normalized counts) ------------------------
# Important !!
  # set counts(sce) <- normcounts(sce) (in normalized counts datasets)

#counts(sce) <- normcounts(sce)

# Packages:
library(FEAST)

#' Pre-processing
#'
#' @name preprocess
#' @param data gene counts matrix
#' @return a pre-processed counts matrix in sce object
#' @examples
#' \dontrun{
#' processed <- preprocess(sce)
#'
#' # extract processed counts
#' processedCounts <- counts(sce)
#' }
#' @import scater
#' @import SingleCellExperiment
#' @importFrom FEAST process_Y
preprocess <- function(datMatrix){
  # Filter ERCC spike-in genes

    if(length(unique((grepl("^ERCC-", rownames(datMatrix)))))==2){
      datMatrix <- datMatrix[-which(grepl("^ERCC-", rownames(datMatrix))),]
  }

  # Filter not expressed genes
  #? filter genes that has total expression sum = 0
  keep_genes <- rowSums(datMatrix > 0) > 0
  #? update the matrix
  datMatrix <- datMatrix[keep_genes, ]

  # Filtering genes with low expression rate based on size factors, as per recommended by FEAST
  filtred <- FEAST::process_Y(datMatrix)

  # return pre-processed
  return(filtred)
} # end function


#' Dimensionality reduction using PCA
#'
#' @name pca
#' @param datMatrix counts matrix
#' @param norm logical parameter to check whether normalized counts or not
#' @return low dimensional components
#' @importFrom stats prcomp
pca <- function(datMatrix, norm = FALSE){

  #? if the counts matrix is normalized (RPKM or FPKM, TPM)
  if(norm == TRUE){
    #? apply log2 +1 transformation
    lg2_tranformation <- log(datMatrix + 1)

    # matrix rotation
    #? transpose the data for (i.e.:genes are the cols and cells are the rows)
    transposed = t(lg2_tranformation)
    #? apply PCA dimentionality reduction
    #? extract x attribute as it contains the low-dimensional components
    pca_results <- stats::prcomp(transposed, center = TRUE, scale. = FALSE)$x
    #? save the low-dimensional components in a dataframe
    pca_results <- data.frame(pca_results)

    #? use the 1st 15 low-dimensional components as recommended by A3
    return(pca_results[, 1:15]) # as recommended by Krzak M et al.
  }

  # matrix rotation

  #? if that the counts matrix is raw
  #? apply PCA directly because this matrix will be normalized by A1 normalization function and this funtion has a log2 + 1 transformation step
  transposed = t(datMatrix)
  pca_results <- stats::prcomp(transposed, center = TRUE, scale. = FALSE)$x
  pca_results <- data.frame(pca_results)


  # take 3% from resulted low-dimensional components
  #? recommended by A2
  perc <- (nrow(pca_results)/100) * 3
  return(pca_results[,1:perc])

} # end function
