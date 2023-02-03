
# Clustering methods: SC3 & k-means ------------------------------------------


#' apply SC3 clustering methodology
#'
#' @name apply_SC3
#' @param datMatrix a gene counts matrix after pre-processing
#' @param addProc a logical parameter for additional processing parameter in SC3
#' @return clustering labels
#' @examples
#' eg1:
#' SC3_estimate_results <- apply_SC3(normalized,
#' addProc = TRUE)
#'
#' eg2:
#' SC3_estimate_resutls <- apply_SC3(normalized,
#' addProc = FALSE)
#'
#' eg3:
#' SC3_estimate_resutls <- apply_SC3(normalized,
#' k = 8)
#'
#' @import SC3
#' @import SingleCellExperiment
apply_SC3 <- function(datMatrix, addProc, k = NULL){

  #? construct a SingleCellExperiment object of counts matrix with low-dimensional components, set counts and log counts
  dataObj = SingleCellExperiment(assays = list(counts = as.matrix(datMatrix),
                                              logcounts = log2(as.matrix(datMatrix) + 1)))
  #? add gene ids in feature_symbol attribute
  rowData(dataObj)$feature_symbol <- rownames(dataObj)

  # addProc = FALSE: no additional processing applied
  #? if parameter addProc which stands for additional processing done by SC3 package = FALSE
  #? then assign parameter gene_filter = FALSE, to cancel additional processing
  if(addProc == FALSE){
    dataObj <- SC3::sc3_prepare(dataObj, gene_filter=FALSE)
  }
  # addProc = TRUE: additional processing is applied
  #? if the addProc = TRUE, then assign the gene_filter = TRUE, to perform additional processing
  else if(addProc == TRUE){
    dataObj <- SC3::sc3_prepare(dataObj, gene_filter=TRUE)
  }

  #? if k=NULL, which means that the user chose the k to be estimated by SC3
  #? then call function sc3_estimate_k to estimate the number of clusters k
  if(is.null(k)){

    dataObj <- SC3::sc3_estimate_k(dataObj)
    #? extract the k value
    k <- metadata(dataObj)$sc3$k_estimation
  }

  #? if the k!=NULL, which means that the user chose a k value
  #? then start the SC3 pipeline directly without k estimation
  dataObj <- SC3::sc3_calc_dists(dataObj)
  dataObj <- SC3::sc3_calc_transfs(dataObj)
  dataObj <- SC3::sc3_kmeans(dataObj, ks = k)
  dataObj <- SC3::sc3_calc_consens(dataObj)

  #? extract the clustering results, and format them in a dataframe
  clusters <- data.frame(colData(dataObj), stringsAsFactors = FALSE)[,1]

  return(clusters)
} # end function


#SC3_estimate_resutls <- apply_SC3(normalized, addProc = TRUE)

#' apply kmeans clustering
#'
#' @name apply_kmeans
#' @param dataMatrix an counts matrix, after pre-processing and dimensionality reduction
#' @param k the number of clusters
#' @return clustering labels
apply_kmeans <- function(dataMatrix, k){

  #? to obtain the same results every time
  set.seed(42) # common seed for reproducibility

  #? run R kmeans fucntion, pass counts matrix with the low-dimensional components, k, and additional parameter suggested by A2
  #? extract the resulted clusters
  km_results = as.data.frame(kmeans(dataMatrix, k, nstart = 1000, iter.max = 10000)$cluster)[,1]
  return(km_results)

} # end function

