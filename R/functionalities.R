
#' @name calc_jaccard_index
#' @param sampData dataframe of samples, labels and clusters
#' @param cluster the passed cluster number
#' @param label the passed true label
#' @return the calculated jaccard index
calc_jaccard_index <- function(sampData, cluster, label){

  options(digits = 4)

  # subset samples that are in the given label and cluster
  labelSamps <- subset(sampData, labels == label)$samples
  clusterSamps <- subset(sampData, clusters == cluster)$samples

  intersect <- length(intersect(clusterSamps, labelSamps))
  uni <- length(union(clusterSamps, labelSamps))

  # calculate jaccard index
  jaccIndex <- intersect/uni
  return(jaccIndex)

} # end function


#' @name jaccardMap
#' @description This function calculate jaccard index between each label and cluster
#' @param samples list of biological samples
#' @param labels list of the biological ground truth labels
#' @param clusters list of assigned clusters
#' @return matrix of clusters x labels filled with their corresponding mapping score
#'
#' @import tibble
#' @importFrom methods is
#' @export
jaccardMap <- function(samples, labels, clusters){

  if(is.null(samples)){
    stop("sample list not found!")
  }
  if(is.null(labels)){
    stop("labels list not found!")
  }
  if(is.null(clusters)){
    stop("clusters list not found")
  }
  if(length(labels) != length(clusters)) {
    stop("two lists have different lengths; no match")
  }

  # aggregate samples, labels and clusters in one matrix
  dataSLC <- tibble(samples, labels, clusters)
  colnames(dataSLC) <- c("samples", "labels", "clusters")

  # ground truth labels
  datLabels <- unique(labels)

  # number of generated clusters
  datClusters <- sort(unique(clusters), decreasing = FALSE)

  # construct a mapping frame
  scoreMatrix <- data.frame(matrix(nrow = length(datClusters), ncol = length(datLabels)))
  colnames(scoreMatrix) <- datLabels
  rownames(scoreMatrix) <- datClusters

  message("Calculating Jaccard Ratio ...")
  # start multi-mapping
  for (c in 1:length(datClusters)) {
    for (l in 1:length(datLabels)) {

      jaccIndex <- calc_jaccard_index(sampData= dataSLC,
                                      label= datLabels[l],
                                      cluster= datClusters[c])
      scoreMatrix[c,datLabels[l]] <- jaccIndex

    } # end 2nd loop
  } # end 1st loop

  message("Done! ...")
  return(scoreMatrix)
} # end function


#' Confusion Matrix Metrics
#'
#' @name calc_confusion_matrix
#' @description calculate the confusion matrix for each entered samples by comparing
#' lables against clusters
#' @param sample is the sample name to search by
#' @param labels is the list of samples true labels
#' @param clusters is the list of samples assigned clusters
#' @return whether the entered sample is a TP, FP or FN
calc_confusion_matrix <- function(sample, labels, clusters){

  `%notin%` <- Negate(`%in%`)
  if ((sample %in% labels) & (sample %in% clusters)) {
    return("TP")

  } else if ((sample %notin% labels) & (sample %in% clusters)) {
    return("FP")

  } else if ((sample %in% labels) & (sample %notin% clusters)) {
    return("FN")

  } else {
    return("Not Found")
  }
} # endFunction


#' Multi-class confusion matrix
#'
#' @name multiclass_handling
#' @description obtain the confusion matrix metrics for multi-class data using the
#' mapping results
#'
#' @param samples is a list of the samples names
#' @param labels is a list of samples true labels
#' @param clusters is a list of the samples assigned clusters
#' @param k is the number of clusters
#' @param mapped is a dataframe of the mapping results
#' @return a dataframe of each class (TP, TN, FP, FN)
#' @importFrom dplyr distinct
#' @import tibble
#'
#' @export
handle_multiclass <- function(samples, labels, clusters, maxmatch.mp){

  mapped <- maxmatch.mp$mapping
  sampLabels <- tibble(samples, labels)
  sampClusters <- tibble(samples, clusters)
  calc_metrics <- data.frame(matrix(nrow = 4, ncol = length(mapped$clusters)))
  row.names(calc_metrics) <- c("TP", "TN", "FP", "FN")
  #colnames(calc_metrics) <- mapped$clusters
  colnames(calc_metrics) <- paste0('(',unlist(mapped$clusters), ')', '-', mapped$labels)

  for (i in 1:nrow(mapped)) {

    label <- subset(sampLabels, labels == mapped[i,1])
    cluster <- subset(sampClusters, clusters == mapped[i,2])
    #merge
    merged <- dplyr::distinct(rbind(label["samples"], cluster["samples"]))

    tmpList <- list()
    for (j in 1:nrow(merged)){

      val <- calc_confusion_matrix(merged$samples[j],
                                   labels = label$"samples",
                                   clusters = cluster$"samples")
      tmpList <- append(tmpList, val)
    }  # end inner loop

    `%notin%` <- Negate(`%in%`)
    calc_metrics["TP",i] <- length(grep("TP", tmpList))
    calc_metrics["TN",i] <- nrow(subset(sampLabels, samples %notin% merged$samples))
    calc_metrics["FP",i] <- length(grep("FP", tmpList))
    calc_metrics["FN",i] <- length(grep("FN", tmpList))
  }

  calc_metrics <- calc_metrics[,order(as.numeric(names(calc_metrics)))]
  return(calc_metrics)
} # end function

