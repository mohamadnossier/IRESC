
#' Measure precision
#'
#' @name precision
#' @param TP true positive
#' @param FP false positive
#' @return precision value
precision <- function(TP, FP){
  prec <- TP/(TP+FP); return(prec)
} # endFunction


#' Measure recall
#'
#' @name recall
#' @param TP true positive
#' @param FN false negative
#' @return  recall value
recall <- function(TP, FN){
  rec <- TP/(TP+FN); return(rec)
} #endFunction


#' Measure specificity
#'
#' @name specificity
#' @param TN true negative
#' @param FP false positive
#' @return specificity value
specificity <- function(TN, FP){
  spec <- TN/(TN+FP); return(spec)
} #endFunction


#' Measure balanced accuracy
#'
#' @param recall is the recall value
#' @param specificty is the specificity value
#' @return balanced accuracy value
balancedAccuracy <- function(recall, specificity){
  BA <- (recall + specificity)/2; return(BA)
} #endFunction


#' Measure F1-Score
#'
#' @name F1Score
#' @param pres is the precision value
#' @param rec is the recall value
#' @return F1-Score value
F1Score <- function(prec, rec){
  f1 <- 2 * ((prec * rec) / (prec + rec))
  return(f1)
} #endFunction



#' Measure performance metrics
#'
#' @name measurePerformance
#' @description This function measure the following performace metrics for each cluster and then
#' take the average of each metric along all clusters:
#'    a. Precision
#'    b. Recall
#'    c. Specificity
#'    d. Balanced Accuracy
#'    e. F1-Score
#' @param pn is the confusion matrix
#' @param k is the number of clusters
#' @return a dataframe of all performance metrics calculated
#'
#' @export
measurePerformance <- function(pn){

  options(digits = 3)
  performance <- data.frame(matrix(nrow = 5, ncol = ncol(pn)))
  row.names(performance) <- c("precision", "recall", "specificity",
                              "Balanced Accuracy", "F1-Score")
  colnames(performance) <- names(pn)

  for (i in 1:ncol(pn)) {

    precisionVal <- precision(pn["TP", i], pn["FP", i])
    recallVal <- recall(pn["TP", i], pn["FN", i])
    specificityVal <- specificity(pn["TN", i], pn["FN", i])
    BA <- balancedAccuracy(recallVal, specificityVal)
    f1_score <- F1Score(precisionVal, recallVal)

    performance[,i] <- c(precisionVal, recallVal, specificityVal,
                         BA, f1_score)
  } # end for

  performance[,"Average"] <- rowMeans(performance)
  return(performance)
} # end function


#' Benchmark with ARI and NMI metrics
#'
#' @name bench_process
#' @param labels data ground truth labels
#' @param clusters resulted clusters
#' @return benchmark barplot fig
#' @importFrom aricode ARI, NMI
#' @export
bench_process <- function(labels, clusters, prf.metrics){

  # calculate adjust rand index
  ari <- ARI(labels, clusters)
  #calculate normalized mutual index
  nmi <- NMI(labels, clusters)

  prf <- data.frame(prf.metrics[,"Average"])
  prf[6:7,] <- c(ari, nmi)
  prf[,2] <- c("precision", "recall", "specificity",
               "balanced accuracy", "F1-Score", "ARI", "NMI")
  colnames(prf) <- c("value", "metrics")

  brp <- benchmark_barplot(prf)

  return(brp)
} # end function
