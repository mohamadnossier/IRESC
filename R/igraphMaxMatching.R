
# igraph maximum matching ----------------------------------------

#' @name igraphMaxMatching
#' @param jaccardIndecies dataframe of calculated jaccard index for clusters and labels
#'
#' @return 1:1 mapping of clusters and labels based on maximum matching
#' @importFrom reshape2 melt
#' @importFrom methods is
#' @import igraph
#'
#' @export
igraphMaxMatching <- function(jaccardIndecies){

    # prepare graph pre-requisites
  from <- rownames(jaccardIndecies)
  to <- names(jaccardIndecies)

  melt.df <- reshape2::melt(as.matrix(jaccardIndecies))
  sortedNodesWeight <- melt.df[order(melt.df$Var1),]$value

  message("Build graph using igraph ...")
  #obj <- setClass("JaccMacMapping",
  #         slots=list(mapping="data.frame", cluster="numeric", label="character"))

  # build graph and apply maximum matching using "igraph" package
  bpGraph <- igraph::graph.full.bipartite(length(from), length(to))
  igraph::V(bpGraph)$name <- c(from, to)
  igraph::E(bpGraph)$weight <- sortedNodesWeight

  maxGraph <- igraph::maximum.bipartite.matching(bpGraph)

  mapWeight <- maxGraph$matching_weight
  matched <- data.frame(maxGraph$matching)

  mapObj <- list(mapping = NULL, weight = NULL, clusters = NULL, labels = NULL)

  # check presence of unbalance
  if(TRUE %in% is.na(maxGraph$matching)){

    NAindex <- which(is.na(matched))
    isNA <- rownames(matched)[NAindex]

    for (i in 1:length(isNA)) {

    if(isNA[i] %in% from){
      mapObj$clusters <- append(mapObj$clusters, isNA[i])
      } # end if

    else if(isNA[i] %in% to){
      mapObj$labels <- append(mapObj$labels, isNA[i])
      } # end if

    } # end for
  } # end if

  mapResult <- na.omit(as.data.frame(maxGraph$matching[from]))
  mapResult[,2] <- rownames(mapResult); colnames(mapResult) <- c("labels", "clusters")

  # add results of mapObj list
  mapObj$mapping <- mapResult
  mapObj$weight <- mapWeight

  message("MAPPING DONE ...")
  # return object
  return(mapObj)
} # end function


