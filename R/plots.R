
# Themes and plots -----------------------------------

# Packages:
library(ggplot2)
library(reshape2)

#a color palette
gg_bxpalette <- c("#D0D0D0", "#969696",
                "#808080", "#5F5F5F", "#4D4D4D")

barpalette2 <- c( "precision"="#E0E0E0",
                  "recall"="#BDBDBD",
                  "specificity"="#9E9E9E",
                  "balanced accuracy"="#757575",
                  "F1-Score"="#616161",
                  "ARI"="#FFFF00", "NMI"="#FF7F00")

#' Melt dataframe
#'
#' @name longFormat
#' @param dfList list of clustering performance dataframes
#' @return re-formated shape suitable for multiple plotting techniques, eg: heatmap, boxplot, etc..
#'
#' @importFrom reshape2 melt
#' @export
longFormat <- function(dfList){

  longform <- data.frame()
  meth <- list()

  for (i in 1:length(methods)) {

    len <- nrow(methods[[i]]) * ncol(methods[[i]])
    longform <- rbind(longform, melt(as.matrix(methods[[i]]), value.name = "value"))
    meth <- append(meth, rep(names(methods)[i], len))
  } # end for

  longform$methods <- unlist(meth)
  colnames(longform) <- c("metric", "clusters", "value", "method")

  return(longform)
} # end function



#' Sankey diagram
#'
#' @name sankey.diagram
#' @param confusioMatrix confusion matrix
#' @param path file path to save sankey diagram in a .html format
#' @import networkD3
#' @import tidyverse
#' @import hrbrthemes
#' @import circlize
#' @import tidyr
sankey.diag <- function(confusionMatrix, path){

  skFormat <- confusionMatrix %>%
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname)

  colnames(skFormat) <- c("source", "target", "value")
  skFormat$target <- paste(skFormat$target, " ", sep="")

  # list all entities in the flow diagram (source and target)
  nodes <- data.frame(name=c(as.character(skFormat$source), as.character(skFormat$target)) %>% unique())

  # with networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  skFormat$IDsource=match(skFormat$source, nodes$name)-1
  skFormat$IDtarget=match(skFormat$target, nodes$name)-1

  # un-comment to activate the color palette
  #ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

  # create the Network (default palette)
  sk <- sankeyNetwork(Links = skFormat, Nodes = nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name",
                      sinksRight=FALSE,  nodeWidth=40, fontSize=13, nodePadding=20)

  saveNetwork(sk, path)
} # end function


#' @importFrom ggplot2 theme element_blank element_text scale_fill_manual labs coord_fixed


#' Construct alluvial flow diagram
#'
#' @name alluvial_flow.plot
#' @param datMatrix jaccard indecies matrix
#' @return alluvial diagram plot
#'
#' @import ggalluvial
#' @export
alluvial_flow.plot <- function(dtaMatrix){

  lngform <- reshape2::melt(as.matrix(dtaMatrix), value.name = "value")
  colnames(lngform) <- c("cluster", "label", "value")

  allu <- ggplot2::ggplot(data = lngform,
                          aes(axis1 = cluster, axis2 = label, y = value)) +

    ggalluvial::geom_alluvium(aes(fill = label),
                              curve_type = "cubic") +

    ggalluvial::geom_stratum(width = 0.4) +
    ggplot2::geom_text(stat = "stratum",
                       aes(label = after_stat(stratum)),
                       family = "sans", fontface = "bold") +

    ggplot2::scale_x_discrete(limits = c("cluster", "label"),
                              expand = c(0.15, 0.05)) +

    ggplot2::scale_fill_viridis_d() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.text = element_text(family = "sans"),
                   legend.title = element_text(face = "bold"))


  return(allu)
} #end function


#' Scatter plot for labels and clusters
#'
#' @name lb_clust_plot
#' @param mtx a dataframe of low dimensional components
#' @param labels cell types or clustering labels
#' @param lgndPos legend position (none by default)
#' @return scatter plot colored by labels
lb_clust_plot <- function(mtx, labels, fgTitle = "", lgndPos = "none"){

  sct_plot <- ggplot2:: ggplot(mtx,
                   aes(x = PC1, y = PC2, color = labels)) +
    ggplot2:: geom_point(size=0.7) +
    # x- and y- axes labels
    ggplot2::labs(x = "PC1",
                  y = "PC2") +
    ggplot2::theme(plot.title = element_text(fgTitle, face="bold", color="black", size=10),
          legend.position = lgndPos,
          legend.title = element_text(size=9, color="black", face = "bold"),

          axis.text.x = element_text(size=8, color="black"),
          axis.text.y = element_text(size=8, color="black"),
          axis.title = element_text(size=9, face="bold", color="black"),

          panel.background = element_rect(fill = "white", colour = "grey80"),
          panel.grid.major = element_line(colour = "grey92", size=0.4),
          panel.grid.minor = element_line(colour = "grey92", size = 0.4)) +

    ggplot2::coord_fixed()

  return(sct_plot)

} # end function


#' Build confusion matrix plot
#'
#' @name confusionMatrix_plot
#' @param datMatrix matrix of size= labels x clusters
#' @return confusion matrix figure
#'
#' @importFrom reshape2 melt
#' @export
confusionMatrix_plot <- function(datMatrix){

  mlt.dat <- reshape2::melt(as.matrix(datMatrix), value.name = "value")

  gg <- ggplot2::ggplot(mlt.dat, aes(x=Var2, y=as.character(Var1), fill=value)) +
    ggplot2::geom_tile() + theme_bw() + # removing legend for `fill`
    ggplot2::geom_text(aes(label=round(value,2)),  colour = "white") + # printing values
    ggplot2::labs(y = 'Predicted Clusters',
                  x = "Actual Labels",
                  title = "Confusion Matrix") +
    ggplot2::theme(axis.title = element_text(family = "sans", size = 10,
                                    colour = "black", face = "bold"),
          axis.text.y = element_text(size=10,color='black', family = "sans"),
          axis.text.x = element_text(size=10, color='black', family = "sans", angle = +90),
          axis.line = element_line(colour = "grey20"),

          legend.title = element_text(size=9, color="black", face = "bold")) +

    ggplot2::coord_equal()

  return(gg)
} #end function

#' Performance heatmap plot
#'
#' @name metrics_heatmap
#' @param melted.df dataframe formulated to build heatmap
#' @return heatmap plot
#' @export
metrics_heatmap <- function(metric.df){

  mlt.df <- melt(as.matrix(metric.df), value.name = "value")
  colnames(mlt.df) <- c("metric", "cluster", "value")

  performance.plot <- ggplot2::ggplot(mlt.df,
                        aes(x = cluster, y = metric, fill = value)) +

    ggplot2::geom_tile(aes(fill = value), colour = "white") +
    ggplot2::geom_text(aes(label = round(value, 2))) + # printing values
    ggplot2::labs(y = 'Metrics',
                  x = "Clusters",
                  title = "Perfromance Metrics") +
    ggplot2::theme(axis.title = element_text(family = "sans", size = 10,
                                    colour = "black", face = "bold"),
          axis.text.y = element_text(size=8,color='black', family = "sans"),
          axis.text.x = element_text(size=8, color='black', family = "sans", angle = +90),
          axis.title.x = element_text(size=11,color='black', family = "sans"),
          axis.title.y = element_text(size=11,color='black', family = "sans"),
          axis.line = element_line(colour = "white"),

          panel.background = element_blank(),
          legend.title = element_text(size=9, color="black", face = "bold")) +
    ggplot2::coord_equal(ratio = 0.4) +
    viridis::scale_fill_viridis(discrete = FALSE)

  return(performance.plot)
} # end function


#' Performance heatmap plot
#'
#' @name performance.heatmap
#' @param melted.df dataframe formulated to build heatmap
#' @param datasetName dataset name to add to title
#' @return heatmap plot
performance.heatmap <- function(melted.df, datasetName){

  hm.plot <- ggplot2::ggplot(melted.df,
                             aes(x = reorder(clusters), y = metric)) +

    ggplot2::geom_tile(aes(fill = value), na.rm = TRUE) +
    ggplot2::geom_text(aes(label = round(value, 2)),
                       size = 3.5, na.rm = TRUE) + # printing values
    ggplot2::labs(y = 'Metrics',
                  x = "Clusters",
                  title = paste0("Clustering Performance | ",  datasetName)) +
    ggplot2::theme(axis.title = element_text(family = "sans", size = 10,
                                             colour = "black", face = "bold"),
                   axis.text.y = element_text(size=10,color='black'),
                   axis.text.x = element_text(size=10, color='black', angle = +30),
                   axis.title.x = element_text(size=11,color='black', family = "sans"),
                   axis.title.y = element_text(size=11,color='black', family = "sans"),
                   axis.line = element_blank(),

                   panel.background = element_rect(fill = "grey"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   legend.title = element_text(size=9, color="black", face = "bold")) +

    ggplot2::coord_fixed(ratio = 0.8) +
    facet_wrap(~method)

  return(hm.plot)
} # end function



#' Visualize clustering performance by specific metric
#'
#' @name bxplot_metric
#' @param lngform melted dataframe of performance metrics
#' @param metric chosen performance metric, default = F1-Score
#' @return box plot visualization of the chosen metric
#' @export
bxplot_metric <- function(lngform, metric = "F1-Score", pal = gg_bxpalette){

  #subset chosen metric
  bxplotMetric <- subset(lngform, metric == metric)

  bxplot <- ggplot2::ggplot(bxplotMetric,
                            aes(x = method, y=value, fill = method)) +
    ggplot2::stat_boxplot(geom = "errorbar",
                          width = 0.25) +
    ggplot2::geom_boxplot(alpha = 0.8,          # Fill transparency
                          colour = "#474747",   # Border color
                          outlier.colour = 1,
                          width = 0.3) +
    ggplot2::geom_hline(yintercept = 0.5,
                        linetype = "longdash", col = 'red') +
    ggplot2::labs(y = '', x = paste0("Clustering Performance by ", metric)) +
    ggplot2::theme(axis.title = element_text(family = "sans", size = (10),
                                             colour = "black", face = "bold"),
                   axis.text.y = element_text(size=10,color='gray29'),
                   axis.text.x = element_text(size=10, color='gray29'),
                   legend.title = element_text(size=9, color="black", face = "bold"),
                   legend.position = "none") +
    ggplot2::scale_fill_manual(values = pal)

return(bxplot)
} # end function


#' Barplot for benchmark
#'
#' @name benchmark_barplot
#' @param dta dataframe of all package metrics alongside ARI and NMI values
#' @param pal colour palette "must be more than 7 colours" (deafult = barpalette2)
#' @return the barplot fig
#' @export
benchmark_barplot <- function(dta,
                              pal=barpalette2){

   br <- ggplot2::ggplot(dta,
                        aes(x=value, y=factor(metrics, levels = metrics),
                            fill=metrics)) +

    ggplot2::geom_bar(stat = "identity", color = "black",
                      width = 0.8, alpha = 0.7) +

    ggplot2::geom_text(aes(label= round(value,2)), vjust=1.5) +
    ggplot2::labs(x = "", y = "Metric") +
    #ggplot2::scale_x_continuous(expand = expansion(mult = c(0,1))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = element_text(size = 8, color = "black", family = "sans",
                                              angle = 90, vjust=0.5, hjust=1),
                   axis.text.y = element_text(size = 8, color = "black",
                                              family = "sans", face = "bold"),
                   axis.title = element_text(size = 9, face = "bold", color = "black"),

                   panel.background = element_blank(),
                   legend.position = "none",
                   aspect.ratio = 2/1) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::coord_flip()

  return(br)
} # end functions

