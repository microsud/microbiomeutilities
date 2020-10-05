#' @title Diversity plot with stats
#' @param x \code{\link{phyloseq-class}} object
#' @param index diversity index. Calculated using microbiome::alpha
#' @param group Grouping variable to compare
#' @param group.colors Colors for plotting groups
#' @param dot.opacity for ggplot alpha to determine opacity for points
#' @param box.opacity for ggplot alpha to determine opacity for box
#' @param violin.opacity for ggplot alpha to determine opacity for violin
#' @param group.order Default is NULL. a list specifing order of x-axis. 
#' E.g. c("H","CRC","nonCRC")
#' @param stats Logical TRUE or FALSE. Calls ggpubr::stat_compare_means. 
#' @param ... params for ggpubr::stat_compare_means
#' @importFrom ggpubr stat_compare_means rotate_x_text
#' @importFrom gghalves geom_half_violin
#' @examples
#' library(microbiomeutilities)
#' library(ggpubr)
#' data("zackular2014")
#' p0 <- zackular2014
#' mycols <- c("brown3", "steelblue", "grey50")
#' p.m <- plot_diversity_stats(p0, group = "DiseaseState", 
#'                             index = "diversity_shannon", 
#'                             group.order = c("H", "CRC", "nonCRC"), 
#'                             group.colors = mycols)
#' print(p.m)
#' @keywords visualization analysis
#' @export
plot_diversity_stats <- function(x, index,
                                 group = NULL,
                                 group.colors = c("brown3", "steelblue"),
                                 dot.opacity = 0.25,
                                 box.opacity = 0.25,
                                 violin.opacity = 0.5,
                                 group.order = NULL, 
                                 stats = TRUE, 
                                 label.format="p.format", 
                                 ...) {
  
  if (stats==TRUE){
    p <- plot_diversity_with_stats(x, index,
                                   group = group,
                                   group.colors = group.colors,
                                   dot.opacity = dot.opacity,
                                   box.opacity = box.opacity,
                                   violin.opacity = violin.opacity,
                                   group.order = group.order,
                                   label.format=label.format)
  } else if (stats==FALSE)
    p <- plot_diversity_without_stats(x, index=index,
                                      group = group,
                                      group.colors = group.colors,
                                      dot.opacity = dot.opacity,
                                      box.opacity = box.opacity,
                                      violin.opacity = violin.opacity,
                                      group.order = group.order,
                                      label.format=label.format)
}
