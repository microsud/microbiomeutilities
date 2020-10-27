#' @title Area plot
#' @description Create an area plot for longitudinal samples 
#' with \code{\link{ggplot2}} package.
#' @param x \code{\link{phyloseq-class}} object.
#' @param xvar Column name to plot on x-axis.
#' @param facet.by Column with variable that has multiple measurements.
#' @param fill.colors brewer.pal(6,"Paired"). Specify colors.
#' @param abund.thres = 0.001 check \code{\link{microbiome}} package aggregate_rare function.
#' @param prev.thres = 0.1 check \code{\link{microbiome}} package aggregate_rare function.
#' @param level Taxonomic level. OTU/ASV level not supported. 
#' @param ncol wrap, specify number of columns.
#' @param nrow wrap, specify number of rows.
#' @return  \code{\link{ggplot}} object.
#' @export 
#' @examples 
#' \dontrun{
#' library(microbiomeutilities)
#' data("hmp2")
#' ps <- hmp2
#' ps.rel <- microbiome::transform(ps, "compositional") 
#' p <- plot_area(ps.rel, xvar="visit_number", 
#'               level = "Phylum",
#'               facet.by = "subject_id",
#'               fill.colors=brewer.pal(6,"Paired"))
#'}           
#' 
#' @keywords visualization
#' 
plot_area <- function(x,
                      xvar=NULL,
                      level = NULL,
                      facet.by = NULL,
                      fill.colors=brewer.pal(6,"Paired"),
                      abund.thres=0.001,
                      prev.thres=0.5,
                      ncol=5,
                      nrow=5){
  
  if(is.null(xvar)){
    stop("xvar cannot be empty")
  }
  if(is.null(facet.by)){
    stop("facet.by cannot be empty")
  }
  x.lev<-xdf<-comp.plt<-facet_var <- NULL
  x.lev <- aggregate_rare(x, 
                          level=level,
                          detection = abund.thres, 
                          prevalence = prev.thres)
  
  xdf <- phy_to_ldf(x.lev, transform.counts = NULL)
  xdf$facet_var <- xdf[,facet.by]
  comp.plt <- ggplot(xdf)
  comp.plt <- comp.plt + 
    geom_area(aes_string(x = xvar, 
                         y = "Abundance", 
                         fill = level)) +
    facet_wrap(~facet_var, 
               scales = "free",
               ncol = ncol,
               nrow = nrow) +
   scale_fill_manual("Taxa", values = fill.colors)
  return(comp.plt + theme_biome_utils() + ylab("Relative abundance"))
  
}
