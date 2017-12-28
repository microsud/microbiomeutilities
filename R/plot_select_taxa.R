#' @title plot_select_taxa plot selected OTUs of interest
#' @description User specifed OTUs are plotted.
#' @details Useful for instances where user is interested only in some OTUs. For example OTUs 
#'          reported to be significantly diferent.
#'          
#' @param x \code{\link{phyloseq-class}} object.
#' @param select.taxa a character list of taxa to be plotted. eg. select.taxa <- c("OTU-370251", "OTU-311173", "OTU-341024", "OTU-179814").
#' @param variableA Variable of interested to be checked. This will also be used to color the plot.
#' @param palette Any of the RColorBrewer plettes.
#' @param plot.type Three optons c("stripchart", "boxplot", "violin")
#' @return  \code{\link{ggplot}} object. This can be further modified using ggpubr.
#' @import ggpubr 
#' @import microbiome
#' @import RColorBrewer
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     library(microbiomeUtilities)
#'     data("biogeogut")
#'     p0 <- biogeogut
#'     
#'     select.taxa <- c("OTU-370251:Endozoicimonaceae", "OTU-311173:Helicobacteraceae")
#'     
#'     p0.f <- format_to_besthit(p0)
#'     p <- select_taxa_plot(ps2.f, select.taxa, "SampleType", "Paired", plot.type = "stripchart")
#'     
#'     ggpar(p, yscale = "log10", ylab = "Abundance (log10)")
#'     
#'           }
#' @keywords utilities
#' 

plot_select_taxa <- function(x, select.taxa, variableA, palette, plot.type){
  
  x.rel <- x.prun <- x.df <- p.box <- p.vio <- p.strp <- NULL
  
  x.rel <- transform(x, "compositional");
  x.prun <- prune_taxa(select.taxa, x.rel);
  x.df <- phy_to_ldf(x.prun, transform.counts = NULL)
  
  if(plot.type == "boxplot"){
    p <- ggboxplot(x.df, variableA, "Abundance", facet.by =  "OTUID", 
                       color = variableA, palette = palette, 
                       legend = "right", add = "jitter", 
                       panel.labs.background = list(fill = "white"))
    
  } else if (plot.type == "violin"){
    p <- ggviolin(x.df, variableA, "Abundance", facet.by =  "OTUID", 
                       color = variableA, palette = palette, 
                       legend = "right", add = "jitter", 
                       panel.labs.background = list(fill = "white"))
    
  } else if (plot.type == "stripchart"){
    p <- ggstripchart(x.df, variableA, "Abundance", facet.by =  "OTUID", 
                      color = variableA, palette = palette, 
                      legend = "right", panel.labs.background = list(fill = "white"))
    
  } 
  
  return(p)
}


