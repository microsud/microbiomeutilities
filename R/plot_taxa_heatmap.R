#' @title Heatmap using phyloseq and pheatmap
#' @description Plot heatmap using \code{\link{phyloseq-class}} object as input.
#' @param x \code{\link{phyloseq-class}} object.
#' @param subset_otu either NA or number of Top OTUs to use for plotting.
#' @param transformation either 'log10', 'clr','Z', 'compositional', or NA
#' @param VariableA main variable of Interest.
#' @param heatcolors is the option for colors in \code{\link{pheatmap}}. Default is to use viridis
#'        inferno
#' @param ... Arguments to be passed \code{\link{pheatmap}}.
#' @return A \code{\link{pheatmap}} plot object.
#' @import tidyr 
#' @import microbiome
#' @import phyloseq
#' @import pheatmap
#' @import viridis
#' @export
#' @author Sudarshan A. Shetty (sudarshanshetty9@gmail.com)
#' @examples \dontrun{
#'   # Example data
#'     library(microbiomeutilities)
#'     library(viridis)
#'     data("biogeogut")
#'     ps0 <- biogeogut
#'
#'     p <- plot_taxa_heatmap(ps, taxonomic.level = "Genus", subset.top = 10,
#'     VariableA = "SampleType",
#'     heatcolors = NA)
#'           }
#' @keywords utilities
plot_taxa_heatmap <- function(x, subset.top, transformation,
                              VariableA, heatcolors = NULL, ...)
{
  
  
  topOTU <- phyobj1 <- phyobj2 <- otu.mat <- meta.tab <- select.meta <- color.heatmap <- NULL
  

  
  if (!is.null(subset.top))
  {
    
    message(paste0("Top ", subset.top, " OTUs selected",
                   sep = " "))
    topOTU <- top_taxa(x, n = subset.top)
    phyobj1 <- prune_taxa(topOTU, x)
    
  } else
  {
    
    message("subset.top is not selected hence will use all taxa values to plot")
    
    phyobj1 <- x
    
  }
  
  
  
  if (transformation == "log10"){
    message("log10, if zeros in data then log10(1+x) will be used")
    phyobj2 <- transform(phyobj1, "log10")
  } else if (transformation == "compositional"){
    phyobj2 <- transform(phyobj1, "compositional")
  } else if(transformation == "Z-OTU"){
    phyobj2 <- transform(phyobj1, "Z")
  } else if(transformation == "clr"){
    phyobj2 <- transform(phyobj1, "clr")
  } else if (!is.null(transformation)){
    phyobj2 <- phyobj1
  }
  
  
  # format the taxonomy to incluse unique names
  #phyobj2 <- format_phyloseq(phyobj2)
  phyobj2 <- suppressWarnings(suppressMessages(format_to_besthit(phyobj2))) 
  
  otu.mat <- abundances(phyobj2)
  meta.tab <- meta(phyobj2)
  
  # choose which variables of interest to include in
  # the heatmap
  select.meta <- subset(meta.tab, select = c(VariableA))
  
  #rownames(otu.mat) <- as.list(tax.lev$Taxa_level)
  
  if (is.null(heatcolors))
  {
    color.heatmap <- inferno(10)
    
  } else
  {
    
    color.heatmap <- heatcolors
    
  }
  
  heatmap <- pheatmap(otu.mat, annotation_col = select.meta,
                    main = "Heatmap", color = color.heatmap, 
                    border_color = "white" ,cluster_cols = FALSE,...)
  return(heatmap)
}
