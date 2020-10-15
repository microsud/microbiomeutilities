#' @title Heatmap using \code{\link{phyloseq-class}} and \code{\link{pheatmap}}
#' @description Plot heatmap using \code{\link{phyloseq-class}} object as input.
#' @param x \code{\link{phyloseq-class}} object.
#' @param subset.top either NA or number of Top OTUs to use for plotting.
#' @param transformation either 'log10', 'clr','Z', 'compositional', or NA
#' @param VariableA main variable of Interest.
#' @param heatcolors is the option for colors in \code{\link{pheatmap}}. Default is to use Spectral
#' @param ... Arguments to be passed \code{\link{pheatmap}}.
#' @return A \code{\link{pheatmap}} plot object.
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @author Sudarshan A. Shetty (sudarshanshetty9@gmail.com)
#' @examples
#'
#' library(microbiomeutilities)
#' library(viridis)
#' library(RColorBrewer)
#' data("zackular2014")
#' ps0 <- zackular2014
#'
#' heat.sample <- plot_taxa_heatmap(ps0,
#'   subset.top = 20,
#'   VariableA = "DiseaseState",
#'   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
#'   transformation = "log10"
#' )
#' @keywords visualization
plot_taxa_heatmap <- function(x, subset.top, 
                              transformation,
                              VariableA, 
                              heatcolors = NULL, ...) {
  topOTU <- phyobj1 <- phyobj2 <- otu.mat <- meta.tab <- select.meta <- color.heatmap <- NULL
  
  
  #x <- ps0
  
  x <- suppressWarnings(suppressMessages(format_to_besthit(x)))
  if (!is.null(subset.top)) {
    message(paste0("Top ", subset.top, " OTUs selected",
                   sep = " "
    ))
    topOTU <- top_taxa(x, n = subset.top)
  } else {
    stop("specify a number/value for subset.top")
  }
  
  
  if (transformation == "log10") {
    message("log10, if zeros in data then log10(1+x) will be used")
    phyobj1 <- prune_taxa(topOTU, x)
    message("First top taxa were selected and \nthen abundances tranformed to log10(1+X)")
    
    prev.tx <- prevalence(phyobj1)
    
    phyobj2 <- transform(phyobj1, "log10")
  } else if (transformation == "compositional") {
    phyobjx <- transform(x, "compositional")
    phyobj2 <- prune_taxa(topOTU, phyobjx)
    
    prev.tx <- prevalence(phyobj2)
    
    message("First converted to compositional \n then top taxa were selected")
  } else if (transformation == "Z-OTU") {
    phyobj1 <- prune_taxa(topOTU, x)
    
    prev.tx <- prevalence(phyobj1)
    
    phyobj2 <- transform(phyobj1, "Z")
    message("First top taxa were selected and \nthen abundances tranformed to Z values")
  } else if (transformation == "clr") {
    phyobj1 <- prune_taxa(topOTU, x)
    
    prev.tx <- prevalence(phyobj1)
    
    phyobj2 <- transform(phyobj1, "clr")
    message("First top taxa were selected and \nthen abundances tranformed to clr")
  } else if (!is.null(transformation)) {
    stop("specify a number for transformation, log10, compositional, Z-OTU, clr")
  }
  
  
  # format the taxonomy to incluse unique names
  # phyobj2 <- format_phyloseq(phyobj2)
  
  
  otu.mat <- abundances(phyobj2)
  meta.tab <- meta(phyobj2)
  
  # choose which variables of interest to include in
  # the heatmap
  select.meta <- subset(meta.tab, select = c(VariableA))
  
  
  if (is.null(heatcolors)) {
    color.heatmap <- brewer.pal(6, "Spectral")
  } else {
    color.heatmap <- heatcolors
  }
  newnames <- NULL
  newnames <- lapply(
    rownames(otu.mat),
    function(x) bquote(italic(.(x))))
  row_df <- NULL
  row_df <- as.data.frame(round(prev.tx*100, 2))
  colnames(row_df) <- c("Prevalence")
  
  
  heatmap <- pheatmap::pheatmap(otu.mat,
                                labels_row = as.expression(newnames),
                                annotation_col = select.meta,
                                #annotation_colors = annotation_colors,
                                annotation_row = row_df,
                                color = color.heatmap, ...)
  
  
  return(list("plot"=heatmap, "tax_tab"=otu.mat))
  
}

