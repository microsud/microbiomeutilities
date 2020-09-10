#' @title Heatmap using \code{\link{phyloseq-class}} and \code{\link{pheatmap}}
#' @description Plot heatmap using \code{\link{phyloseq-class}} object as input.
#' @param x \code{\link{phyloseq-class}} object.
#' @param subset.top either NA or number of Top OTUs to use for plotting.
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
#'   heatcolors = brewer.pal(100, "Blues"),
#'   transformation = "log10"
#' )
#' @keywords utilities
plot_taxa_heatmap <- function(x, subset.top, transformation,
                              VariableA, heatcolors = NULL, ...) {
  topOTU <- phyobj1 <- phyobj2 <- otu.mat <- meta.tab <- select.meta <- color.heatmap <- NULL



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
    phyobj2 <- transform(phyobj1, "log10")
  } else if (transformation == "compositional") {
    phyobjx <- transform(x, "compositional")
    phyobj2 <- prune_taxa(topOTU, phyobjx)
    message("First converted to compositional \n then top taxa were selected")
  } else if (transformation == "Z-OTU") {
    phyobj1 <- prune_taxa(topOTU, x)
    phyobj2 <- transform(phyobj1, "Z")
    message("First top taxa were selected and \nthen abundances tranformed to Z values")
  } else if (transformation == "clr") {
    phyobj1 <- prune_taxa(topOTU, x)
    phyobj2 <- transform(phyobj1, "clr")
    message("First top taxa were selected and \nthen abundances tranformed to clr")
  } else if (!is.null(transformation)) {
    stop("specify a number for transformation, log10, compositional, Z-OTU, clr")
  }


  # format the taxonomy to incluse unique names
  # phyobj2 <- format_phyloseq(phyobj2)
  phyobj2 <- suppressWarnings(suppressMessages(format_to_besthit(phyobj2)))

  otu.mat <- abundances(phyobj2)
  meta.tab <- meta(phyobj2)

  # choose which variables of interest to include in
  # the heatmap
  select.meta <- subset(meta.tab, select = c(VariableA))

  # rownames(otu.mat) <- as.list(tax.lev$Taxa_level)

  if (is.null(heatcolors)) {
    color.heatmap <- inferno(10)
  } else {
    color.heatmap <- heatcolors
  }

  heatmap <- pheatmap(otu.mat,
    annotation_col = select.meta,
    main = "Heatmap", color = color.heatmap, ...
  )
  return(heatmap)
}
