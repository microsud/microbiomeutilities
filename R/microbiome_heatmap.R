#' @title Heatmap using phyloseq and pheatmap
#' @description Plot heatmap using \code{\link{phyloseq-class}} object as input.
#' @param phyobj \code{\link{phyloseq-class}} object.
#' @param subset_otu either NA or number of Top OTUs to use for plotting.
#' @param transformation either 'log10', 'clr','Z', 'compositional', or NA
#' @param VariableA main variable of Interest.
#' @param heatcolors is the option for colors in \code{\link{pheatmap}}. Default is to use viridis
#'        inferno
#' @param ... Arguments to be passed (for \code{\link{pheatmap}} pheatmap).
#' @return A \code{\link{pheatmap}} plot object
#' @export
#' @author Sudarshan A. Shetty (sudarshanshetty9@gmail.com)
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data(DynamicIBD)
#'     ps0 <- DynamicIBD
#'     colnames(tax_table(ps)) <- c("Domain", "Phylum", "Order", "Class", "Family", "Rank6" "Rank7")
#'
#'     p <- microbiome_heatmap(ps, subset_otu = 10,
#'     VariableA = "ibd_subtype",
#'     heatcolors = NA)
#'           }
#' @keywords utilities
microbiome_heatmap <- function(phyobj, subset_otu, transformation,
                               VariableA, heatcolors, ...)
{


  topOTU <- phyobj1 <- phyobj2 <- otu.mat <- meta.tab <- tax.tab <- taxDF <- tax.lev <- select.meta <- color.heatmap <- NULL
    if (!is.na(subset_otu))
    {

      message(paste0("Top ", subset_otu, " OTUs selected",
                     sep = " "))
      topOTU <- top_taxa(phyobj, n = subset_otu)
      phyobj1 <- prune_taxa(topOTU, phyobj)

    } else
    {

      message("subset otu is NA hence will use input values to plot")

      phyobj1 <- phyobj

    }

  message("log10, if zeros in data then log10(1+x) will be used")

  if (transformation == "log10"){
    phyobj2 <- transform(phyobj1, "log10")
  } else if (transformation == "compositional"){
    phyobj2 <- transform(phyobj1, "compositional")
  } else if(transformation == "Z-OTU"){
    phyobj2 <- transform(phyobj1, "Z")
  } else if(transformation == "clr"){
    phyobj2 <- transform(phyobj1, "clr")
  } else if (!is.na(transformation)){
    phyobj2 <- phyobj1
  }


  # format the taxonomy to incluse unique names
  phyobj2 <- format_phyloseq(phyobj2)

  otu.mat <- abundances(phyobj2)
  meta.tab <- meta(phyobj2)

  # get the taxonomy data
  tax.tab <- tax_table(phyobj2)
  taxDF <- as.data.frame(tax.tab)

  # We will merege all the column into one except the
  # Doamin as all is bacteria in this case
  if (ncol(tax_table(phyobj2)) == 6){
    tax.lev <- tidyr::unite(taxDF, Taxa_level, c("Family", "Genus"), sep = " ", remove = TRUE)
  } else {
    tax.lev <- tidyr::unite(taxDF, Taxa_level, c("Genus",
                                                 "Species"), sep = " ", remove = TRUE)
  }

  tax.lev$Taxa_level <- gsub(pattern = "[a-z]__",
                             replacement = "", tax.lev$Taxa_level)


  # choose which variables of interest to include in
  # the heatmap
  select.meta <- subset(meta.tab, select = c(VariableA))

  rownames(otu.mat) <- as.list(tax.lev$Taxa_level)

  if (is.na(heatcolors))
  {
    color.heatmap <- inferno(10)

  } else
  {

    color.heatmap <- heatcolors

  }
  heatmap <- pheatmap::pheatmap(otu.mat, annotation_col = select.meta,
                     main = "Heatmap", color = color.heatmap, ...)
  return(heatmap)
}
