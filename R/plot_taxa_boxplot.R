#' @title Taxonomic Composition Plot boxplot
#' @description Plot taxon abundance for samples.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxonomic.level Merge the OTUs (for phyloseq object) into a higher taxonomic level. This has to be one from colnames(tax_table(x)).
#' @param top.otu Top number of taxa to plot.
#' @param VariableA Specify main variable of interest. This should be one of the variables in sample_variables(x).
#' @param title title for the plot
#' @param color any of the palette supported by ggpubr/RColorBrewer packages or  user specified as c("red", "blue").
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiomeutilities)
#'     data("DynamicsIBD")
#'     ps1 <- DynamicsIBD
#'     colnames(tax_table(ps1)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species" )
#'     pn <- plot_taxa_boxplot(ps1,
#'     taxonomic.level = "Phylum",
#'     top.otu = 3, VariableA = "ibd_subtype",
#'     type = "dotplot",title = "rel plot")
#'     print(pn)
#'
#'           }
#' @keywords utilities

plot_taxa_boxplot <- function(x, taxonomic.level, top.otu, VariableA, title, color = NULL){

  if (!is.null(x@phy_tree)){

    message("For plotting purpuses the phy_tree will be removed")
    x@phy_tree = NULL
  }
  else {

    message("The phy_tree slot is empty, easy to make the plot")
  }

  taxic <- as.data.frame(x@tax_table);
  otudf <- as.data.frame(otu_table(x));
  taxic$OTU <- row.names(otudf);
  taxmat <- as.matrix(taxic);
  new.tax <- tax_table(taxmat);
  tax_table(x) <- new.tax;


  # Merge the taxa at a higher taxonomic level

  if (!taxonomic.level == "OTU") {

    x <- aggregate_taxa(x, taxonomic.level, top = top.otu)
  }

  x1 <- transform(x, "compositional")
  x.df0 <- suppressWarnings(suppressMessages(psmelt(x1))) 
  p <- ggboxplot(x.df0, x = taxonomic.level, y = "Abundance", fill = VariableA, palette = color)
  p <- p + ylab("Relative Abundance") + ggtitle(title) + theme(axis.text.x = element_text(face="italic", angle = 90))
  p <- ggpar(p, legend = "right")
  return(p)
}


