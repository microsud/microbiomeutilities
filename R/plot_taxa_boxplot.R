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
#'     library(RColorBrewer)
#'     data("zackular2014")
#'     ps1 <- zackular2014
#'     pn <- plot_taxa_boxplot(ps1,
#'     taxonomic.level = "Phylum",
#'     top.otu = 3, VariableA = "DiseaseState",
#'     title = "Rel plot", color = "Set2")
#'     print(pn)
#'     }
#'
#' @keywords utilities

plot_taxa_boxplot <- function(x, taxonomic.level, top.otu, VariableA, title, color = NULL){

  Abundance <- NULL
  if (!is.null(x@phy_tree)){

    message("For plotting purpuses the phy_tree will be removed")
    x@phy_tree = NULL
  }
  else {

    message("The phy_tree slot is empty, easy to make the plot")
  }

  taxic <- as.data.frame(x@tax_table);
  # using abundances function from microbiome as sometime otu_table can have taxa_are_rows FALSE in phyloseq
  otudf <- as.data.frame(abundances(x));
  taxic$OTU <- row.names(otudf);
  taxmat <- as.matrix(taxic);
  new.tax <- tax_table(taxmat);
  tax_table(x) <- new.tax;


  # Merge the taxa at a higher taxonomic level

  if (!taxonomic.level == "OTU") {

    x <- aggregate_top_taxa(x, level = taxonomic.level, top = top.otu)

  }

  x1 <- transform(x, "compositional")

  x.df0 <- suppressWarnings(suppressMessages(phy_to_ldf(x1, transform.counts = NULL)))


  p <- ggplot(x.df0, aes(x =x.df0[,taxonomic.level],
                      y=Abundance,
                      fill = x.df0[,VariableA]))

  p <- p + geom_boxplot(position = position_dodge(1)) +
    geom_point(position = position_dodge(width = 0.75),
               aes(group = x.df0[, VariableA]))

  p <- p + ggtitle(title) + theme_bw() +
    theme(axis.text.x = element_text(face ="italic",
                                     angle = 90))

  p <- p + ylab("Relative abundance") + xlab(taxonomic.level) +
    scale_fill_brewer(VariableA,
                      palette = color)

  return(p)

}


