#' @title Taxonomic Composition Plot boxplot
#' @description Plot taxon abundance for samples.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxonomic.level Merge the OTUs (for phyloseq object) into a higher taxonomic level. This has to be one from colnames(tax_table(x)).
#' @param top.otu Top number of taxa to plot.
#' @param VariableA Specify main variable of interest. This should be one of the variables in sample_variables(x).
#' @param keep.other TRUE or FALSE. Default is FALSE. This will not plot taxa group as Other
#' @param title title for the plot
#' @param color any of the palette supported by ggpubr/RColorBrewer packages or  user specified as c("red", "blue").
#' @param group.order Default is NULL. a list specifing order of x-axis. E.g. c("H","CRC","nonCRC")
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiomeutilities)
#' library(RColorBrewer)
#' data("zackular2014")
#' ps1 <- zackular2014
#' pn <- plot_taxa_boxplot(ps1,
#'   taxonomic.level = "Phylum",
#'   top.otu = 3, VariableA = "DiseaseState",
#'   title = "Rel plot", color = "Set2"
#' )
#' print(pn)
#' }
#'
#' @keywords utilities

plot_taxa_boxplot <- function(x, taxonomic.level, top.otu, 
                              keep.other = FALSE, VariableA, 
                              title, color = NULL,
                              group.order = NULL) {
  Abundance <- Taxa <- RelAbun <- NULL
  if (!is.null(x@phy_tree)) {
    message("For plotting purpuses the phy_tree will be removed")
    x@phy_tree <- NULL
  }
  else {
    message("The phy_tree slot is empty, easy to make the plot")
  }


  # taxic <- as.data.frame(x@tax_table)
  taxic <- tax_table(x) %>%
    as("matrix") %>%
    as.data.frame()
  # using abundances function from microbiome as sometime otu_table can have taxa_are_rows FALSE in phyloseq
  # otudf <- as.data.frame(abundances(x)
  otudf2 <- abundances(x) %>%
    as("matrix") %>%
    as.data.frame()
  taxic$OTU <- row.names(otudf2)
  taxmat <- as.matrix(taxic)
  new.tax <- tax_table(taxmat)
  tax_table(x) <- new.tax


  # Merge the taxa at a higher taxonomic level

  if (!taxonomic.level == "OTU") {
    x <- aggregate_top_taxa2(x, taxonomic.level, top = top.otu)
  }

  x1 <- transform(x, "compositional")

  x.df0 <- suppressWarnings(suppressMessages(phy_to_ldf(x1, transform.counts = NULL)))
  x.df0$RelAbun <- as.numeric(x.df0$Abundance * 100)
  x.df0$Taxa <- x.df0[, taxonomic.level]

  if (keep.other == FALSE) {
    x.df0 <- subset(x.df0, Taxa != "Other")
  }
  if (!is.null(group.order)) {
    x.df0[, VariableA] <- factor(x.df0[, VariableA],
                              levels = group.order
    )
  }
  
  p <- ggplot(x.df0, aes(
    x = x.df0[, VariableA],
    y = RelAbun,
    fill = x.df0[, VariableA]
  ))

  p <- p + geom_boxplot() +
    geom_jitter(aes(group = x.df0[, VariableA]), alpha = 0.25)

  p <- p + ggtitle(title) + theme_bw() + 
    facet_wrap(~Taxa, scales = "free")

  p <- p + ylab("Relative abundance (%)") + xlab(taxonomic.level) +
    scale_fill_brewer(VariableA,
      palette = color
    )

  return(p + xlab(""))
}
