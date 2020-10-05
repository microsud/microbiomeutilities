#' @title Taxonomic Composition Plot boxplot
#' @description Plot taxon abundance for samples.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxonomic.level Merge the OTUs (for phyloseq object) into a higher taxonomic level. This has to be one from colnames(tax_table(x)).
#' @param top.otu Top number of taxa to plot.
#' @param group Specify main variable of interest. This should be one of the variables in sample_variables(x).
#' @param keep.other TRUE or FALSE. Default is FALSE. This will not plot taxa group as Other
#' @param title title for the plot
#' @param group.colors Colors for plotting groups
#' @param dot.opacity For ggplot alpha to determine opacity for points
#' @param box.opacity For ggplot alpha to determine opacity for box
#' @param dot.size For ggplot alpha to determine size for points
#' @param add.violin Loical. If half violoin to the added. Default=TRUE
#' @param violin.opacity If add.violin=TRUE, opacity for violin. 
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
#' ps0 <- zackular2014
#' mycols <- c("brown3", "steelblue", "grey50")
#' pn <- plot_taxa_boxplot(ps0,
#'                         taxonomic.level = "Phylum",
#'                         top.otu = 6, 
#'                         group = "DiseaseState",
#'                         title = "Relative abudance plot",
#'                         keep.other = FALSE,
#'                         group.order = c("H","CRC","nonCRC"),
#'                         group.colors = mycols)
#' print(pn + theme_biome_utils())
#' }
#'
#' @keywords visualization

plot_taxa_boxplot <- function(x, taxonomic.level, 
                              top.otu, 
                              keep.other = FALSE, 
                              group, 
                              title, 
                              group.colors = NULL,
                              group.order = NULL,
                              add.violin= TRUE,
                              violin.opacity=0.25,
                              box.opacity=0.25,
                              dot.opacity=0.25,
                              dot.size=2) {
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
  #x.df0$RelAbun <- as.numeric(x.df0$Abundance * 100)
  x.df0$Taxa <- x.df0[, taxonomic.level]

  if (keep.other == FALSE) {
    x.df0 <- subset(x.df0, Taxa != "Other")
  }
  if (!is.null(group.order)) {
    x.df0[, group] <- factor(x.df0[, group],
                              levels = group.order
    )
  }
  
  p <- ggplot(x.df0, aes(
    x = x.df0[, group],
    y = Abundance,
    fill = x.df0[, group]
  ))

  p <- p + geom_boxplot(aes_string(fill=group), 
                        width = 0.2, 
                        outlier.shape = NA, 
                        alpha = box.opacity) +
    geom_jitter(aes_string(group = x.df0[, group], 
                           color=group), 
                alpha = dot.opacity, size=dot.size)
  if(add.violin==TRUE){
    p  <- p + geom_half_violin(position = position_nudge(x = 0.15, y = 0), 
                               alpha = violin.opacity, side = "r")
  }
  p <- p + ggtitle(title) + theme_bw() + 
    facet_wrap(~Taxa, scales = "free")

  p <- p + ylab("Relative abundance (%)") + xlab(taxonomic.level) +
    scale_fill_manual(group,
      values = group.colors
    ) + 
    scale_color_manual(group,
                      values = group.colors
    ) + 
    scale_y_continuous(labels = scales::percent)

  return(p + xlab(""))
}
