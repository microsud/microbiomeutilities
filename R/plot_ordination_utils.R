#' @title Plot species loading with ordinations
#' @description This function extends the plot_ordination function of \code{\link{phyloseq}}to highlight the top taxa loadings
#' on the species ordination.
#' @details This function is useful for visualizing specifc taxa that could be important in explaining variations in ordinations.
#' @param x \code{\link{phyloseq-class}} object
#' @param ordiObject Output of ordinate from package phyloseq. Only NMDS/CCA and Bray supported.
#' @param plot.arrow If arrow should be plotted for species either TRUE or FALSE.
#' @param scale.arrow If arrow is plotted a constant to multiply axis values for clearing visualisations.
#' @param top.taxa Top varying taxa to plot, default is 5.
#' @param color.opt Variable of interest from metadata.
#' @import ggplot2
#' @import ggrepel
#' @import phyloseq
#' @import microbiome
#' @return plot
#' @export
#' @examples \dontrun{
#' 
#'     library(microbiomeutilities)
#'     library(RColorBrewer)
#'     data("zackular2014")
#'     ps1 <- zackular2014
#'     ps2 <- tax_glom(ps1, "Genus")
#'     ps2f <- format_to_besthit(ps2)     
#'     orddi <- ordinate(ps2f, method = "CCA", distance = "bray")
#'     p <- plot_ordination_utils(ps2f, orddi, 
#'     color="DiseaseState", plot.arrow = TRUE, 
#'     scale.arrow = 1.3, top.taxa = 5)
#'     print(p)
#' }
#' @keywords utilities
#'
plot_ordination_utils <- function(x,
                                  ordiObject,
                                  color.opt = NULL,
                                  plot.arrow = TRUE,
                                  scale.arrow = NULL,
                                  top.taxa = 5)
{
  p.base <-
    y.axis <-
    x.axis <-
    pdf.sam <-
    pdf.tax <-
    pdf.tax2 <-
    id.type <-
    best_hit <-
    select.top <- dif.tax.ord <- pdf.tax3 <- plot.ord.load <- NULL
  p.base <- plot_ordination(x, ordiObject, color = color.opt,
                            type = "split")
  y.axis <- p.base$labels$y[1]
  x.axis <- p.base$labels$x[1]
  pdf.base <- p.base$data
  # unique(pdf$id.type)
  pdf.sam <- subset(pdf.base, id.type == "Samples")
  pdf.tax <- subset(pdf.base, id.type == "Taxa")
  rownames(pdf.tax) <- pdf.tax$best_hit
  pdf.tax2 <- pdf.tax[, 1:2]
  select.top <- min(top.taxa, dim(pdf.tax2)[1])
  diff.taxa.ord <- names(sort(apply(abs(pdf.tax2), 1, max),
                              decreasing = T))[1:select.top]
  pdf.tax3 <- subset(pdf.tax, best_hit %in% diff.taxa.ord)
  plot.ord.load <-
    ggplot() + geom_point(data = pdf.sam,
                          aes(pdf.sam[,
                                      1], pdf.sam[, 2], color = pdf.sam[, color.opt]))
  plot.ord.load <- plot.ord.load + geom_text(data = pdf.tax3,
                                             aes(scale.arrow * pdf.tax3[, 1], scale.arrow * pdf.tax3[,
                                                                                                     2], label = best_hit)) + ggrepel::geom_text_repel()
  if (plot.arrow == TRUE)
  {
    plot.ord.load <- plot.ord.load + geom_segment(
      data = pdf.tax3,
      aes(
        x = 0,
        xend = scale.arrow * pdf.tax3[, 1],
        y = 0,
        yend = scale.arrow * pdf.tax3[, 2]
      ),
      arrow = arrow(length = unit(0.4,
                                  "cm")),
      color = "grey50"
    )
    plot.ord.load <- plot.ord.load + theme_bw() + xlab(x.axis) +
      ylab(y.axis)
    return(plot.ord.load)
  } else
  {
    plot.ord.load <- plot.ord.load + theme_bw() + xlab(x.axis) +
      ylab(y.axis)
    return(plot.ord.load)
  }
}
