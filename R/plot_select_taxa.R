#' @title A boxplot for user specified list of taxa
#' @description User specifed OTUs are plotted.
#' @details Useful for instances where user is interested only in some OTUs. For example OTUs
#'          reported to be significantly diferent.
#'
#' @param x \code{\link{phyloseq-class}} object.
#' @param select.taxa a character list of taxa to be plotted. eg. select.taxa <- c("OTU-370251", "OTU-311173", "OTU-341024").
#' @param variableA Variable of interested to be checked. This will also be used to color the plot.
#' @param palette Any of the RColorBrewer plettes.
#' @param plot.type Three optons c("stripchart", "boxplot", "violin")
#' @param group.order Default is NULL. a list specifing order of x-axis. E.g. c("H","CRC","nonCRC")
#' @return  \code{\link{ggplot}} object. This can be further modified using ggpubr.
#' @export
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiome)
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' p0.f <- format_to_besthit(p0)
#' select.taxa <- c("OTU-d__denovo31:Dorea", "OTU-d__denovo24:Blautia")
#' p <- plot_select_taxa(p0.f, select.taxa, "DiseaseState", "Paired", plot.type = "stripchart")
#' print(p)
#' }
#' @keywords utilities
#'

plot_select_taxa <- function(x, select.taxa, variableA, 
                             palette, plot.type,
                             group.order = NULL) {
  x.rel <- x.prun <- x.df <- p.box <- p.vio <- p.strp <- NULL

  x.rel <- microbiome::transform(x, "compositional")

  x.prun <- prune_taxa(select.taxa, x.rel)

  x.df <- phy_to_ldf(x.prun, transform.counts = NULL)
  
  make_pairs <- function(x) {
    if (is.character(x) == TRUE) {
      # message("is char")
      var.lev <- unique(x)
    } else if (is.factor(x) == TRUE) {
      # message("is fac")
      var.lev <- levels(x)
    }
    # make a pairwise list that we want to compare.
    lev.pairs <- combn(seq_along(var.lev), 2, simplify = FALSE, FUN = function(i) var.lev[i])
    return(lev.pairs)
  }
  
  if (!is.null(group.order)) {
    x.df[, variableA] <- factor(x.df[, variableA],
                              levels = group.order
    )
  }
  
  
  if (plot.type == "boxplot") {
    p <- ggpubr::ggboxplot(x.df, variableA, "Abundance",
      facet.by = "OTUID",
      color = variableA, palette = palette,
      legend = "right", add = "jitter",
      panel.labs.background = list(fill = "white")
    )
  } else if (plot.type == "violin") {
    p <- ggpubr::ggviolin(x.df, variableA, "Abundance",
      facet.by = "OTUID",
      color = variableA, palette = palette,
      legend = "right", add = "jitter",
      panel.labs.background = list(fill = "white")
    )
  } else if (plot.type == "stripchart") {
    p <- ggpubr::ggstripchart(x.df, variableA, "Abundance",
      facet.by = "OTUID",
      color = variableA, palette = palette,
      legend = "right", panel.labs.background = list(fill = "white")
    ) 
  }

  return(p)
}
