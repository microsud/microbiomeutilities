#' @title A boxplot for user specified list of taxa
#' @description User specifed OTUs are plotted.
#' @details Useful for instances where user is interested only in some OTUs. For example OTUs
#'          reported to be significantly diferent. This can also be used at higher taxonomic
#'          levels, output from phyloseq::tax_glom or microbiome::aggregate_taxa.
#'
#' @param x \code{\link{phyloseq-class}} object.
#' @param select.taxa a character list of taxa to be plotted. eg. select.taxa <- c("OTU-370251", "OTU-311173", "OTU-341024").
#' @param group Grouping variable to compare
#' @param group.colors Colors for plotting groups
#' @param dot.opacity For ggplot alpha to determine opacity for points
#' @param box.opacity For ggplot alpha to determine opacity for box
#' @param add.violin Loical. If half violoin to the added. Default=TRUE
#' @param violin.opacity If add.violin=TRUE, opacity for violin.
#' @param group.order Default is NULL. a list specifing order of x-axis.
#' @param panel.arrange panels "grid" or "wrap" ggplot's facet_XXX
#' @param ncol if wrap, specify number of columns.
#' @param nrow if wrap, specify number of rows.
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
#' mycols <- c("brown3", "steelblue", "grey50")
#' p <- plot_listed_taxa(p0.f, select.taxa,
#'   group = "DiseaseState",
#'   add.violin = TRUE,
#'   group.colors = mycols
#' )
#' print(p)
#' }
#' @keywords visualization
#'

plot_listed_taxa <- function(x,
                             select.taxa,
                             group,
                             group.colors,
                             dot.opacity = 0.25,
                             box.opacity = 0.25,
                             group.order = NULL,
                             add.violin = TRUE,
                             violin.opacity = 0.25,
                             panel.arrange = "grid",
                             ncol = NULL,
                             nrow = NULL) {
  x.rel <- x.prun <- x.df <- p.box <- p.vio <- p.strp <- NULL

  x.rel <- microbiome::transform(x, "compositional")

  x.prun <- prune_taxa(select.taxa, x.rel)

  x.df <- phy_to_ldf(x.prun, transform.counts = NULL)

  if (!is.null(group.order)) {
    x.df[, group] <- factor(x.df[, group],
      levels = group.order
    )
  }

  p <- ggplot(
    data = x.df,
    aes_string(
      x = group, y = "Abundance",
      fill = group
    )
  ) +
    geom_point(aes_string(x = group, color = group),
      position = position_jitter(width = 0.15),
      size = 1, alpha = dot.opacity
    ) +
    geom_boxplot(
      width = 0.2, outlier.shape = NA,
      alpha = box.opacity
    ) +
    guides(fill = FALSE, color = FALSE)
  if (add.violin == TRUE) {
    p <- p + geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      alpha = violin.opacity, side = "r"
    )
  }

  p <- p +
    scale_fill_manual(values = group.colors) +
    scale_colour_manual(values = group.colors)

  if (panel.arrange == "grid") {
    # Make seperate samples based on main varaible
    p <- p + facet_grid(~OTUID, scales = "free")
  } else if (panel.arrange == "wrap") {
    p <- p + facet_wrap(~OTUID,
      scales = "free",
      ncol = ncol,
      nrow = nrow
    )
  }

  p <- p +
    theme(axis.title.x = element_blank())

  return(p + theme_biome_utils())
}
