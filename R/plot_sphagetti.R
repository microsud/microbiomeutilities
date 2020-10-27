#' @title Spaghetti Plots
#' @param x \code{\link{phyloseq-class}} object
#' @param select.taxa a character list of taxa to be plotted. eg. select.taxa <- c("OTU-370251", "OTU-311173", "OTU-341024").
#' @param group Column describing sample that have multiple measurements. For e.g. Subject_Ids
#' which has all subject_ids listed.
#' @param xvar X-axis variable. This can be visit_number for participants in subject_id.
#' @param line.bg.color Line color for background. Default ="#8d99ae",
#' @param bg.opacity Line opacity for background.
#' @param focus.color Line color for focus being plotted.
#' @param focus.line.size Line size for focus being plotted.
#' @param line.size Line size for background.
#' @param ncol if wrap, specify number of columns.
#' @param nrow if wrap, specify number of rows.
#' @param plot.var One of "by_sample":Many sample one taxa or "by_taxa":Many taxa one sample
#' @return  \code{\link{ggplot}} object. This can be further modified using ggpubr.
#' @export
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiomeutilities)
#' data("hmp2")
#' pseq <- hmp2 # Ren
#' pseq.rel <- microbiome::transform(pseq, "compositional")
#' pseq.relF <- format_to_besthit(pseq.rel)
#' # Choose one participant
#' phdf.s <- subset_samples(pseq.relF, subject_id ==
#'   "Participant_1")
#' # Choose top 12 taxa to visualize
#' ntax <- top_taxa(phdf.s, 12)
#' phdf.s <- prune_taxa(ntax, phdf.s)
#'
#' plot_spaghetti(phdf.s,
#'   plot.var = "by_taxa",
#'   select.taxa = ntax,
#'   xvar = "visit_number",
#'   line.bg.color = "#8d99ae",
#'   focus.color = "#555b6e",
#'   ncol = 3,
#'   nrow = 4,
#'   line.size = 0.2
#' )
#'
#' print(p)
#' }
#' @keywords visualization
#'

plot_spaghetti <- function(x,
                           plot.var = "by_taxa",
                           select.taxa = NULL,
                           xvar = NULL,
                           group = NULL,
                           line.bg.color = "#8d99ae",
                           bg.opacity = 0.5,
                           focus.color = "brown3",
                           ncol = NULL,
                           nrow = NULL,
                           focus.line.size = 0.5,
                           line.size = 1) {
  if (plot.var == "by_sample") {
    if (is.null(group)) {
      stop("if plot.var is 'by_sample', specify variable in 'group' argument ")
    }

    plot_spaghetti_sample(
      x = x,
      select.taxa = select.taxa,
      xvar = xvar,
      group = group,
      line.bg.color = line.bg.color,
      bg.opacity = bg.opacity,
      focus.color = focus.color,
      ncol = ncol,
      nrow = nrow,
      focus.line.size = focus.line.size,
      line.size = line.size
    )
  } else if (plot.var == "by_taxa") {
    plot_spaghetti_taxa(
      x = x,
      select.taxa = select.taxa,
      xvar = xvar,
      group = group,
      bg.opacity = bg.opacity,
      line.bg.color = line.bg.color,
      focus.color = focus.color,
      ncol = ncol,
      nrow = nrow,
      focus.line.size = focus.line.size,
      line.size = line.size
    )
  }
}
