#' @title Coefficient of variations
#' @description Plots CV for OTUs/ASVs.
#' @details Check if there are spurious OTUs/ASVs.
#' @param x \code{\link{phyloseq-class}} object.
#' @param plot.type scatter or hist (histogram)
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @importFrom graphics hist
#' @importFrom stats aggregate median sd
#' @importFrom microbiome abundances
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiome)
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' p <- plot_taxa_cv(p0, plot.type = "hist")
#' print(p)
#' }
#'
#' @keywords utilities
#'
plot_taxa_cv <- function(x, plot.type) {
  MeanAbun <- CV <- Phylum <- NULL
  cal_cv <- function(x) abs(sd(x) / mean(x))
  x.mean.rel <- apply(abundances(x), 1, function(x) mean(x))
  # head(ps.mean.rel)
  x.cvs <- apply(abundances(x), 1, function(x) cal_cv(x))
  # head(ps.cvs)

  # tax.x.df <- as.data.frame(tax_table(x))
  tax.x.df <- tax_table(x) %>%
    as("matrix") %>%
    as.data.frame() %>%
    mutate(
      MeanAbun = x.mean.rel,
      CV = x.cvs
    )
  # tax.x.df$MeanAbun <- cbind(x.mean.rel)
  # tax.x.df$CV <- cbind(x.cvs)
  h <- hist(tax.x.df$MeanAbun)

  # head(tax.ps.df)
  if (plot.type == "scatter") {
    cvplot <- ggplot(tax.x.df,
      aes(MeanAbun, CV),
      label = NA
    ) +
      geom_point(aes(color = Phylum)) +
      scale_x_continuous(breaks = h$breaks) +
      theme_bw() +
      xlab("Abundance")
  } else {
    cvplot <- ggplot(tax.x.df, aes(CV)) +
      geom_histogram() +
      theme_bw()
  }

  return(cvplot)
}
