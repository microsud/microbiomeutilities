#' @title Convert \code{\link{phyloseq-class}} object to long data format
#' @description An alternative to psmelt function from \code{\link{phyloseq-class}} object.
#' @param x \code{\link{phyloseq-class}} object
#' @param transform.counts Data transform to be used in plotting
#' (but not in sample/taxon ordering). The options are 'Z-OTU', 'Z-Sample',
#' 'log10' and 'compositional'. See the \code{\link{transform}} function
#' @return A data frame in long format with appropriate transfomation if requested
#' @import tidyr
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @import microbiome
#' @export
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiomeutilities)
#' data("zackular2014")
#' pseq <- zackular2014
#' pseq_df <- phy_to_ldf(pseq, transform.counts = NULL)
#' }
#' @keywords utilities
phy_to_ldf <- function(x, transform.counts) {
  if (is.null(transform.counts)) {
    x <- x
  } else if (transform.counts == "log10") {
    x <- microbiome::transform(x, "log10")
  } else if (transform.counts == "Z-OTU") {
    x <- microbiome::transform(x, "Z", "OTU")
  } else if (transform.counts == "Z-Sample") {
    x <- microbiome::transform(x, "Z", "Sample")
  } else if (transform.counts == "compositional") {
    x <- microbiome::transform(x, "compositional", "OTU")
  } else {
    stop("Please provide appropriate transformation")
  }

  message("An additonal column Sam_rep with sample names is created for reference purpose")
  meta_df <- microbiome::meta(x)
  meta_df$Sam_rep <- rownames(meta_df)
  # tax_df <- data.frame(tax_table(x)) %>%
  #  rownames_to_column("OTUID")
  tax_df <- tax_table(x) %>%
    as("matrix") %>%
    as.data.frame() %>%
    rownames_to_column("OTUID")

  otu_df <- data.frame(abundances(x),
    check.names = FALSE
  ) %>% rownames_to_column("OTUID")
  suppressWarnings(suppressMessages(otu_df %>%
    left_join(tax_df) %>%
    gather_(
      "Sam_rep",
      "Abundance", setdiff(
        colnames(otu_df),
        "OTUID"
      )
    ) %>%
    left_join(meta_df)))
}
