#' @title Find samples dominated by specific taxa
#' @description Finding the samples dominated by user provided taxa in a phyloseq object.
#' This is useful especially if user suspects a taxa to be contaminant and wishes to identify
#' which samples are dominated by the contaminant taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxa this should match the rownames of otu_table(x)
#' @param relative Logical. If TRUE will transform input to relative abundance. Default=FALSE
#' @return A character with sample ids.
#' @examples
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' p0.f <- aggregate_taxa(p0, "Genus")
#' bac_dom <- find_samples_taxa(p0.f, taxa = "g__Bacteroides")
#' # get samples dominated by g__Bacteroides
#' ps.sub <- prune_samples(sample_names(p0.f) %in% bac_dom, p0.f)
#' @export
#' @keywords utilities
#'
find_samples_taxa <- function(x, taxa = NULL, relative = FALSE) {
  abund <- y <- NULL
  if (is.null(taxa)) {
    stop("Please specific name of taxa")
  }

  if (relative == TRUE) {
    abund <- abundances(x, transform = "compositional")
  }

  abund <- abundances(x)

  if (taxa %in% taxa_names(x)) {
    y <- apply(abund, 2, function(x) {
      order(x, decreasing = TRUE)[1] == grep(taxa,
        rownames(abund),
        perl = TRUE
      )
    })

    return(colnames(abund[, y]))
  } else {
    stop(paste0(taxa, " does not match any rownames of otu_table!!!"))
  }
}
