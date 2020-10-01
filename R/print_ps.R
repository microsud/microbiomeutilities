#' @title Overview of \code{\link{phyloseq-class}}
#' @description Prints an overview \code{\link{phyloseq-class}}.
#' @param x \code{\link{phyloseq-class}} object
#' @return Prints information about the \code{\link{phyloseq-class}} object.
#' @import utils
#' @export
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @examples
#' library(microbiomeutilities)
#' data("zackular2014")
#' pseq <- zackular2014
#' print_ps(pseq)
#' @keywords utilities

print_ps <- function(x) {
  comp <- NULL
  ave <- round(sum(sample_sums(x)) / nsamples(x), 2)
  message(paste0("01] object is ", class(x)))
  message(paste0("02] ntaxa = ", ntaxa(x)))
  message(paste0("03] nsamples = ", nsamples(x)))
  message(paste0("04] nsamplesvariables = ", length(sample_variables(x))))
  message(paste0("05] nranks = ", length(colnames((tax_table(x))))))
  message(paste0("06] Min. number of reads = ", min(sample_sums(x))))
  message(paste0("07] Max. number of reads = ", max(sample_sums(x))))
  message(paste0("08] Total number of reads = ", sum(sample_sums(x))))
  message(paste0("09] Average number of reads = ", ave))
  message(paste0("10] Median number of reads = ", median(sample_sums(x))))
  message(paste0("11] Sparsity = ", length(which(abundances(x) ==
    0)) / length(abundances(x))))
  message(paste0(
    "12] Number of singletons = ",
    length(taxa_sums(x)[taxa_sums(x) <=
      1])
  ))
  message(paste0(
    "13] % of taxa that are singletons 
  (i.e. exactly one read detected across all samples) = ",
    mean(taxa_sums(x) == 1) * 100
  ))
}
