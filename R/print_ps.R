#' @title Overview of \code{\link{phyloseq-class}}
#' @description Prints an overview \code{\link{phyloseq-class}}.
#' @param x \code{\link{phyloseq-class}} object
#' @return Prints information about the \code{\link{phyloseq-class}} object.
#' @import utils
#' @export
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @examples
#'     library(microbiomeutilities)
#'     data("biogeogut")
#'     pseq <- biogeogut
#'     print_ps(pseq)
#'
#' @keywords utilities

print_ps <- function(x){

  comp<- NULL

  message(paste0("object is ", class(x)))
  comp <- length(which(colSums(abundances(x)) > 1))
  if (comp == 0) {
    message("Compositional = Yes")

  }

  else {

    message("Compositional = No")
  }

  message(paste0("ntaxa = ", ntaxa(x)))
  message(paste0("nsamples = ", nsamples(x)))
  message(paste0("nsamplesvariables = ", length(sample_variables(x))))
  message(paste0("nranks = ", length(colnames((tax_table(x))))))

}



















