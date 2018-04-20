#' @title Give taxa summary at specified taxonomic level
#' @description Data frame with mean, max, median standard deviation of relative abundance.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Taxonomic level for which summary is required
#' @return returns a data frame with relative abundance summary.
#' @import utils
#' @export
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @examples \dontrun{
#'     # Example data
#'     library(microbiomeutilities)
#'     data("biogeogut")
#'     p0 <- biogeogut
#'     p0.rel <- microbiome::transform(p0, "compositional")
#'     tx.sum1 <- taxa_summary(p0, "Phylum")
#'
#'     tx.sum2 <- taxa_summary(p0.rel, "Phylum")
#'     }
#'
#' @keywords utilities

taxa_summary <- function(x, level) {

  pobj <- taxdf <- pobj.ag <- otudf <- outputdf
  pobj <- x

  taxdf <- as.data.frame(pobj@tax_table)
  taxdf$OTU <- rownames(tax_table(pobj))
  tax_table(pobj) <- tax_table(as.matrix(taxdf))

  pobj.ag <- microbiome::aggregate_taxa(pobj, level)
  com <- all(sample_sums(pobj.ag) == 1)
    if(com == TRUE){

    message("Data provided is compositional \n will use values directly")

    otudf2 <- as.data.frame(abundances(pobj.ag))

  } else {

    message("Data provided is not compositional \n will first transform")
    otudf2 <- as.data.frame(abundances(pobj.ag, "compositional"))

  }

  output=NULL
  for(j in 1:nrow(otudf2)){
    x2=as.numeric(otudf2[j,])
    mx.rel=max(x2)
    mean.rel=mean(x2)
    med.rel=median(x2)
    Std.dev=sd(x2)

    output=rbind(output,c(row.names(otudf2)[j],mx.rel, mean.rel,med.rel, Std.dev))
  }

  outputdf <- as.data.frame(output)
  colnames(outputdf) <- c("Taxa", "Max.Rel.Ab", "Mean.Rel.Ab", "Median.Rel.Ab", "Std.dev")
  return(outputdf)

}

