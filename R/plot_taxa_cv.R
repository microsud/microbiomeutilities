#' @title Coefficient of variations
#' @description Plots CV for OTUs/ASVs.
#' @details Check if there are spurious OTUs/ASVs.
#' @param x \code{\link{phyloseq-class}} object.
#' @param plot.type scatter or hist (histogram)
#' @return A \code{\link{ggplot}} plot object.
#' @import tidyr
#' @import dplyr
#' @import microbiome
#' @import phyloseq
#' @export
#' @examples
#'     library(microbiome)
#'     library(microbiomeutilities)
#'     data("biogeogut")
#'     p0 <- biogeogut
#'     p <- plot_taxa_cv(p0, plot.type = "hist")
#'     print(p)
#'
#' @keywords utilities
#'
plot_taxa_cv <- function(x, plot.type){

  cal_cv <- function(x) abs(sd(x)/mean(x))
  x.mean.rel <- apply(otu_table(x), 1, function(x) mean(x))
  #head(ps.mean.rel)
  x.cvs <- apply(otu_table(x), 1, function(x) cal_cv(x))
  #head(ps.cvs)

  #ps.dfa <- data.frame(ps.mean.rel)
  #ps.dfb <- data.frame(ps.cvs)
  #ps.dfa$CV <- ps.dfb$ps.cvs
  #colnames(ps.dfa)[1] <- "Mean_abun"
  #head(ps.dfa)
  tax.x.df <- as.data.frame(tax_table(x))
  tax.x.df$MeanAbun <- cbind(x.mean.rel)
  tax.x.df$CV <- cbind(x.cvs)
  h <- hist(tax.x.df$MeanAbun)

  #head(tax.ps.df)
  if(plot.type == "scatter")
  cvplot <- ggplot(tax.x.df,
                   aes(MeanAbun, CV),
                   label = NA) + geom_point(aes(color = Phylum)) + scale_x_continuous(breaks= h$breaks) +
    theme_bw() + xlab("Abundance")

  else{

    cvplot <- ggplot(tax.x.df, aes(CV)) + geom_histogram() + theme_bw()
  }


    return(cvplot)
}




