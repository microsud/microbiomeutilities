#' @title Summarize the percent taxa classification for \code{\link{phyloseq-class}}
#' @description Summarize the percent taxa classification for \code{\link{phyloseq-class}}.
#' @param x \code{\link{phyloseq-class}} object
#' @return Table with information on percent OTUs classified.
#' @export
#' @author Contact: Sudarshan A. Shetty \email{sudarshanshetty9@@gmail.com}
#' @examples
#' \dontrun{
#' library(microbiomeutilities)
#' data("zackular2014")
#' pseq <- zackular2014
#' percent_classified(pseq)
#' }
#'
#' @keywords utilities
percent_classified <- function(x) {
  ta <- ld <- ldna <- lp <- lpna <- lc <- lcna <- lo <- lona <- lf <- lfna <- lg <- lgna <- ls <- lsna <- df1 <- df2 <- NULL

  ta <- as.data.frame.matrix(tax_table(x))

  message("Only patterns such as [g__] or similar is expected. [g__<empty>] or [g__unclassified] not considered\n
          please convert for eg. g__unclassified to uniform [g__] or NAs")

  tax_table(x)[is.na(tax_table(x)[, 1])] <- "k__"
  tax_table(x)[is.na(tax_table(x)[, 2])] <- "p__"
  tax_table(x)[is.na(tax_table(x)[, 3])] <- "c__"
  tax_table(x)[is.na(tax_table(x)[, 4])] <- "o__"
  tax_table(x)[is.na(tax_table(x)[, 5])] <- "f__"
  tax_table(x)[is.na(tax_table(x)[, 6])] <- "g__"

  if (ncol(ta) == 7) {
    tax_table(x)[is.na(tax_table(x)[, 7])] <- "s__"
  }



  ############################### DOMAIN######################

  if (colnames(ta[1]) == "Domain") {
    ld <- length(grep("k__$", ta$Domain, value = TRUE))
    ldna <- length(anyNA(ta$Domain))

    lev1 <- paste0(signif(((nrow(ta)) - (ld + ldna)) / (nrow(ta)) *
      100), " %")
  } else if (colnames(ta[1]) == "Kingdom") {
    ld <- length(grep("k__$", ta$Kingdom, value = TRUE))
    ldna <- length(anyNA(ta$Kingdom))

    lev1 <- paste0(signif(((nrow(ta)) - (ld + ldna)) / (nrow(ta)) *
      100), " %")
  } else {
    stop(paste("First rank name must be either Kingdom or Domain now it is =", rank_names(x)[1]))
  }

  ############################## PHYLUM#######################


  lp <- length(grep("p__$", ta$Phylum, value = TRUE))
  lpna <- length(anyNA(ta$Phylum))

  lev2 <- paste0(signif(((nrow(ta)) - (lp + lpna)) / (nrow(ta)) *
    100), " %")

  ############################ CLASS########################

  lc <- length(grep("c__$", ta$Class, value = TRUE))
  lcna <- length(anyNA(ta$Class))

  lev3 <- paste0(signif(((nrow(ta)) - (lc + lcna)) / (nrow(ta)) *
    100), " %")

  ############################ ORDER#########################

  lo <- length(grep("o__$", ta$Order, value = TRUE))
  lona <- length(anyNA(ta$Order))

  lev4 <- paste0(signif(((nrow(ta)) - (lo + lona)) / (nrow(ta)) *
    100), " %")


  ########################## FAMILY###########################

  lf <- length(grep("f__$", ta$Family, value = TRUE))
  lfna <- length(anyNA(ta$Family))

  lev5 <- paste0(signif(((nrow(ta)) - (lf + lfna)) / (nrow(ta)) *
    100), " %")

  ######################### GENUS###########################

  lg <- length(grep("g__$", ta$Genus, value = TRUE))
  lna <- length(anyNA(ta$Genus))

  lev6 <- paste0(signif(((nrow(ta)) - (lg + lna)) / (nrow(ta)) *
    100), " %")

  ########################## SPECIES###########################

  if (ncol(ta) == 7) {
    ls <- length(grep("s__$", ta$Species, value = TRUE))
    lsna <- length(anyNA(ta$Species))

    lev7 <- paste0(signif(((nrow(ta)) - (ls + lsna)) / (nrow(ta)) *
      100), " %")


    Taxonomic_Levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTUs/ASVs")
    Percent_Classification <- c(lev1, lev2, lev3, lev4, lev5, lev6, lev7, nrow(ta))
    df1 <- as.data.frame(cbind(Taxonomic_Levels, Percent_Classification))
    return(df1)
  } else {
    Taxonomic_Levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTUs/ASVs")
    Percent_Classification <- c(lev1, lev2, lev3, lev4, lev5, lev6, nrow(ta))

    df2 <- as.data.frame(cbind(Taxonomic_Levels, Percent_Classification))

    return(as.data.frame.matrix(df2))
  }


  #####################################################
}
