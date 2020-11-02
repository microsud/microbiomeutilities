#' @title Taxa abundance summary by group
#' @description Taxa abundance summary by group. Useful for description of microbiome.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Taxonomic level uses microbiome::aggregate_taxa, if NULL with return OTU/ASVs
#' level stats.
#' @param group Provide overview by groups. Default=NULL and will return values for entire
#' dataset, akin to taxa_summary.
#' @param transform Either "compositional" or "counts". Default= compositional 
#' @return A data frames/ grouped tibble
#' @examples
#' \dontrun{
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' get_group_abundances(p0, level = "Phylum", group = "DiseaseState")
#' }
#' @export
#' @keywords utilities
#'
get_group_abundances <- function(x, level, group, transform = "compositional") {
  
  phy_tab <- mean_abundance <- Abundance <- OTUID <- sd_abundance <- NULL
  

  if (!is.null(level)) {
    x <- aggregate_taxa(x, level = level)
  }
  
  if (isFALSE(any(transform ==c("compositional", "counts")))){
    stop("transform can be only compositional or counts ")
  } else if (transform =="compositional"){
    phy_tab <- phy_to_ldf(x, transform.counts = "compositional")
  } else {
    phy_tab <- phy_to_ldf(x, transform.counts = NULL)
  }
  
  if(is.null(group)){
    message("No group specified, return values for entire data")
    phy_tab <- phy_tab %>%
      group_by(OTUID) %>%
      summarise(
        mean_abundance = mean(Abundance, na.rm = T),
        sd_abundance = sd(Abundance, na.rm = T)
      )
    return(phy_tab)
  } else {
    group2 <- sym(group)
    phy_tab <- phy_tab %>%
      group_by(!!group2, OTUID) %>%
      summarise(
        mean_abundance = mean(Abundance, na.rm = T),
        sd_abundance = sd(Abundance, na.rm = T)
      )
    return(phy_tab)
  }
  
  
}
