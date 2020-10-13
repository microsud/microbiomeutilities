#' @title Taxa abundance summary by group
#' @description Taxa abundance summary by group. Useful for description of microbiome.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Taxonomic level uses microbiome::aggregate_taxa
#' @param group Provide overview by groups. Default=NULL 
#' @param transform Default compositional 
#' @return A data frames/ grouped tibble
#' @examples
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' get_group_abundances(p0, level = "Phylum", group="DiseaseState")
#' @export
#' @keywords utilities
#' 
get_group_abundances <- function(x, level, group, transform="compositional"){
  
  phy_tab <- mean_abundance <- sd_abundance <- NULL
  group <- sym(group)
  
  if (!is.null(level)) {
    x <- aggregate_taxa(x, level = level)
  }
  
  phy_tab <- phy_to_ldf(x, transform.counts = transform) %>% 
    group_by(!! group, OTUID) %>% 
    summarise(mean_abundance = mean(Abundance, na.rm = T),
              sd_abundance = sd(Abundance, na.rm = T))
  return(phy_tab)
}







