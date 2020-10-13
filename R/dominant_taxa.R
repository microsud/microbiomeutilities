#' @title Dominant Taxa
#' @description Identify dominant taxa in each sample and give overview.
#' @details Identifies the dominant taxa in each sample and gives an overview of frequency
#' and percent sample that are dominated by each taxon. Can be group wise or overall.
#' @param x \code{\link{phyloseq-class}} object
#' @param level Taxonomic level uses microbiome::aggregate_taxa
#' @param group Provide overview by groups. Default=NULL 
#' @return A list of two data frames/tibbles
#' @examples
#' library(microbiomeutilities)
#' library(dplyr)
#' data("zackular2014")
#' p0 <- zackular2014
#' x.d <- dominant_taxa(p0,level = "Genus", group="DiseaseState")
#' head(x.d$dominant_overview)
#' @export
#' @keywords utilities
dominant_taxa <- function (x, level = NULL, group=NULL) {
  
  sams <- taxs <- out_dat <- meta_dat <- sample_id <-dominant_taxa <- NULL
  rel.freq <- rel.freq.pct <- NULL
  if (!is.null(level)) {
    x <- aggregate_taxa(x, level = level)
  }
  sams <- apply(abundances(x), 2, which.max)
  taxs <- taxa(x)[apply(abundances(x), 2, which.max)]
  out_dat <- data.frame(sample_id = names(sams),
                        dominant_taxa = taxs)
  
  meta_dat <- meta(x) 
  meta_dat$sample_id <- rownames(meta_dat)
  meta_dat <- meta_dat %>% 
    left_join(out_dat)
  
  if(is.null(group)){
    df <-  meta_dat %>% 
      group_by(dominant_taxa) %>% tally() %>% 
      mutate(rel.freq = round(100 * n/sum(n), 1),
             rel.freq.pct = paste0(round(100 * n/sum(n), 0), "%")) %>% 
      arrange(desc(n))
  } else {
    group <- sym(group)
    #dominant_taxa <-"dominant_taxa"
    df <-  meta_dat %>% 
      group_by(!! group, dominant_taxa) %>% tally() %>% 
      mutate(rel.freq = round(100 * n/sum(n), 1),
             rel.freq.pct = paste0(round(100 * n/sum(n), 0), "%")) %>% 
      arrange(desc(n))
  }
  
  return(list(dominant_overview = df, all_data=meta_dat)) 
  
}
