#' @title Plasticity
#' @description Calculated difference in microbiota composition for each individual between 
#' two timepoints.
#' @details Using a beta diversity metrics or correlation matrix to identify variability in 
#' microbiota of an individual. The code is slight modification from Grembi et. al. see ref below.
#' This is useful for instance if one wants to quantiify changes in microbiota before and after 
#' a treatment, dietary modulation, antibiotic treatment, etc. The choice of index is important. 
#' For example, Bray-Curtis dissimilarity, the higher values mean higher plasticity/variability. 
#' On the contrary, higer spearman correlation values mean lower plasticity.
#' @param x \code{\link{phyloseq-class}} object
#' @param dist.method Any of the methods supported by phyloseq::distance or correlation method cor()
#' @param participant.col Column name with participant IDs 
#' @importFrom stats cor
#' @return plot
#' @examples
#' library(microbiome)
#' library(microbiomeutilities)
#' library(dplyr)
#' library(ggpubr)
#' data(peerj32)
#' pseq <- peerj32$phyloseq
#' pseq.rel <- microbiome::transform(pseq, "compositional")
#' pl <- plasticity(pseq.rel, participant.col="subject")
#' 
#' @references
#' \itemize{
#' \item{}{Grembi, J.A., Nguyen, L.H., Haggerty, T.D. et al. Gut microbiota plasticity is 
#' correlated with sustained weight loss on a low-carb or low-fat dietary intervention. 
#' Sci Rep 10, 1405 (2020).https://www.nature.com/articles/s41598-020-58000-y
#' }
#' }
#' @export
plasticity <- function(x, dist.method="bray", participant.col){
  
  mdata <- subject <- sampleids <- x1 <- abund_tab <- matrix_dist <-stab_tab <- NULL
  subject_2 <- subject_1 <- stability <- plasticity.value <- S1 <- S2  <- NULL
  
  mdata <- meta(x)
  mdata$subject <- as.factor(mdata[,participant.col])
  mdata$sampleids <- rownames(mdata)
  
   
  mdata <- mdata %>%
    data.frame() %>%
    mutate_at(
      .funs = as.character, 
      .vars = c(sampleids, subject)
    )
  
  x1 <- prune_taxa(taxa_sums(x) > 0 , x)
  #dist.method = "pearson"
  #x1 <- pseq.rel
  if(dist.method %in% c("pearson", "kendall", "spearman")){
    
    abund_tab <- abundances(x1)
    matrix_dist <- cor(abund_tab,abund_tab, method= dist.method,use= "na.or.complete")
    
  } else{
    
    matrix_dist <- phyloseq::distance(x1, dist.method)
    
  }
  
  stab_tab <- reshape2::melt(
    as.matrix(matrix_dist), 
    varnames = c("S1", "S2"),
    value.name = dist.method) %>%
    mutate_if(is.factor, as.character) %>% 
    left_join(mdata, by = c("S1" = "sampleids")) %>%
    left_join(mdata, by = c("S2" = "sampleids"),suffix = c("_1", "_2")) %>% 
    filter(S1 != S2) %>% 
    filter(subject_1 == subject_2) %>%
    mutate(subject = subject_1) %>%
    dplyr::select(-subject_1, -subject_2) %>% 
    distinct(subject, .keep_all = TRUE) 
  return(stab_tab)
}


