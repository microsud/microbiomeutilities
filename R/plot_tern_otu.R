#' @title Ternary plot OTU
#' @description Create a ternary plot \code{\link{ggtern}} package.
#' @details Plots the mean relative abundance of taxa in 3 groups being compared.
#' @param x \code{\link{phyloseq-class}} object.
#' @param abund.thres = 0.0001 check \code{\link{microbiome}} package core function.
#' remove taxa that are dectected at 0.0001 in less than prev.thres of samples.
#' @param prev.thres = 0.1 check \code{\link{microbiome}} package core function.
#' @param level.color = "Phylum" Taxonomic level to color.
#' @param group Grouping variable to compare, for this plot there has to be three 
#' groups in the data 
#' @param dot.size for ggplot alpha to determine size for points
#' @param dot.opacity for ggplot alpha to determine opacity for points
#' @param level.palette color for level.color
#' @return  \code{\link{ggplot}} object. 
#' @examples
#' # library(microbiome)
#' # library(microbiomeutilities)
#' # library(dplyr)
#' # data("zackular2014")
#' # p0 <- zackular2014
#' # plot_tern_otu(p0, group="DiseaseState", 
#' # abund.thres=0.0001, prev.thres=0.25)
#' @keywords visualization
plot_tern_otu <- function(x, 
                          abund.thres=0.0001,
                          prev.thres=0.1,
                          group=NULL,
                          level.color="Phylum",
                          dot.size = 3,
                          dot.opacity = 0.25,
                          level.palette=brewer.pal(12,"Paired")) {
  
  p0.mr <- p0.mr.cr <- df <- taxa_df <- OTUID <- Sam_rep <- taxa_mean <- NULL 
  Abundance <- value <-NULL
  
  if(is.null(group)){
    message("Ternary plot requires 3 groups to compare")
    stop("Please provide group variable to compare")
  }
  
  if(length(unique(meta(x)[,group])) < 3){
    message("Ternary plot requires 3 groups to compare")
    stop(paste0(group, " variable has less that 3 groups to compare"))
  }
  
  if(length(unique(meta(x)[,group])) > 3){
    message("Ternary plot requires 3 groups to compare")
    stop(paste0(group, " variable has more than 3 groups to compare"))
  }
  
  
  p0.mr <- transform(x, "compositional")
  p0.mr.cr <- core(p0.mr,detection = abund.thres, 
                   prevalence = prev.thres)
  
  tax_table(p0.mr.cr)[,colnames(tax_table(p0.mr.cr))] <- 
    gsub(tax_table(p0.mr.cr)[,colnames(tax_table(p0.mr.cr))],
         pattern="[a-z]__",replacement="")
  
  tax_table(p0.mr.cr)[is.na(tax_table(p0.mr.cr)[,level.color]),level.color] <- "Other"
  tax_table(p0.mr.cr)[tax_table(p0.mr.cr)[,level.color]=="",level.color] <- "Other"
  
  p0.mr <- merge_samples(p0.mr.cr, group = group)
  
  df <- phy_to_ldf(p0.mr, NULL)
  
  # extract taxonomy 
  taxa_df <- tax_table(p0.mr) %>% 
    as("matrix") %>% 
    as.data.frame()
  taxa_df$OTUID <- rownames(taxa_df)
  
  taxa_df[,level.color][which(taxa_df[,level.color] == "")] <- "Other"
  
  #taxa_df$level_var <- taxa_df[,level]
  
  taxa_mean <- df %>% 
    group_by(OTUID, Sam_rep) %>% 
    summarise(value=mean(Abundance))
  
  taxa_tern <- smas <- NULL
  
  taxa_tern <- taxa_mean %>%
    group_by(Sam_rep, OTUID) %>%
    #mutate(index=row_number()) %>%
    pivot_wider(names_from=Sam_rep,values_from=value) %>% 
    left_join(taxa_df, by="OTUID")
  
  if(length(level.palette) < length(unique(taxa_df[,level.color]))){
    stop(paste0("Please provide ", length(unique(taxa_df[,level.color])), 
                " color options to match number of ", level.color, "in data"))
  }
  
  smas <- unique(df$Sam_rep)
  
  p <- ggtern(data=taxa_tern, aes_string(x=smas[1], y=smas[2], z=smas[3])) + 
    #theme_rgbw() +
    geom_point(aes_string(color= level.color), 
               alpha=dot.opacity, 
               show.legend=T, 
               size=dot.size) +
    #scale_size(range=c(0, 6)) + 
    geom_mask() + 
    scale_colour_manual(values=level.palette) +
    theme_bw() +
    theme(axis.text=element_blank(), 
          axis.ticks=element_blank()) 
  return(p)
}

  





