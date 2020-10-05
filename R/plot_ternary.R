#' @title Ternary plot OTU
#' @description Create a ternary plot \code{\link{ggtern}} package.
#' @details Plots the mean relative abundance of taxa in 3 groups being compared.
#' @param x \code{\link{phyloseq-class}} object
#' @param abund.thres = 0.0001 check \code{\link{microbiome}} package core function
#' remove taxa that are dectected at 0.0001 in less than prev.thres of samples
#' @param prev.thres = 0.1 check \code{\link{microbiome}} package core function
#' @param level = "Genus" Taxonomic level. If OTU/ASV level specify="lowest"
#' Does not support phylum level aggregation 
#' @param level.color = "Phylum" Taxonomic level to color
#' @param level.palette color for level.color
#' @param group Grouping variable to compare, for this plot there has to be three 
#' groups in the data 
#' @param dot.size for ggplot alpha to determine size for points
#' @param dot.opacity for ggplot alpha to determine opacity for points
#' @return  \code{\link{ggplot}} object. 
#' @importFrom ggtern ggtern geom_mask
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' library(microbiome)
#' library(microbiomeutilities)
#' library(dplyr)
#' data("zackular2014")
#' p0 <- zackular2014
#' plot_ternary(p0, group="DiseaseState", 
#' abund.thres=0.0001, level= "Genus", prev.thres=0.25)
#' @keywords visualization

plot_ternary <- function(x, 
                         abund.thres=0.0001,
                         prev.thres=0.1,
                         group=NULL,
                         level= "lowest",
                         level.color="Phylum",
                         dot.size = 3,
                         dot.opacity = 0.25,
                         level.palette=brewer.pal(12,"Paired")){
  
  p0.agg <- Sam_rep <- Abundance <- value <- taxa_df <- OTUID <- taxa_mean <- NULL 
  p0.mr.cr <- level_var <- smas <- taxa_tern <- NULL 
  
  if(level=="lowest"){
    
   p <- plot_tern_otu(x, 
                  abund.thres= abund.thres,
                  prev.thres= prev.thres,
                  group=group,
                  level.color=level.color,
                  dot.size = dot.size,
                  dot.opacity = dot.opacity,
                  level.palette=level.palette)
   return(p)
    
  } else {
    
    p0.mr <- transform(x, "compositional")
    tax_table(p0.mr)[,colnames(tax_table(p0.mr))] <- 
      gsub(tax_table(p0.mr)[,colnames(tax_table(p0.mr))],
           pattern="[a-z]__",replacement="")
    
    p0.agg <- aggregate_taxa(p0.mr, level=level)
    p0.mr.cr <- core(p0.agg,detection = abund.thres, 
                     prevalence = prev.thres)
    #DT::datatable(tax_table(p0.mr.cr))
    
    tax_table(p0.mr.cr)[is.na(tax_table(p0.mr.cr)[,level]),level] <- "Other"
    tax_table(p0.mr.cr)[tax_table(p0.mr.cr)[,level]=="",level] <- "Other"
    #tax_table(p0.mr.cr)[tax_table(p0.mr.cr)[,level]=="",level] <- "Other"
    
    p0.mr <- merge_samples(p0.mr.cr, group = group)
    
    # extract taxonomy 
    taxa_df <- tax_table(p0.mr) %>% 
      as("matrix") %>% 
      as.data.frame()
    taxa_df$OTUID <- rownames(taxa_df)
    
    taxa_df[,level.color][which(taxa_df[,level.color] == "")] <- "Other"
    taxa_df$level_var <- taxa_df[,level]
    
    df <- phy_to_ldf(p0.mr, NULL)
    df$level_var <- df[,level]
    # variable goes to Sam_rep column
    taxa_mean <- df %>% 
      group_by(level_var, Sam_rep) %>% 
      summarise(value=mean(Abundance))
    
    taxa_tern <- taxa_mean %>%
      group_by(Sam_rep, level_var) %>%
      #mutate(index=row_number()) %>%
      pivot_wider(names_from=Sam_rep,values_from=value) %>% 
      left_join(taxa_df, by="level_var")
    
    #taxa_tern$size <- (apply(taxa_tern[2:4], 1, mean))*times  
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
  
}


