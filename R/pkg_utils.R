#' @title Custom theme for microbiomeutilities pkg
#' @description Opiniated elegant theme.
#' @export
#' @keywords utilities
theme_biome_utils <- function(){
  theme_bw() +
    theme(#axis.text = element_text(size = 16), 
          #axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          #panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          #plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          #plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          #legend.text = element_text(size = 12),          
          #legend.title = element_blank(),                              
          #legend.position = c(0.95, 0.15), 
          strip.background = element_rect(fill="white"),
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

#' @title Pool Taxa
#' @description Creates a list of dataframes at different taxonomic levels.
#' @param x phyloseq object.
#' @examples 
#' # library(phyloseq)
#' # library(microbiome)
#' #data("zackular2014")
#' #sub.sm <- sample(sample_names(zackular2014), 20)
#' # pseq <- prune_samples(sub.sm,zackular2014)
#' # taxa_pooler_mcola(pseq)
#' @keywords utilities
taxa_pooler_mcola <- function(x) {
  
  tax_table(x)[,colnames(tax_table(x))] <- gsub(tax_table(x)[,colnames(tax_table(x))],
                                                pattern="[a-z]__",replacement=" ")
  lev_tax <- pool_tabs <- otu_tab <- OTU <- NULL
  message("Creating a list of tables for ")  
  message(cat(rank_names(x))) 
  lev_tax <- rank_names(x)
  #lev_tax <- lev_tax[-length(lev_tax)]
  pool_tabs <- c()
  for(i in lev_tax){
    tax_table(x)[tax_table(x)[,i]==" ",i] <- NA
    otu_tab <- t(abundances(aggregate_taxa(x, i)))
    pool_tabs[[i]] <- otu_tab
  }
  OTU <- length(lev_tax) + 1
  pool_tabs[[OTU]] <- t(abundances(x))
  return(pool_tabs)
}

#' @title Distribution of taxa
#' @description Plots distribution of taxa.
#' @param x Phyloseq object.
#' @param color.level Taxonomic level to color.
#' @param color.taxa vector of colors specified by user. 
#' Default is brewer.pal(12,"Paired")
#' @return ggplot2 object.
#' @importFrom RColorBrewer brewer.pal
#' @examples  
#' library(microbiomeutilities)
#' data("zackular2014")
#' pseq <- zackular2014
#' p <- taxa_distribution(pseq)
#' p
#' @export
#' @keywords utilities
taxa_distribution <- function(x, color.level="Phylum", 
                              color.taxa=brewer.pal(12,"Paired")) {
  tax_table(x)[,colnames(tax_table(x))] <- gsub(tax_table(x)[,colnames(tax_table(x))],
                                                pattern="[a-z]__",replacement="")
  taxasums <- taxatable <- tax_plot1 <- NULL
  tax_table(x)[tax_table(x)[,color.level]=="",color.level] <- "Unclassified"
  taxasums <- rowSums(abundances(x))
  taxatable <- as.data.frame.matrix(tax_table(x))
  tax_plot1 <- ggplot(taxatable, aes(x = taxasums, 
                                     color = taxatable[,color.level])) + 
    geom_line(size = 0.5, stat = "density") +
    xlab(paste0(color.level, " Counts")) + 
    scale_x_log10() + 
    scale_color_manual(color.level, values = color.taxa) +
    theme_biome_utils()
  return(tax_plot1)
}

