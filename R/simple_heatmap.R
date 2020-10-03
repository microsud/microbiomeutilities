#' @title Simple Heatmap 
#' @description Create a simple heatmap with \code{\link{ggplot2}} package.
#' @details Wrapper converts \code{\link{phyloseq-class}} object to long data frame
#' and generates a heatmap.
#' @param x \code{\link{phyloseq-class}} object.
#' @param group.facet Variable to make facet/panel the plot.
#' @param group.order Default is NULL. a list specifing order of x-axis. 
#' E.g. c("H","CRC","nonCRC")
#' @param abund.thres = 0.01 check \code{\link{microbiome}} package aggregate_rare function.
#' @param prev.thres = 0.1 check \code{\link{microbiome}} package aggregate_rare function.
#' @param level = "Genus" Taxonomic level. OTU/ASV level not supported. Check plot_taxa_heatmap
#' @param scale.color Scale the tiles colors "log10" or "sqrt"
#' @param color.fill User specified color vectors. 
#' @param na.fill Color to fill NAs. e.g. "white"
#' @param remove.other Rare clubbed as Other to be removed. Logical TRUE/FALSE. 
#' @param taxa.arrange Arrange the order of taxa. User can supply a list of vectors. 
#' @param panel.arrange panels "grid" or "wrap" ggplot's facet_XXX
#' @param ncol if wrap, specify number of columns.
#' @param nrow if wrap, specify number of rows.
#' @return  \code{\link{ggplot}} object. 
#' @export
#' @examples
#' library(microbiome)
#' library(microbiomeutilities)
#' library(dplyr)
#' data("zackular2014")
#' p0 <- zackular2014
#' p0.rel <- transform(p0, "compositional")
#' p <- simple_heatmap(p0.rel, group.facet = "DiseaseState",
#'                    group.order = c("H", "CRC", "nonCRC"),
#'                    abund.thres = 0.01,
#'                    prev.thres = 0.1,
#'                    level = "Genus",
#'                    scale.color = "log10",
#'                    na.fill = "white",
#'                    color.fill = NULL,
#'                    taxa.arrange=NULL,
#'                    remove.other=TRUE,
#'                    panel.arrange="wrap",
#'                    ncol=2,
#'                    nrow=2)
#'
#' print(p)
#' @keywords utilities
simple_heatmap <- function(x, group.facet = "DiseaseState",
                           group.order = c("H", "CRC", "nonCRC"),
                           abund.thres = 0.01,
                           prev.thres = 0.1,
                           level = "Genus",
                           scale.color = "log10",
                           na.fill = "white",
                           color.fill = NULL,
                           taxa.arrange=NULL,
                           panel.arrange=NULL,
                           remove.other=TRUE,
                           ncol=NULL,
                           nrow=NULL) {
  
  ps0.gen <- ps_df <- sum.ab <- ord.tx <- NULL
  vec_colors <- p.heat <- Abundance <- NULL
  
  ps0.gen <- aggregate_rare(x, 
                            detection = abund.thres, 
                            prevalence = prev.thres, 
                            level = level)
  
  tax_table(ps0.gen)[,colnames(tax_table(ps0.gen))] <- 
    gsub(tax_table(ps0.gen)[,colnames(tax_table(ps0.gen))],
         pattern="[a-z]__",replacement="")
  
  tax_table(ps0.gen)[is.na(tax_table(ps0.gen)[,level]),level] <- "Other"
  tax_table(ps0.gen)[tax_table(ps0.gen)[,level]=="",level] <- "Other"
  
  ps_df <- phy_to_ldf(ps0.gen, NULL) 
  if (!is.null(group.order)) {
    ps_df[, group.facet] <- factor(ps_df[, group.facet],
                             levels = group.order
    )
  }
  
  ps_df$group_plx <- ps_df[,group.facet]
  ps_df$taxa <- ps_df[,level]
  if (is.null(taxa.arrange)) {
    sum.ab <- ps_df %>% 
      group_by(taxa) %>% 
      summarise(sum.ab = sum(Abundance)) %>% 
      arrange(sum.ab)
    ord.tx <- sum.ab$taxa
  } else {
    ord.tx <- taxa.arrange
  }
  
  ps_df$taxa <- factor(ps_df$taxa, levels=ord.tx)
  
  ## Get colorpalette for colorscale or set default
  if (!is.null(color.fill)) {
    vec_colors <- color.fill
  } else {
    vec_colors <- c("#e63946","#a8dadc","#1d3557")
  }
  
  if(remove.other==TRUE){
    # remove unknown
    ps_df <- subset(ps_df, taxa!= "Other")
  } else {
    ps_df <- ps_df
  }
  
  p.heat <- ggplot(ps_df, 
                   aes_string(x = "Sam_rep", y = "taxa")) + 
    geom_tile(aes(fill = Abundance)) + 
    theme_bw() +
    # Make bacterial names italics
    theme(axis.text.y = element_text(colour = 'black',
                                     size = 10, 
                                     face = 'italic')) 
  if(is.null(panel.arrange)){
    p.heat <- p.heat
  } else if(panel.arrange == "grid"){
    # Make seperate samples based on main varaible 
    p.heat <- p.heat + facet_grid(~group_plx, scales = "free")
  } else if(panel.arrange == "wrap"){
    p.heat <- p.heat + facet_wrap(~group_plx, scales = "free",
                                  ncol=ncol,
                                  nrow=nrow)
  }
  # Make seperate samples based on main varaible 
  p.heat <- p.heat + 
    ylab("Taxa") + 
    #Clean the x-axis
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_blank(),
          strip.background = element_rect(fill="white")) +
    # Clean the facet label box
    theme(legend.key = element_blank(),
          strip.background = element_rect(colour="black", fill="white"))
  
  p.heat <- p.heat + scale_fill_gradientn(colours = vec_colors, 
                                          trans = scale.color, 
                                          na.value = na.fill)
  p.heat
  #scale_fill_distiller("Rel. Abundance (log10 + 1)", palette = "RdYlBu") + rremove("x.text")
  return(p.heat)
}
