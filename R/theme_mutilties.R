#' @title simple theme optimised for visual appeal and clarity
#' @description simple theme optimised for visual appeal and clarity.
#' @import scales
#' @examples \dontrun{
#'
#' library(microbiomeUtilities)
#' data("biogeogut")
#' p0 <- biogeogut

#' select.taxa <- c("OTU-182933:Blautia", "OTU-183089:f__Clostridiaceae")
#' tax_table(p0.f)
#'
#' p0.f <- format_to_besthit(p0)
#' p <- plot_select_taxa(p0.f, select.taxa, "SampleType", "Paired", plot.type = "stripchart")
#' p + theme_microutilities() + scale_colour_Publication()
#' }
#' @keywords utilities


theme_mutilties <- function(){

  theme_bw(base_size=12, base_family="")
  theme(
    text=element_text(family = "Helvetica", face = "plain",
                      color = "black", size = 16,
                      hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                      margin = margin(), debug = FALSE),

    plot.background = element_rect(colour = NA),
    plot.title = element_text(size = rel(1.2), colour = "black"),
    panel.background = element_rect(colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "#1a1a1a", fill=NA, size = 0.25),
    axis.ticks = element_line(colour = "#1a1a1a", size = 0.15),
    axis.text = element_text(size = rel(1.2), colour = "black"),
    axis.title.x = element_text(size = rel(1.0), colour = "black"),
    axis.title.y = element_text(size = rel(1.0), colour = "black", angle=90, vjust =2),
    axis.text.x = element_text(angle = 90),
    axis.title = element_text(size = rel(1)),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size = rel(1.0), colour = "black"),
    legend.key=element_blank(),
    legend.title = element_text(hjust = 0, colour = "black"),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(colour = "black", size = rel(1.0)),
    complete = TRUE)

}




scale_color_mutilties <- function(...){
requireNamespace("scales")
  discrete_scale("colour","Publication",manual_pal(values = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02",
                                                              "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                                                              "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
                                                              "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8",
                                                              "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                                                              "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                                                              "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                                                              "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")), ...)

}


scale_fill_mutilties <- function(...){
  requireNamespace("scales")
  discrete_scale("fill","Publication",manual_pal(values = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02",
                                                            "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                                                            "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
                                                            "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8",
                                                            "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                                                            "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                                                            "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                                                            "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")), ...)

}


