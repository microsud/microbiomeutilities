#' @title Creat a plot for alpha diversities calculated using the \code{\link{microbiome}} package
#' @description Utility plot function for diversity measures calcualted by \code{\link{microbiome}} package.
#' @details Useful for instances where user is interested only in some OTUs. For example OTUs
#'          reported to be significantly diferent.
#'
#' @param x \code{\link{phyloseq-class}} object.
#' @param type Either diversities (Diversity Index) or dominance	(Dominance Index) or evenness	(Evenness Index).
#' @param index.val see global function in \code{\link{microbiome}} package.
#' @param variableA Variable of interested to be checked. This will also be used to color the plot
#' @param palette Any of the RColorBrewer plettes.
#' @param plot.type Three optons c("stripchart", "boxplot", "violin")
#' @return  \code{\link{ggplot}} object. This can be further modified using ggpubr.
#' @import ggpubr
#' @import microbiome
#' @import RColorBrewer
#' @importFrom reshape2 melt
#' @export
#' @examples
#'     library(microbiome)
#'     library(microbiomeutilities)
#'     data("biogeogut")
#'     p0 <- biogeogut
#'     p <- plot_alpha_diversities(p0,
#'     type = "diversities",
#'     index.val = "all",
#'     plot.type = "stripchart",
#'     variableA = "SampleType",
#'     palette = "jco")
#'
#'     print(p)
#'
#' @keywords utilities
#' index.val <- c("shannon, simpson")
plot_alpha_diversities <- function(x, type, index.val = "all", plot.type, variableA, palette){

  x1 <- x
  meta_df <- meta(x1)
  meta_df$sam_rep_nw <- rownames(meta_df)
  if (type == "diversities"){

    adiv <- diversities(x, index = index.val)
    adiv$sam_rep_nw <- rownames(adiv)
  }
  else if (type == "dominance"){
    adiv <- dominance(x, index = index.val)
    adiv$sam_rep_nw <- rownames(adiv)
  }
  else if (type == "evenness"){
    message("This will take some time")
    adiv <- evenness(x, index = index.val)
    adiv$sam_rep_nw <- rownames(adiv)

  } else if (type == "global"){
    adiv <- global(x, index = index.val)
    adiv$sam_rep_nw <- rownames(adiv)
  }
    adiv.nw <- melt(adiv)
    colnames(adiv.nw) <- c("sam_rep_nw","Diversity","div.val")
    meta_df_nw <- melt(meta_df)
    meta_adiv <- merge.data.frame(meta_df_nw, adiv.nw, by = "sam_rep_nw")

    if (plot.type == "boxplot"){
      p <- ggboxplot(meta_adiv, x = variableA, y = "div.val",
                     facet.by = "Diversity", add = "jitter",
                     fill = variableA,
                     palette = palette)

      p2 <- facet(p + theme_bw(), facet.by = "Diversity",
                 short.panel.labs = FALSE,
                 scales = "free",  # Allow long labels in panels
                 panel.labs.background = list(fill = "white"))

    } else if (plot.type == "stripchart"){
      p <- ggstripchart(meta_adiv, x = variableA, y = "div.val",
                     facet.by = "Diversity",
                     color = variableA,
                     palette = palette)

      p2 <- facet(p + theme_bw(), facet.by = "Diversity",
                 short.panel.labs = FALSE,
                 scales = "free",  # Allow long labels in panels
                 panel.labs.background = list(fill = "white"))
    } else if (plot.type == "violin"){
      p <- ggviolin(meta_adiv, x = variableA, y = "div.val",
                        facet.by = "Diversity", add = "jitter",
                        fill = variableA,
                        palette = palette)

      p2 <- facet(p + theme_bw(), facet.by = "Diversity",
                 short.panel.labs = FALSE,
                 scales = "free",  # Allow long labels in panels
                 panel.labs.background = list(fill = "white"))
    }

    return(p2)

  }

