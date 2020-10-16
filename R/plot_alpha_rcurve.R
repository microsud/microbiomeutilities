#' @title Rarefaction curves for alpha diversity indices
#' @description Calculates alpha diversity idenx at varying sampling units (sequencing depth).
#' @param x \code{\link{phyloseq-class}} object
#' @param index Default:: "observed",
#' @param subsamples Default: c(100,1000, 2000, 3000, 4000, 5000)
#' @param lower.conf Default: 0.025. If type=CI
#' @param upper.conf Default: 0.975.
#' @param group Default: NULL
#' @param linetype.main For ggplot line type for line by group. Default: 1
#' @param line.opacity.main For ggplot alpha to determine opacity for line by group. Default: 0.5 
#' @param linetype.type For ggplot line type for line CI or SD. Default: 2
#' @param line.opacity.type For ggplot line type to determine opacity for line CI or SD. Default: 0.25
#' @param label.min TRUE or FALSE. Default: TRUE
#' @param label.size Label min size 
#' @param label.color Label min color
#' @param type Either CI (confidence interval) or SD (Standard deviation) Default: CI
#' @examples
#' \dontrun{
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' # e.g. to make range of 
#' # subsamples <- seq(0, 5000, by=100)[-1]
#' p <- plot_alpha_rcurve(p0, index="observed", 
#' lower.conf = 0.025, upper.conf = 0.975, 
#' group="DiseaseState") + 
#' scale_color_brewer(palette = "Paired") + 
#' scale_fill_brewer(palette = "Paired")
#' print(p )
#' }
#' @keywords visualization analysis
#' @export
plot_alpha_rcurve <- function(x, 
                              index="observed",
                              subsamples = c(100,1000, 2000, 3000, 4000, 5000),
                              lower.conf = 0.025,
                              upper.conf = 0.975,
                              group = NULL,
                              linetype.main = 1,
                              line.opacity.main = 0.5,
                              linetype.type = 2,
                              line.opacity.type = 0.25,
                              type = "CI",
                              label.min=TRUE,
                              label.size = 3,
                              label.color="grey70"){
  
  depth_df <- ps.rar <- adiv <- adiv_nw <- depth_dfx <- NULL
  mean.measure <-  sd.measure  <- sub_sample <- d <- NULL
  
  
  
  depth_df <- c()
  for (d in subsamples) {
    #ps.rar <- rarefy_even_depth(x, sample.size = d, verbose = FALSE)
    #adiv <- suppressMessages(alpha(ps.rar, index=index))
    #adiv <- as.data.frame(adiv)
    #colnames(adiv) <- index
    #adiv_nw <- cbind(adiv, meta(ps.rar))
    #print(colnames(adiv_nw) )
    #adiv_nw$sub_sample <- d
    #rownames(adiv_nw) < seq(1:nrow(adiv_nw))
    xd <- rarefy_util(x, d, index=index)
    depth_df <- rbind(depth_df,xd)
  }
  
  #unique(depth_df$DiseaseState)
  #unique(depth_df$sub_sample)
  group <- sym(group)
  index <- sym(index)
  
  depth_dfx <- depth_df %>% 
    group_by(sub_sample,!!group) %>% 
    summarise(mean.measure = mean(!!index, na.rm = TRUE),
              sd.measure = sd(!!index, na.rm = TRUE),
              lower.sd = (mean.measure - sd.measure),
              upper.sd = (mean.measure + sd.measure),
              lower.ci = quantile(!!index, lower.conf, na.rm = TRUE),
              upper.ci = quantile(!!index, upper.conf, na.rm = TRUE)) %>% 
    as.data.frame()
  
  
  p <- ggplot(depth_dfx, aes(sub_sample, mean.measure)) + 
    geom_line(aes_string(color=group), 
              linetype=linetype.main, 
              alpha=line.opacity.main) 
  if (label.min==TRUE){
    p <- p + geom_vline(xintercept = min(subsamples), linetype= linetype.type) +
      geom_text(aes(x=min(subsamples), 
                    label= paste0("Min. depth- ", min(subsamples), " reads/sample"), 
                    y=mean(mean.measure)), 
                colour=label.color, angle=90, vjust = 1.2, size=label.size)
  }
   
  
  if (type=="CI"){
    p <- p + geom_ribbon(aes_string(ymin="lower.ci", 
                                    ymax="upper.ci",
                                    fill = group), 
                         linetype=linetype.type, alpha=line.opacity.type)
  } else if (type=="SD"){
    p <- p + geom_errorbar(aes_string(ymin="lower.sd", 
                                      ymax="upper.sd",
                                      color = group), 
                           linetype=linetype.type, alpha=line.opacity.type)
  }
  p <- p +
    theme_biome_utils() +
    ylab(index) +
    xlab("No. of reads")
  
  return(p)
}




