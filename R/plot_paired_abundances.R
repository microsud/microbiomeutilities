#' @title A paired-boxplot for user specified list of taxa
#' @description User specified taxa are plotted.
#' @details Useful for instances where user is interested only in some taxa and thier change after
#'          an intervention. This can also be used at higher taxonomic
#'          levels, output from phyloseq::tax_glom or microbiome::aggregate_taxa.
#'
#' @param x \code{\link{phyloseq-class}} object.
#' @param select.taxa a character list of taxa to be plotted. eg. select.taxa <- c("OTU-370251", "OTU-311173", "OTU-341024").
#' @param group Grouping variable to compare. x axis, eg. before-after, t1-t2.
#' @param group.colors Colors for plotting groups.
#' @param dot.opacity For ggplot alpha to determine opacity for points.
#' @param dot.size For ggplot point size.
#' @param jitter.width Value to avoid over plotting by moving points.
#' @param add.box Logical. If boxplot to the added. Default=TRUE
#' @param box.opacity For ggplot alpha to determine opacity for box.
#' @param add.violin Logical. If half violin to the added. Default=TRUE
#' @param violin.opacity If add.violin=TRUE, opacity for violin.
#' @param group.order Default is NULL. a list specifying order of x-axis.
#' @param ncol If 2 or more taxa to plot, specify number of columns.
#' @param nrow If 2 or more taxa to plot, specify number of rows.
#' @param line Variable to use for lines. E.g. "subject" before-after
#' @param line.down Line Color for change when negative. Decreased abundance.
#' @param line.stable Line Color for no change.
#' @param line.up Line Color for change when positive. Increased abundance.
#' @param line.na.value "grey50" for no/missing observations.
#' @param line.guide "none" to plot guide for line.
#' @param line.opacity Line opacity.
#' @param line.size Size of line to plot.
#'
#' @return  \code{\link{ggplot}} object. This can be further modified using ggpubr.
#' @importFrom tidyr pivot_longer
#' @export
#' @examples
#' \dontrun{
#' library(microbiome)
#' library(microbiomeutilities)
#' library(gghalves)
#' library(tidyr)
#' data(peerj32) # Source: https://peerj.com/articles/32/
#' pseq <- peerj32$phyloseq # Ren
#' pseq.rel <- microbiome::transform(pseq, "compositional")
#' select.taxa <- c("Akkermansia", "Dialister")
#' group.colors <- c("brown3", "steelblue", "grey70")
#' p <- plot_paired_abundances(pseq.rel,
#'   select.taxa = select.taxa,
#'   group = "time",
#'   group.colors = group.colors,
#'   dot.opacity = 0.25,
#'   dot.size = 2,
#'   group.order = NULL,
#'   line = "subject"
#' )
#' p
#' }
#' @keywords visualization

plot_paired_abundances <- function(x,
                                   select.taxa = NULL,
                                   group = NULL,
                                   group.colors = NULL,
                                   dot.opacity = 0.25,
                                   dot.size = 2,
                                   add.box = FALSE,
                                   box.opacity = 0.25,
                                   group.order = NULL,
                                   add.violin = TRUE,
                                   violin.opacity = 0.25,
                                   ncol = NULL,
                                   nrow = NULL,
                                   line = NULL,
                                   line.down = "#7209b7",
                                   line.stable = "#8d99ae",
                                   line.up = "#14213d",
                                   line.na.value = "grey50",
                                   line.guide = "legend",
                                   line.opacity = 0.25,
                                   line.size = 1,
                                   jitter.width = 0) {
  Abundance <- change <- change.order <- change.sign <- linevar <- NULL
  # check arguments
  if (is.na(line) | is.null(line)) {
    stop(" 'line' argument cannot be empty")
  }

  if (is.na(group) | is.null(group)) {
    stop(" 'group' argument cannot be empty")
  }
  #x <- pseq.rel
  xmeta <- meta(x)
  for (i in select.taxa) {
    df.tx <- as.data.frame(abundances(x)[i, ])
    colnames(df.tx) <- i
    xmeta <- cbind(xmeta, df.tx)
  }

  if(length(unique(xmeta[, group])) > 2){
    stop("Only two group comparison e.g. before n after")
  }
  
  # check if factor else convert to factor
  if (!is.factor(xmeta[, group])) {
    xmeta$group <- factor(as.character(xmeta[, group]))
  }

  # convert to wide format
  xmeta_lf <- xmeta %>%
    pivot_longer(
      cols = all_of(select.taxa),
      names_to = "taxa",
      values_to = "Abundance"
    )

  xmeta_lf$linevar <- factor(xmeta_lf[[line]])

  x.grp <- sym(group)

  df2 <- suppressWarnings(xmeta_lf %>%
    arrange(taxa,linevar, !!x.grp) %>%
    group_by(taxa,linevar) %>%
    summarise(change = diff(Abundance)))

  xmeta_lf_2 <- suppressWarnings(xmeta_lf %>%
                     arrange(taxa,linevar, !!x.grp) %>%
                     group_by(taxa,linevar))
  
  xmeta_lf_2 <- xmeta_lf_2 %>% 
    left_join(df2)
  
  #xmeta_lf$change <- df2$change[match(xmeta_lf$linevar, df2$linevar)]
  
  xmeta_lf_2$change.sign <- sign(xmeta_lf_2$change)
  
  if (!is.null(group.order)) {
    xmeta_lf_2[, group] <- factor(xmeta_lf_2[, group],
                                  levels = group.order
    )
  }
  xmeta_lf_2 <- xmeta_lf_2 %>%
    mutate(change.order = ifelse(change.sign == 1, "Up",
      ifelse(change.sign == -1, "Down",
        "Stable"
      )
    ))
  # start plotting
  p <- ggplot(
    data = xmeta_lf_2,
    aes_string(
      x = "group",
      y = "Abundance",
      fill = "group"
    )
  ) +
    geom_point(aes_string(
      x = "group",
      fill = "group"
    ),
    position = position_jitter(width = jitter.width),
    size = dot.size,
    alpha = dot.opacity,
    shape = 21
    ) +
    geom_line(aes(
      group = linevar,
      color = change.order
    ),
    size = line.size,
    alpha = line.opacity
    #lty = line.type
    )

  p <- p + scale_fill_manual(values = group.colors)

  if (add.box == TRUE) {
    p <- p + geom_boxplot(
      width = 0.2,
      outlier.shape = NA,
      alpha = box.opacity
    )
  }

#geom_half_violin(
#  data = iris %>% filter(Species=="versicolor"), 
#  aes(x = Species, y = Sepal.Length, fill = Species), side="r")
  if (add.violin == TRUE) {
    p <- p + geom_half_violin(data = xmeta_lf_2 %>% 
                                filter(group==unique(xmeta_lf_2$group)[1]),
      position = position_nudge(x = -0.15, y = 0),
      alpha = violin.opacity,
      side = "l"
    ) +
      geom_half_violin(data = xmeta_lf_2 %>% 
                         filter(group==unique(xmeta_lf_2$group)[2]),
                       position = position_nudge(x = 0.15, y = 0),
                       alpha = violin.opacity,
                       side = "r"
      )
  }

  if (length(select.taxa) >= 2) {
    p <- p + facet_wrap(~taxa,
      scales = "free",
      ncol = ncol,
      nrow = nrow
    )
  }
  p <- p + scale_color_manual(
    values = c(
      Up = line.up,
      Down = line.down,
      Stable = line.stable
    ),
    guide = line.guide
  )

  p <- p + theme_biome_utils() + xlab(group)

  return(p)
}
