#' @title Diversity plot with stats
#' @param x \code{\link{phyloseq-class}} object
#' @param index diversity index. Calculated using microbiome::alpha
#' @param group Grouping variable to compare
#' @param group.colors Colors for plotting groups
#' @param dot.opacity for ggplot alpha to determine opacity for points
#' @param box.opacity for ggplot alpha to determine opacity for box
#' @param violin.opacity for ggplot alpha to determine opacity for violin
#' @param group.order Default is NULL. a list specifing order of x-axis.
#' @param dot.size size of point
#' E.g. c("H","CRC","nonCRC")
#' @param ... params for ggpubr::stat_compare_means
#' @importFrom ggpubr stat_compare_means rotate_x_text
#' @importFrom gghalves geom_half_violin
#' @keywords utilties
plot_diversity_with_stats <- function(x, index,
                                      group = NULL,
                                      group.colors = c("brown3", "steelblue"),
                                      dot.opacity = 0.25,
                                      box.opacity = 0.25,
                                      violin.opacity = 0.5,
                                      dot.size = 2,
                                      group.order = NULL, ...) {

  # x <- ps.c
  adiv <- alpha(x, index = index)
  df_met <- cbind(adiv, meta(x))

  if (!is.null(group.order)) {
    df_met[, group] <- factor(df_met[, group],
      levels = group.order
    )
  }

  plt <- ggplot(
    data = df_met,
    aes_string(group, index, fill = group)
  ) +
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      alpha = violin.opacity, side = "r"
    ) +
    geom_point(aes_string(y = index, color = group),
      position = position_jitter(width = 0.15),
      size = dot.size, alpha = dot.opacity
    ) +
    geom_boxplot(
      width = 0.2, outlier.shape = NA,
      alpha = box.opacity
    ) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values = group.colors) +
    scale_colour_manual(values = group.colors)

  if (length(unique(df_met[, group])) == 2) {
    plt <- plt + stat_compare_means(
      label = "p.format",
      tip.length = 0.05,
      method = "wilcox.test", ...
    )
  } else if (length(unique(df_met[, group])) > 2) {
    comps <- make_pairs(df_met[, group])
    plt <- plt + stat_compare_means(
      comparisons = comps,
      label = "p.format",
      tip.length = 0.05,
      method = "wilcox.test", ...
    )
  }

  plt + theme_biome_utils() + rotate_x_text()
}



##################################################################################
#' @title Diversity plot without stats
#' @param x \code{\link{phyloseq-class}} object
#' @param index diversity index. Calculated using microbiome::alpha
#' @param group Grouping variable to compare
#' @param group.colors Colors for plotting groups
#' @param dot.opacity for ggplot alpha to determine opacity for points
#' @param box.opacity for ggplot alpha to determine opacity for box
#' @param dot.size size of point
#' @param violin.opacity for ggplot alpha to determine opacity for violin
#' @param group.order Default is NULL. a list specifing order of x-axis.
#' E.g. c("H","CRC","nonCRC")
#' @param ... params for ggpubr::stat_compare_means
#' @importFrom ggpubr stat_compare_means rotate_x_text
#' @importFrom gghalves geom_half_violin
plot_diversity_without_stats <- function(x, index,
                                         group = NULL,
                                         group.colors = c("brown3", "steelblue"),
                                         dot.opacity = 0.25,
                                         box.opacity = 0.25,
                                         violin.opacity = 0.5,
                                         dot.size = 2,
                                         group.order = NULL, ...) {

  # x <- ps.c
  adiv <- alpha(x, index = index)
  df_met <- cbind(adiv, meta(x))

  if (!is.null(group.order)) {
    df_met[, group] <- factor(df_met[, group],
      levels = group.order
    )
  }

  plt <- ggplot(
    data = df_met,
    aes_string(group, index, fill = group)
  ) +
    geom_half_violin(
      position = position_nudge(x = 0.15, y = 0),
      alpha = violin.opacity, side = "r"
    ) +
    geom_point(aes_string(y = index, color = group),
      position = position_jitter(width = 0.15),
      size = dot.size, alpha = dot.opacity
    ) +
    geom_boxplot(
      width = 0.2, outlier.shape = NA,
      alpha = box.opacity
    ) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values = group.colors) +
    scale_colour_manual(values = group.colors)
  plt <- plt + theme_biome_utils() + rotate_x_text()
  return(plt)
}
