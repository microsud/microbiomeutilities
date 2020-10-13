#' @title Mean Abundance-Prevalence relation
#' @description Plots Mean Abundance-Prevalence for taxa. Mean abundance, mean prevalence,
#' and upper and lower confidence interval for each taxa is calculated by random subsampling.
#' @details Check if there are spurious OTUs/ASVs.
#' @param x \code{\link{phyloseq-class}} object
#' @param lower.conf Lower confidence interval =0.025
#' @param upper.conf Upper confidence interval =0.975
#' @param bs.iter Number of bootstrap iterations =99
#' @param color taxa level to color. preferablly at phylum
#' @param point.opacity  Numeric for ggplot alpha. Defualt is 0.5
#' @param label.core Logical default is FALSE
#' @param label.size If label_core is TRUE specify text size.
#' Default is NULL
#' @param label.opacity Numeric for ggplot alpha. Defualt is NULL
#' @param nudge.label Argument to pass to ggrepel::geom_text_repel Default is NULL
#' @param mean.abund.thres If label_core is TRUE specify mean
#' abundance threshold. Default is NULL
#' @param mean.prev.thres If label_core is TRUE specify mean
#' prevalence threshold. Default is NULL
#' @param log.scale Plot log10 scale. Default is TRUE
#' abundance criteria. Default is NULL
#' @param ... Arguments to pass to sample() function.
#' @return A \code{\link{ggplot}} plot object.
#' @importFrom stats quantile
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples
#' \dontrun{
#' # Example data
#' library(microbiomeutilities)
#' asv_ps <- zackular2014
#' asv_ps <- microbiome::transform(asv_ps, "compositional")
#' asv_ps <- core(asv_ps, detection = 0.0001, prevalence = 0.5)
#' asv_ps <- format_to_besthit(asv_ps)
#' set.seed(2349)
#' p_v <- plot_abund_prev(asv_ps, size = 20, replace = TRUE) +
#'   geom_vline(xintercept = 0.75, lty = "dashed", alpha = 0.7) +
#'   geom_hline(yintercept = 0.01, lty = "dashed", alpha = 0.7) +
#'   scale_color_brewer(palette = "Paired")
#' p_v
#' }
#'
#' @keywords visualization analysis

plot_abund_prev <- function(x, lower.conf = 0.025,
                            upper.conf = 0.975,
                            bs.iter = 99,
                            color = "Phylum",
                            point.opacity = 0.5,
                            label.core = FALSE,
                            label.size = NULL,
                            label.opacity = NULL,
                            mean.abund.thres = NULL,
                            mean.prev.thres = NULL,
                            log.scale = TRUE,
                            nudge.label = NULL,
                            ...) {
  psx <- rand_sams <- ps.sub <- sub.sum <- txvp <- sxi <- NULL
  s <- ab <- sx <- cis <- taxsp_lc <- taxsp_uc <- cis_df <- NULL
  abx <- cis_ab <- tax_list <- cis_ab_df <- abx_abun <- NULL
  Mean.Rel.Ab <- MeanAbun <- Taxa <- NULL

  ci_ab_prev_tax <- tax_df <- core_df <- NULL
  message("Make sure to set.seed")
  psx <- x
  s <- c()
  ab <- c()
  for (i in seq_len(bs.iter)) {
    rand_sams <- sample(sample_names(psx), ...)
    ps.sub <- prune_samples(sample_names(psx) %in% rand_sams, psx)
    sub.sum <- abun_summary(ps.sub)
    # ps.sub <- prune_taxa(taxa_sums(ps.sub) > 0, ps.sub)
    txvp <- prevalence(ps.sub, detection = 0, sort = TRUE, count = F)
    # rownames(txvp) <- txvp$Taxa
    s[[i]] <- txvp

    ab[[i]] <- sub.sum
  }
  sx <- dplyr::bind_rows(s)
  sx <- t(sx) %>% as.data.frame()
  # head(sx)
  cis <- c()
  for (tax in rownames(sx)) {
    taxsp_lc <- quantile(sx[tax, ], lower.conf, na.rm = TRUE)
    taxsp_uc <- quantile(sx[tax, ], upper.conf, na.rm = TRUE)
    cis <- rbind(cis, c(tax, taxsp_lc, taxsp_uc))
  }
  cis_df <- as.data.frame(cis)
  colnames(cis_df) <- c("Taxa", "PrevLowerCI", "PrevUpperCI")
  sx$meanPrev <- rowMeans(sx)
  meanPrev <- sx[, "meanPrev"]
  sxi <- cbind(meanPrev, cis_df)
  # head(sxi)

  abx <- dplyr::bind_rows(ab) %>%
    as.data.frame() %>%
    group_by(Taxa)

  cis_ab <- NULL
  tax_list <- unique(abx$Taxa)
  for (taxa_ls in tax_list) {
    abx.sub <- subset(abx, Taxa == taxa_ls)
    taxsp_ab_lc <- quantile(abx.sub[, "Mean.Rel.Ab"], lower.conf, na.rm = TRUE)

    taxsp_ab_uc <- quantile(abx.sub[, "Mean.Rel.Ab"], upper.conf, na.rm = TRUE)

    cis_ab <- rbind(cis_ab, c(taxa_ls, taxsp_ab_lc, taxsp_ab_uc))
  }
  cis_ab_df <- as.data.frame(cis_ab)

  colnames(cis_ab_df) <- c("Taxa", "MeanAbunLowerCI", "MeanAbunUpperCI")

  abx_abun <- abx %>%
    group_by(Taxa) %>%
    summarise(MeanAbun = mean(Mean.Rel.Ab))

  ci_ab_prev <- cis_ab_df %>%
    left_join(sxi) %>%
    right_join(abx_abun)

  tax_df <- tax_table(psx) %>%
    as("matrix") %>%
    as.data.frame()
  tax_df$Taxa <- rownames(tax_df)

  ci_ab_prev_tax <- ci_ab_prev %>%
    left_join(tax_df)
  # head(ci_ab_prev)

  ci_ab_prev_tax$MeanAbunLowerCI <- as.numeric(ci_ab_prev_tax$MeanAbunLowerCI)
  ci_ab_prev_tax$MeanAbunUpperCI <- as.numeric(ci_ab_prev_tax$MeanAbunUpperCI)
  ci_ab_prev_tax$PrevLowerCI <- as.numeric(ci_ab_prev_tax$PrevLowerCI)
  ci_ab_prev_tax$PrevUpperCI <- as.numeric(ci_ab_prev_tax$PrevUpperCI)

  # hist(ci_ab_prev_tax$meanPrev)
  # hist(ci_ab_prev_tax$MeanAbun)


  p <- ggplot(ci_ab_prev_tax, aes(meanPrev, MeanAbun)) +
    geom_point(aes_string(color = color), alpha = 0.5, size = 2) +
    geom_linerange(aes_string(
      ymin = "MeanAbunLowerCI",
      ymax = "MeanAbunUpperCI",
      color = color
    ), alpha = point.opacity) +
    geom_linerange(aes_string(
      xmin = "PrevLowerCI",
      xmax = "PrevUpperCI",
      color = color
    ), alpha = point.opacity)

  if (label.core == TRUE) {
    core_df <- subset(
      ci_ab_prev_tax,
      meanPrev >= mean.prev.thres & MeanAbun >= mean.abund.thres
    )
    p <- p + ggrepel::geom_text_repel(
      data = core_df,
      aes(label = Taxa),
      size = label.size,
      alpha = label.opacity,
      force = 0.5,
      nudge_x = nudge.label,
      direction = "y",
      hjust = 1,
      segment.size = 0.2
    )
  }

  if (log.scale == TRUE) {
    p <- p + scale_y_log10() + scale_x_log10() +
      xlab("Prevalance (log10)") +
      ylab("Mean abundance (log10)")
    return(p + theme_biome_utils())
  } else {
    p <- p + xlab("Prevalance") + ylab("Mean abundance")
    return(p + theme_biome_utils())
  }
}
