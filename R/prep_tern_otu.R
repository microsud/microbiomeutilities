#' @title Create table for ternary plot OTU
#' @description Create a table for ternary plot ggtern package.
#' @details Plots the mean relative abundance of taxa in 3 groups being compared.
#' @param x \code{\link{phyloseq-class}} object.
#' @param abund.thres = 0.0001 check \code{\link{microbiome}} package core function.
#' remove taxa that are detected at 0.0001 in less than prev.thres of samples.
#' @param prev.thres = 0.1 check \code{\link{microbiome}} package core function.
#' @param group Grouping variable to compare, for this plot there has to be three
#' groups in the data
#' @return  tibble object.
#' @examples
#' # library(microbiome)
#' # library(microbiomeutilities)
#' # library(dplyr)
#' # data("zackular2014")
#' # p0 <- zackular2014
#' # prep_tern_otu(p0, group="DiseaseState",
#' # abund.thres=0.0001, prev.thres=0.25)
#' @keywords utilities visualization
prep_tern_otu <- function(x,
                          abund.thres = 0.0001,
                          prev.thres = 0.1,
                          group = NULL) {
  p0.mr <- p0.mr.cr <- df <- taxa_df <- OTUID <- Sam_rep <- taxa_mean <- NULL
  Abundance <- value <- NULL

  if (is.null(group)) {
    message("Ternary plot requires 3 groups to compare")
    stop("Please provide group variable to compare")
  }

  if (length(unique(meta(x)[, group])) < 3) {
    message("Ternary plot requires 3 groups to compare")
    stop(paste0(group, " variable has less that 3 groups to compare"))
  }

  if (length(unique(meta(x)[, group])) > 3) {
    message("Ternary plot requires 3 groups to compare")
    stop(paste0(group, " variable has more than 3 groups to compare"))
  }


  p0.mr <- transform(x, "compositional")
  p0.mr.cr <- core(p0.mr,
    detection = abund.thres,
    prevalence = prev.thres
  )

  tax_table(p0.mr.cr)[, colnames(tax_table(p0.mr.cr))] <-
    gsub(tax_table(p0.mr.cr)[, colnames(tax_table(p0.mr.cr))],
      pattern = "[a-z]__", replacement = ""
    )

  p0.mr <- merge_samples(p0.mr.cr, group = group)

  df <- phy_to_ldf(p0.mr, NULL)

  # extract taxonomy
  taxa_df <- tax_table(p0.mr) %>%
    as("matrix") %>%
    as.data.frame()
  taxa_df$OTUID <- rownames(taxa_df)

  # taxa_df[,level.color][which(taxa_df[,level.color] == "")] <- "Other"

  # taxa_df$level_var <- taxa_df[,level]

  taxa_mean <- df %>%
    group_by(OTUID, Sam_rep) %>%
    summarise(value = mean(Abundance))

  taxa_tern <- smas <- NULL

  taxa_tern <- taxa_mean %>%
    group_by(Sam_rep, OTUID) %>%
    # mutate(index=row_number()) %>%
    pivot_wider(names_from = Sam_rep, values_from = value) %>%
    left_join(taxa_df, by = "OTUID")

  return(taxa_tern)
}
