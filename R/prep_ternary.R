#' @title Create table for Ternary plot
#' @description Create a table for ternary plot ggtern R package.
#' @details Plots the mean relative abundance of taxa in 3 groups being compared.
#' @param x \code{\link{phyloseq-class}} object
#' @param abund.thres = 0.0001 check \code{\link{microbiome}} package core function
#' remove taxa that are dectected at 0.0001 in less than prev.thres of samples
#' @param prev.thres = 0.1 check \code{\link{microbiome}} package core function
#' @param level = "Genus" Taxonomic level. If OTU/ASV level specify="lowest"
#' Does not support phylum level aggregation
#' @param group Grouping variable to compare, for this plot there has to be three
#' groups in the data
#' @return  Tibble object.
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#' library(microbiome)
#' library(microbiomeutilities)
#' library(dplyr)
#' data("zackular2014")
#' p0 <- zackular2014
#' prep_ternary(p0, group = "DiseaseState", abund.thres = 0.0001, level = "Genus", prev.thres = 0.25)
#' @keywords visualization

prep_ternary <- function(x,
                         abund.thres = 0.0001,
                         prev.thres = 0.1,
                         group = NULL,
                         level = "lowest") {
  p0.agg <- Sam_rep <- Abundance <- value <- taxa_df <- OTUID <- taxa_mean <- NULL
  p0.mr.cr <- level_var <- smas <- taxa_tern <- NULL

  if (level == "lowest") {
    p <- prep_tern_otu(x,
      abund.thres = abund.thres,
      prev.thres = prev.thres,
      group = group
    )
    return(p)
  } else {
    # x <- p0
    p0.mr <- transform(x, "compositional")
    tax_table(p0.mr)[, colnames(tax_table(p0.mr))] <-
      gsub(tax_table(p0.mr)[, colnames(tax_table(p0.mr))],
        pattern = "[a-z]__", replacement = ""
      )

    p0.agg <- aggregate_taxa(p0.mr, level = level)
    p0.mr.cr <- core(p0.agg,
      detection = abund.thres,
      prevalence = prev.thres
    )
    # DT::datatable(tax_table(p0.mr.cr))

    p0.mr <- merge_samples(p0.mr.cr, group = group)

    # extract taxonomy
    taxa_df <- tax_table(p0.mr) %>%
      as("matrix") %>%
      as.data.frame()
    taxa_df$OTUID <- rownames(taxa_df)
    # taxa_df$unique <- taxa_df[,unique]

    df <- phy_to_ldf(p0.mr, NULL)
    # df$level_var <- df[,unique]
    # variable goes to Sam_rep column
    taxa_mean <- df %>%
      group_by(unique, Sam_rep) %>%
      summarise(value = mean(Abundance))

    taxa_tern <- taxa_mean %>%
      group_by(Sam_rep, unique) %>%
      # mutate(index=row_number()) %>%
      pivot_wider(names_from = Sam_rep, values_from = value) %>%
      left_join(taxa_df, by = "unique")

    return(taxa_tern)
  }
}
