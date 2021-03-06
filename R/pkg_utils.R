#' @title Custom theme for microbiomeutilities pkg
#' @description Opiniated elegant theme.
#' @export
#' @keywords utilities
theme_biome_utils <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = element_text(size = 12),
      # axis.title = element_text(size = 18),
      axis.line.x = element_line(color = "#242424"),
      axis.line.y = element_line(color = "#242424"),
      # panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "#242424"),
      strip.background = element_rect(fill = "white", color = "#242424"),
      strip.text = element_text(color = "#353535"),
      legend.key = element_blank(),
      legend.text = element_text(color = "#353535"),
      legend.background = element_rect(colour = NA, fill = NA)
    )
}

#' @title Pool Taxa
#' @description Creates a list of dataframes at different taxonomic levels.
#' @param x \code{\link{phyloseq-class}} object
#' @examples
#' # library(phyloseq)
#' # library(microbiome)
#' # data("zackular2014")
#' # sub.sm <- sample(sample_names(zackular2014), 20)
#' # pseq <- prune_samples(sub.sm,zackular2014)
#' # taxa_pooler_mcola(pseq)
#' @keywords utilities
taxa_pooler_mcola <- function(x) {
  tax_table(x)[, colnames(tax_table(x))] <- gsub(tax_table(x)[, colnames(tax_table(x))],
    pattern = "[a-z]__", replacement = " "
  )
  lev_tax <- pool_tabs <- otu_tab <- OTU <- NULL
  message("Creating a list of tables for ")
  message(cat(rank_names(x)))
  lev_tax <- rank_names(x)
  # lev_tax <- lev_tax[-length(lev_tax)]
  pool_tabs <- c()
  for (i in lev_tax) {
    tax_table(x)[tax_table(x)[, i] == " ", i] <- NA
    otu_tab <- t(abundances(aggregate_taxa(x, i)))
    pool_tabs[[i]] <- otu_tab
  }
  OTU <- length(lev_tax) + 1
  pool_tabs[[OTU]] <- t(abundances(x))
  return(pool_tabs)
}

#' @title Distribution of taxa
#' @description Plots distribution of taxa.
#' @param x \code{\link{phyloseq-class}} object
#' @param color.level Taxonomic level to color
#' @param color.taxa vector of colors specified by user
#' Default is brewer.pal(12,"Paired")
#' @return ggplot2 object
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' library(microbiomeutilities)
#' data("zackular2014")
#' pseq <- zackular2014
#' p <- taxa_distribution(pseq)
#' p
#' @export
#' @keywords utilities
taxa_distribution <- function(x, color.level = "Phylum",
                              color.taxa = brewer.pal(12, "Paired")) {
  taxasums <- taxatable <- tax_plot1 <- NULL

  tax_table(x)[, colnames(tax_table(x))] <- gsub(tax_table(x)[, colnames(tax_table(x))],
    pattern = "[a-z]__", replacement = ""
  )

  tax_table(x)[tax_table(x)[, color.level] == "", color.level] <- "Unclassified"
  taxasums <- rowSums(abundances(x))
  taxatable <- tax_table(x) %>%
    as("matrix") %>%
    as.data.frame()
  taxatable$taxasums <- taxasums

  tax_plot1 <- ggplot(taxatable, aes_string(
    x = "taxasums",
    color = color.level
  )) +
    geom_line(size = 0.5, stat = "density") +
    xlab(paste0(color.level, " counts (log10)")) +
    scale_x_log10() +
    scale_color_manual(color.level, values = color.taxa)
  return(tax_plot1)
}

# @title Summarize abundance
# @param x \code{\link{phyloseq-class}} object
# @keywords internal
#'
abun_summary <- function(x) {
  Max.Rel.Ab <- Mean.Rel.Ab <- MeanAbun <- Median.Rel.Ab <- Std.dev <- Taxa <- NULL
  otudf2 <- as.data.frame(abundances(x))

  output <- NULL
  for (j in 1:nrow(otudf2)) {
    x2 <- as.numeric(otudf2[j, ])
    mx.rel <- max(x2)
    mean.rel <- mean(x2)
    med.rel <- median(x2)
    std.dev <- sd(x2)

    output <- rbind(output, c(
      row.names(otudf2)[j],
      as.numeric(mx.rel),
      as.numeric(mean.rel),
      as.numeric(med.rel),
      as.numeric(std.dev)
    ))
  }

  # head(output)
  outputdf <- as.data.frame(output, stringsAsFactors = F)
  colnames(outputdf) <- c("Taxa", "Max.Rel.Ab", "Mean.Rel.Ab", "Median.Rel.Ab", "Std.dev")
  outputdf <- outputdf %>%
    mutate_at(vars(Max.Rel.Ab, Mean.Rel.Ab, Median.Rel.Ab, Std.dev), as.numeric)
}


#' @title Make pairs
#' @description Creates a combination of variables
#' for use with ggpubr::stat_compare_means.
#' @param x list of vector to compare
#' @examples
#' # library(microbiomeutilities)
#' # data("zackular2014")
#' # pseq <- zackular2014
#' # comps <- make_pairs(meta(pseq)$DiseaseState)
#' @keywords utilities
#' @export
make_pairs <- function(x) {
  if (is.character(x) == TRUE) {
    # message("is char")
    var.lev <- unique(x)
  } else if (is.factor(x) == TRUE) {
    # message("is fac")
    var.lev <- levels(x)
  }
  # make a pairwise list that we want to compare.
  lev.pairs <- combn(seq_along(var.lev), 2, simplify = FALSE, FUN = function(i) var.lev[i])
  return(lev.pairs)
}


#   NOT WORKING
# mat: the matrix passed to pheatmap
#   names_fun: either rownames or colnames
#   list_names: vector of names that should appear in boldface
#' #@param mat otu table rows as taxa names
#' #@param names_fun either rownames or colnames
#' #@param list_names vector of names that should appear in boldface
#make_font_face <- function(mat, names_fun, list_names, bold=TRUE) {
#  #names.tx <- names_fun(mat)
#  names.tx <- NULL
#  names.tx <- lapply(
#    rownames(mat),
#    function(x) bquote(italic(.(x))))
  
#  ids <- list_names %>% match(names_fun(mat))
#  if(bold==TRUE){
#    ids %>%
#      lapply(
#        function(i)
#          list_names[i] <<-
#          bquote(bold(.(names_fun(mat)[i]))) %>%
#          as.expression()
#      )
#    names.tx
#  }
#}


#' @keywords internal 
rarefy_util <- function(ps, sample.size, index){
  
  ps.rar <- vegan::rrarefy(t(abundances(ps)),sample.size)
  
  adiv <- suppressMessages(alpha(t(ps.rar), index=index))
  
  adiv <- as.data.frame(adiv)
  
  colnames(adiv) <- index
  
  adiv_nw <- cbind(adiv, meta(ps))
  
  adiv_nw$sub_sample <- sample.size
  
  rownames(adiv_nw) < seq(1:nrow(adiv_nw))
  
  return(adiv_nw)
}
