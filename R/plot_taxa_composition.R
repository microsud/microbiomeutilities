#' @title Taxonomic Composition Plot
#' @description Plot taxon abundance for samples. It is a legacy function from \code{\link{microbiome}}.
#' @param x \code{\link{phyloseq-class}} object
#' @param taxonomic.level Merge the OTUs (for phyloseq object) into a higher taxonomic level. This has to be one from colnames(tax_table(x)).
#' @param sample.sort Order samples. Various criteria are available:
#'   \itemize{
#'     \item NULL or 'none': No sorting
#'     \item A single character string: indicate the metadata field to be used for ordering
#'     \item A character vector: sample IDs indicating the sample ordering.
#'     \item 'neatmap' Order samples based on the neatmap approach. See \code{\link{neatsort}}. By default, 'NMDS' method with 'bray' distance is used. For other options, arrange the samples manually with the function.
#'   }
#' @param otu.sort Order taxa. Same options as for the sample.sort argument but instead of metadata, taxonomic table is used. Also possible to sort by 'abundance'.
#' @param x.label Specify how to label the x axis. This should be one of the variables in sample_variables(x).
#' @param plot.type Plot type: 'barplot' or 'lineplot'.
#' @param verbose verbose.
#' @param transform Data transform to be used in plotting (but not in sample/taxon ordering). The options are 'Z-OTU', 'Z-Sample', 'log10' and 'compositional'. See the \code{\link{transform}} function.
#' @param mar Figure margins.
#' @param average_by Variable to group.
#' @param palette The number and palette \code{\link{RColorBrewer}} has to be specified e.g brewer.pal(12, "Paired").
#' @param ... Arguments to be passed (for \code{\link{neatsort}} function)
#' @return A \code{\link{ggplot}} plot object.
#' @export
#' @examples \dontrun{
#'     # Example data
#'     library(microbiome)
#'     library(microbiomeutilities)
#'     data("biogeogut")
#'     pseq <- biogeogut
#'     plot_taxa_composition(pseq, taxonomic.level = "Phylum")
#'           }
#' @keywords utilities
plot_taxa_composition <- function (x, sample.sort = NULL, taxonomic.level, transform, otu.sort = NULL, palette = brewer.pal(12, "Paired"), x.label = "sample",
                               plot.type = "barplot", average_by, verbose = FALSE, mar = c(5, 12, 1,
                                                                               1), ...)
{
  Sample <- Abundance <- Taxon <- horiz <- value <- scales <- ID <- meta <- OTU <- taxic <- otu.df <- taxmat <- new.tax <- NULL
  if (!is.null(x@phy_tree)) {
    x@phy_tree = NULL
  }

  taxic <- x@tax_table
  taxic <- as.data.frame.matrix(taxic)
  otu.df <- abundances(x)
  otu.df <- as.data.frame.matrix(otu.df)
  taxic$OTU <- row.names(otu.df)
  taxmat <- as.matrix(taxic)
  new.tax <- tax_table(taxmat)
  tax_table(x) <- new.tax
  xorig <- x

  # Merge the taxa at a higher taxonomic level
  if (!taxonomic.level == "OTU") {
    if (verbose) {
      message("Aggregating the taxa.")
    }
    x <- aggregate_taxa(x, taxonomic.level)
  }

  if (verbose) {
    message("Check data transforms.")
  }

  if (is.null(transform)) {
    x <- x
  } else if (transform == "Z-OTU") {
    x <- microbiome::transform(x, "Z", "OTU")
  } else if (transform == "Z-Sample") {
    x <- microbiome::transform(x, "Z", "Sample")
  } else if (transform == "compositional") {
    x <- microbiome::transform(x, "compositional")
  } else {
    x <- microbiome::transform(x, transform)
  }

  abu <- abundances(x)
  group <- NULL
  if (!is.null(average_by)) {
    dff <- as.data.frame(t(abu))
    dff$group <- sample_data(x)[[average_by]]
    if (is.numeric(dff$group)) {
      dff$group <- factor(dff$group, levels = sort(unique(dff$group)))
    }
    dff <- dff %>% filter(!is.na(group))
    dff$group <- droplevels(dff$group)
    av <- aggregate(. ~ group, data = dff, mean)
    rownames(av) <- as.character(av$group)
    av$group <- NULL
    abu <- t(av)
  }
  if (is.null(sample.sort) || sample.sort == "none" || !is.null(average_by)) {
    sample.sort <- colnames(abu)
  }
  else if (length(sample.sort) == 1 && sample.sort %in% names(sample_data(x)) &&
           is.null(average_by)) {
    sample.sort <- rownames(sample_data(x))[order(sample_data(x)[[sample.sort]])]
  }
  else if (all(sample.sort %in% sample_names(x)) & is.null(average_by)) {
    sample.sort <- sample.sort
  }
  else if (length(sample.sort) == 1 && sample.sort == "neatmap") {
    sample.sort <- neatsort(x, method = "NMDS", distance = "bray",
                            target = "sites", first = NULL)
  }
  else if (!sample.sort %in% names(sample_data(x))) {
    warning(paste("The sample.sort argument", sample.sort,
                  "is not included in sample_data(x). \n            Using original sample ordering."))
    sample.sort <- sample_names(x)
  }
  if (is.null(otu.sort) || otu.sort == "none") {
    otu.sort <- taxa(x)
  }
  else if (length(otu.sort) == 1 && otu.sort == "abundance") {
    otu.sort <- rev(names(sort(rowSums(abu))))
  }
  else if (length(otu.sort) == 1 && otu.sort %in% names(tax_table(x))) {
    otu.sort <- rownames(sample_data(x))[order(tax_table(x)[[otu.sort]])]
  }
  else if (all(otu.sort %in% taxa(x))) {
    otu.sort <- otu.sort
  }
  else if (length(otu.sort) == 1 && otu.sort == "neatmap") {
    otu.sort <- neatsort(x, method = "NMDS", distance = "bray",
                         target = "species", first = NULL)
  }
  if (verbose) {
    message("Prepare data.frame.")
  }
  dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
  names(dfm) <- c("OTU", "Sample", "Abundance")
  dfm$Sample <- factor(dfm$Sample, levels = sample.sort)
  dfm$OTU <- factor(dfm$OTU, levels = otu.sort)

  colourCount = length(unique(dfm$OTU))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
  getPalette = colorRampPalette(palette)

  if (x.label %in% colnames(sample_data(x)) & is.null(average_by)) {
    meta <- sample_data(x)
    dfm$xlabel <- as.vector(unlist(meta[as.character(dfm$Sample),
                                        x.label]))
    if (is.factor(meta[, x.label])) {
      lev <- levels(meta[, x.label])
    }
    else {
      lev <- unique(as.character(unname(unlist(meta[,
                                                    x.label]))))
    }
    dfm$xlabel <- factor(dfm$xlabel, levels = lev)
  }
  else {
    dfm$xlabel <- dfm$Sample
  }
  if (verbose) {
    message("Construct the plots")
  }
  if (plot.type == "barplot") {
    dfm <- dfm %>% arrange(OTU)
    p <- ggplot(dfm, aes(x = Sample, y = Abundance, fill = OTU))
    p <- p + geom_bar(position = "stack", stat = "identity")
    p <- p + scale_x_discrete(labels = dfm$xlabel, breaks = dfm$Sample)
    p <- p + ylab("Abundance") +  scale_fill_manual(taxonomic.level, values = getPalette(colourCount))
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              vjust = 0.5, hjust = 0))
    p <- p + guides(fill = guide_legend(reverse = FALSE))
  }
  else if (plot.type == "lineplot") {
    dfm <- dfm %>% arrange(OTU)
    p <- ggplot(dfm, aes(x = Sample, y = Abundance, color = OTU,
                         group = OTU))
    p <- p + geom_point()
    p <- p + geom_line() + scale_color_brewer(guide = guide_legend(title = taxonomic.level))
    p <- p + scale_x_discrete(labels = dfm$xlabel, breaks = dfm$Sample)
    if (!is.null(transform) && transform == "compositional") {
      suppressMessages(p <- p + ylab("Relative abundance (%)") +  scale_color_manual(taxonomic.level, values = getPalette(colourCount)))
    }
    else {
      p <- p + ylab("Abundance")
    }
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              vjust = 0.5, hjust = 0))
    suppressMessages(p <- p + guides(fill = guide_legend(reverse = FALSE)) +
                       scale_color_manual(taxonomic.level, values = getPalette(colourCount)))
  }
  p + theme_bw() + theme(axis.text.x = element_text(face ="italic", angle = 90))
}
