#' @title Plotting core microbiota on ordinations
#' @description This function will plot the ordaination along with highligthing the core microbes
#' on the species ordination.
#' @details This function is useful for visualizing core taxa in a 2D ordination plot.
#' @param x \code{\link{phyloseq-class}} object
#' @param ordiObject Output of ordinate from package phyloseq. Only NMDS and Bray supported.
#' @param prevalences Prevalences as supported by microbiome package.
#' @param detections Detections as supported by microbiome package.
#' @param min.prevalence Minimum prevalence value to plot.
#' @param color.opt Variable of interest from metadata.
#' @param shape Variable of interest from metadata.
#' @param Samples c("TRUE" or "FALSE")
#' @importFrom grDevices colorRampPalette
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggpubr ggarrange
#' @return plot
#' @export
#' @examples
#' \dontrun{
#' 
#' library(microbiomeutilities)
#' library(RColorBrewer)
#' data("zackular2014")
#' p0 <- zackular2014
#' ps1 <- format_to_besthit(p0)
#' ps1 <- subset_samples(ps1, DiseaseState == "H")
#' ps1 <- prune_taxa(taxa_sums(ps1) > 0, ps1)
#' prev.thres <- seq(.05, 1, .05)
#' det.thres <- 10^seq(log10(1e-4), log10(.2), length = 10)
#' pseq.rel <- microbiome::transform(ps1, "compositional")
#' ord.bray <- ordinate(pseq.rel, "NMDS", "bray")
#' 
#' p <- plot_ordiplot_core(pseq.rel, ord.bray,
#'   prev.thres, det.thres,
#'   min.prevalence = 0.9,
#'   color.opt = "DiseaseState", shape = NULL, Sample = TRUE
#' )
#' p
#' }
#' @keywords utilities
#'
plot_ordiplot_core <-
  function(x,
             ordiObject,
             prevalences,
             detections,
             min.prevalence,
             color.opt,
             shape,
             Samples = c(TRUE, FALSE)) {
    ordi <-
      spps <-
      OTU_ID <-
      prevdf0 <-
      prevdf1 <-
      subprevdf <-
      core <-
      NMDS1 <-
      NMDS2 <-
      Taxa_level <- tax.unit <- min.det <- list.sp <- tax <- Core <- NULL

    # Extract the data.frame or ordination
    proja <- plot_ordination(x, ordiObject, justDF = T)

    # Plot the sample ordination
    p0 <- plot_ordination(x, ordiObject, color = color.opt, shape = NULL)

    if (ordiObject$converged == TRUE) {
      message("Convergence reached in the ordination object provided")
    } else {
      message("Caution: Convergence was not reached in the ordination object provided")
    }

    p0 <- p0 + ggtitle(paste0("Samples NMDS plot ", "Stress: ", round(ordiObject$stress, 5)))

    # Extract species values

    spps <- as.data.frame(vegan::scores(ordiObject,
      display = "species"
    ))
    # Add names
    spps$OTU_ID <- rownames(spps)

    # Plot the core
    p1 <- plot_core(x,
      plot.type = "heatmap",
      prevalences = prevalences,
      detections = detections,
      min.prevalence = min.prevalence,
      horizontal = TRUE,
      colours = rev(brewer.pal(5, "Spectral"))
    )+ theme_bw()

    # get plot object
    prevdf0 <- p1$data

    prevdf1 <- as.data.frame(prevdf0)

    # get the list of OTUs
    list.sp <- prevdf1[, 1]

    # get the taxonomy data
    tax <- tax_table(x)
    tax <- as.data.frame(tax)

    # add the OTus to last column
    tax$OTU <- rownames(tax)

    # select taxonomy of only
    # those OTUs that are used in the plot
    tax2 <- dplyr::filter(tax, rownames(tax) %in% list.sp)

    # We will merge all the column into one except the Doamin as all is bacteria in this case

    Genus <- NULL
    Phylum <- NULL
    tax.unit <- tidyr::unite(tax2,
      Taxa_level,
      c("Genus"),
      sep = "_;",
      remove = TRUE
    )

    prevdf0[, 1] <- tax.unit$Taxa_level
    # p1$data <- prevdf0

    min.det <- min(prevdf1[, 2])

    min.prev <- min(prevdf1[, 3])

    p1 <- p1 + theme(axis.text.x = element_text(
      vjust = 0.5,
      hjust = 1,
      face = "italic"
    )) + theme_bw() + ggtitle(paste0(
      "Core taxa at minimum abundance of ",
      min.det, " and minimum prevalence ",
      min.prevalence
    )) 


    subprevdf <- subset(prevdf1, prevdf1[, 2] <= min.det)

    spps$Core <- match(spps$OTU_ID, subprevdf$Taxa, nomatch = 0)

    spps$Core <- as.logical(spps$Core)
    taxdf <- as.data.frame(tax_table(x))
    taxdf$OTU_ID <- rownames(taxdf)
    df1 <- merge(spps, taxdf, by.x = "OTU_ID", by.y = "OTU_ID")
    df2 <- mutate(df1,
      Core = ifelse(Core == TRUE, paste(Genus),
        df1$Core
      )
    )

    p2 <- ggplot(
      df2,
      aes(
        x = NMDS1, y = NMDS2,
        label = NA
      )
    ) + theme_bw() +
      geom_point(aes(color = Phylum), alpha = 0.5) +
      geom_text_repel(
        data = subset(
          x = df2,
          subset = Core != FALSE
        ),
        aes(label = Genus),
        alpha = 0.9,
        size = 3,
        fontface = "italic"
      ) + ggtitle("NMDS ordination for taxa")

    if (Samples == TRUE) {
      p.core <- ggarrange(p1, ggarrange(p0, p2), nrow = 2)
    } else {
      p.core <- ggarrange(p1, p2, nrow = 2)
    }

    return(p.core)
  }
