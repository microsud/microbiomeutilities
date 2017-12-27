#' @title Plotting core microibota on ordinations
#' @description This function will plot the ordaination alsong with highligthing the core microbes
#' on the species ordination.
#' @details Most commonly it is observed that the taxonomy file has classification until a given
#'          taxonomic level.
#'          Hence, to avoid loss of OTU information while using the function tax_glom() for merging #'          at a specific taxonomic level.
#'          we will fill the empty cells with the maximum classification available along with the
#'          #' OTU number. This code is a slight modification.
#'          the code from  \pkg{ampvis} \code{\link{phyloseq-class}}. Here, we directly take the
#'          #' phyloseq object as input and make the necessary formatting.
#' @param x \code{\link{phyloseq-class}} object
#' @param ordiObject Output of ordinate from package phyloseq. Only NMDS and Bray supported.
#' @param coreplot TRUE or FALSE for core heatmap from microbiome package.
#' @param prevalences Prevalences as supported by microbiome package.
#' @param detections Detections as supported by microbiome package.
#' @param color Variable of interest from metadata.
#' @param shape Variable of interest from metadata..
#' @return plot
#' @export
#' @examples \dontrun{
#'   # Example data
#'     data(DynamicsIBD)
#'     ps <- DynamicsIBD
#'     ps1 <- format_phyloseq(ps)
#'     px <- plot_ordiplot_core(ps1, ordiObject = ordi, coreplot = TRUE, prevalences = prev.thres,
#'     detections = det.thres, color = "ibd_subtype", shape = NULL)
#'           }
#' @keywords utilities
#'
plot_ordiplot_core <-
  function(x,
           ordiObject,
           coreplot = c(TRUE, FALSE),
           prevalences,
           detections,
           color.opt,
           shape) {
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


    proja <- plot_ordination(x, ordiObject, justDF = T)
    p0 <- plot_ordination(x, ordiObject, color = color.opt, shape = NULL)

    p0 <- p0 + ggtitle("Sites")
    spps  <-
      as.data.frame(vegan::scores(ordiObject, display = "species"))

    spps$OTU_ID <- rownames(spps)

    p1 <- plot_core(
      x,
      plot.type = "heatmap",
      prevalences = prevalences,
      detections = detections,
      min.prevalence = .90,
      horizontal = TRUE,
      colours = rev(brewer.pal(5, "Spectral"))
    )
    prevdf0 <- p1$data

    prevdf1 <- as.data.frame(prevdf0)

    # get the list of OTUs
    list.sp <- prevdf1[, 1]
    # check the OTU ids
    # print(list)

    # get the taxonomy data
    tax <- tax_table(x)
    tax <- as.data.frame(tax)

    # add the OTus to last column
    tax$OTU <- rownames(tax)

    # select taxonomy of only
    # those OTUs that are used in the plot
    tax2 <- dplyr::filter(tax, rownames(tax) %in% list.sp)


    # We will merege all the column into one except the Doamin as all is bacteria in this case
    tax.unit <-
      tidyr::unite(tax2,
                   Taxa_level,
                   c("Genus"),
                   sep = "_;",
                   remove = TRUE)
    prevdf0[, 1] <- tax.unit$Taxa_level
    p1$data <- prevdf0
    p1 <-
      p1 + theme(axis.text.x = element_text(
        vjust = 0.5,
        hjust = 1,
        face = "italic"
      ))

    min.det <- min(prevdf1[, 2])
    subprevdf <- subset(prevdf1, prevdf1[, 2] <= min.det)

    spps$Core <- match(spps$OTU_ID, subprevdf$Taxa, nomatch = 0)
    spps$Core = as.logical(spps$Core)
    taxdf <- as.data.frame(tax_table(x))
    taxdf$OTU_ID <- rownames(taxdf)
    df1 <- merge(spps, taxdf, by.x = "OTU_ID", by.y = "OTU_ID")
    df2 <-
      mutate(df1, Core = ifelse(Core == TRUE, paste(Genus), df1$Core))
    p2 <- ggplot(df2, aes(x = NMDS1, y = NMDS2, label = NA)) + theme_bw() +
      ggrepel::geom_text_repel(
        data = subset(x = df2, subset =  Core != FALSE),
        aes(label = Genus),
        alpha = 0.5,
        size = 3,
        fontface = "italic"
      ) +
      geom_point(aes(color = Phylum), alpha = 0.7) + ggtitle("Detection limit set is min " , subtitle = min.det)


    if (coreplot==TRUE) {
      p3 <- ggarrange(p1, ggarrange(ncol = 2, p0, p2), nrow = 2)
      return(p3)
    } else {
      p4 <- ggarrange(p0, p2)
      return(p4)
    }
  }
