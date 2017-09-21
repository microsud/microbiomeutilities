#' @title Microbiome analysis pipeline
#' @description This is a single workflow for prelimanry investigation of microbiome data.
#' @param biom biom object object
#' @param mapping Metadata variable to check for groups based sequencing effort.
#' @param taxonomy either density or histogram plot
#' @param treefilename tree file
#' @param out_dir otuput directory
#' @param VariableA Main variable of interest.
#' @param VariableB Secondary variable of interest
#' @param filterpseq TRUE or FALSE
#' @param Distmethod Distance method 'bray', 'unifrac', 'wunifrac', 'jsd'
#' @param OrdinatMethod Ordination method 'NMDS', 'MDS', 'PCoA'
#' @param QC Default TRUE QC plots
#' @param filterCount Filter OTUs below this count number
#' @param filterPrev Filter OTUs not detected in less than this percent of samples
#' @param transformation  'compositional', 'Z', 'log10', 'hellinger', 'clr' choose to tranform your counts data defaults to FALSE
#' @param col.palette 'Spectral', 'Set2', 'Paired' use any of the RColorbrewer palette depending on number of groups you have in VaraibleA
#'
#' @author Contact: SUdarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @return A \code{\link{ggplot}} plot object
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data(DynamicIBD)
#'     ps0 <- DynamicIBD
#'     p <- plot_ReadDistribution(ps0, groups='ibd_subtype', plot.type= 'density')
#'           }
#' @keywords utilities


run_microbiome_pipeline <- function(biom, mapping, taxonomy, treefilename,
                                    out_dir, VariableA, VariableB, Distmethod, OrdinatMethod,
                                    QC, filterCount, filterPrev, transformation, col.palette, filterpseq)
{

  ps0 <- qc_plot <- ps1 <- ps2 <- Variance.plot.a <- Variance.plot.b <- otu_tab_ps2 <- phy_tree_ps2 <- df.pd <- metadf <- metatab <- tab <- alpha_div <- alpha.plot <- ord <- NULL
  oripath <- getwd()
  setwd(out_dir)
  dir.create("QC")
  dir.create("AlphaDiversity")
  dir.create("BetaDiversity")
  dir.create("Others")

  message("Building phyloseq object, may take time depending on size of the file")
  ps0 <- read_phyloseq(otu.file = biom, metadata.file = mapping, taxonomy.file = taxonomy, type = "biom")
  print(ps0);
  treefile_p0 <- read.tree(treefilename)
  qc_plot1 <- plot_ReadDistribution(ps0, groups = VariableA,
                                    plot.type = "histogram")
  ggsave("./QC/ReadDistribution.pdf")

  # Investigate library sizes
  message("Investigating library sizes")
  pdf("./QC/LibrarySizePerSample.pdf")
  barplot(sort(sample_sums(ps0)), horiz = TRUE, las = 2, sub = "Check if there is large differeence in library sizes")
  hist(sort(sample_sums(ps0)), las = 2)
  dev.off()

  message("Investigating OTU counts distribution")
  pdf("./QC/Distribution_OTU_Counts.pdf")
  barplot(sort(taxa_sums(ps0)), horiz = TRUE, las = 2, sub = "Check if there is large differeence in otu counts")
  hist(taxa_sums(ps0), las = 2, main = "raw")
  dev.off()
  ps1 <- prune_taxa(taxa_sums(ps0) > 0, ps0)

  message("Checking for variance in OTU counts")
  # variance in data

  Variance.plot.a <- qplot(log10(apply(otu_table(ps1), 1, var)),
                           xlab = "log10(variance)", main = "Variance in OTUs") +
    ggtitle("before filtering") + theme_minimal()
  ggsave("./QC/Variance before filtering.pdf")
  if (filterpseq == TRUE)
  {
    message("Filtering OTUs with less than 10 counts in at least 5 % of the dataset")
    ps2 = filter_taxa(ps1, function(x) sum(x > filterCount) >
                        (filterPrev * length(x)), TRUE)
    message("below is the infor of filtered phyobject")
    print(ps2)
  } else
  {
    ps2 <- ps1
  }

  # Filtering to remove spurious counts
  Variance.plot.b <- qplot(log10(apply(otu_table(ps1), 1, var)),
                           xlab = "log10(variance)", main = "Variance in OTUs") +
    ggtitle("after filtering") + theme_minimal()

  ggsave("./QC/Variance After filtering.pdf")

  # Alpha diversities

  message("Calculating diversity indices")
  #otu_tab_ps2 <- abundances(ps2)
  #phy_tree_ps2 <- ps2@phy_tree

  #df.pd <- pd(t(otu_tab_ps2), otu_tab_ps2, include.root = F)

  metadf <- meta(ps2)

  #metatab <- cbind(metadf, df.pd)  #,df.pd)
  #write.csv("./AlphaDiversity/Alpha_diversity_metadata.csv")

  #alpha.plot <- ggboxplot(metatab, VariableA, "PD", fill = VariableB,
  # palette = "rickandmorty") + ggtitle("Phylogenetic diversity")
  #print(alpha.plot)
  #ggsave("./Figures/Phylogenetic_alpha_diversity.pdf")

  # All
  alpha_div <- plot_richness(ps2, color = VariableA, shape = VariableB,
                             measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
  alpha_div <- alpha_div + geom_boxplot()
  ggsave("./AlphaDiversity/Non-phylogenetic_alpha_diversity.pdf",
         height = 6, width = 18)

  # Beta diversity
  message("Analysing beta diversitiy")

  if (transformation != FALSE)
  {

    ps3 <- microbiome::transform(ps2, transform = transformation)
  } else
  {
    ps3 <- ps2
  }

  ord <- ordinate(ps3, distance = Distmethod, method = OrdinatMethod)

  ordi_plot <- plot_ordination(ps3, ord, color = VariableA,
                               shape = VariableB)
  ordi_plot <- ordi_plot + scale_color_brewer(palette = col.palette) +
    theme_bw() + geom_point(size = 2) + ggtitle(paste("./BetaDiversity/",
                                                      OrdinatMethod, " based on ", Distmethod))
  print(ordi_plot)

  ggsave(filename = paste("./BetaDiversity/", OrdinatMethod,
                          " based on ", OrdinatMethod, ".pdf", sep = ""), plot = ordi_plot,
         height = 8, width = 12)

  message("Plotting composition barplot")

  taxic <- as.data.frame(ps2@tax_table)
  colourCount = length(unique(taxic$Phylum))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
  getPalette = colorRampPalette(brewer.pal(8, col.palette))  # change the palette as well as the number of colors will change according to palette.
  guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15,
                                                                         face = "italic", colour = "Black", angle = 0)))
  comp.plot <- plot_composition(ps2, taxonomic.level = "Phylum",
                                transform = "compositional")
  comp.plot <- comp.plot + theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Relative abundance") + guide_italics #+
    #scale_fill_brewer(palette = getPalette(colorCount))

  if (nrow(metadf) > 30)
  {
    ggsave("./Others/compositionbarplot.pdf", plot = comp.plot,
           height = 8,
           width = 28)
  } else
  {
    ggsave("./Others/compositionbarplot.pdf", plot = comp.plot,
           height = 8,
           width = 18)
  }
  message("setting back to this" , setwd(oripath));
  setwd(oripath)
  print("done!")
}



#############plot_read_dist
#' @title Distribution of reads
#' @description Plots distribution of reads.
#' @param x \code{\link{phyloseq-class}} object
#' @param groups Metadata variable to check for groups based sequencing effort.
#' @param plot.type either density or histogram plot
#
#' @author Contact: SUdarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @return A \code{\link{ggplot}} plot object
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data(DynamicIBD)
#'     ps0 <- DynamicIBD
#'     p <- plot_ReadDistribution(ps0, groups="ibd_subtype", plot.type= "density")
#'           }
#' @keywords utilities

plot_ReadDistribution <-
  function(x,
           groups,
           plot.type = c("density", "histogram")) {
    df <- pdfa <- pdfb <- Reads_per_sample <- NULL
    title <- "Distribution of reads in the dataset"
    df <-
      data.table(
        as(sample_data(x), "data.frame"),
        Reads_per_sample = sample_sums(x),
        keep.rownames = TRUE
      )
    if (plot.type == "density") {
      p.dfa <-
        ggplot(df, aes(x = Reads_per_sample, fill = factor(groups))) + geom_density(alpha = 0.5, fill = "steelblue") + facet_wrap(groups) +
        theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
      print(p.dfa + ggtitle(title))
    } else if (plot.type == "histogram") {
      p.dfb <-
        ggplot(df, aes(x = Reads_per_sample, fill = factor(groups))) + geom_histogram(alpha = 0.5) + facet_wrap(groups) +
        theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
      print(p.dfb + ggtitle(title))
    }

    print("Done plotting")
  }
