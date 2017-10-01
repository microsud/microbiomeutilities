#' @title Microbiome analysis pipeline
#' @description This is a single workflow for prelimanry investigation of microbiome data.
#' @param biom biom object object
#' @param mapping Metadata variable to check for groups based sequencing effort.
#' @param taxonomy either density or histogram plot.
#' @param treefilename tree .
#' @param type "biom", "mothur", "simple" simple is for *.csv file.
#' @param out_dir otuput directory
#' @param VariableA Main variable of interest.
#' @param VariableB Secondary variable of interest.
#' @param filterpseq TRUE or FALSE.
#' @param Distmethod Distance method 'bray', 'unifrac', 'wunifrac', 'jsd'.
#' @param OrdinatMethod Ordination method 'NMDS', 'MDS', 'PCoA'.
#' @param QC Default TRUE QC plots.
#' @param filterCount Filter OTUs below this count number.
#' @param filterPrev Filter OTUs not detected in less than this percent of samples
#' @param transformation  'compositional', 'Z', 'log10', 'hellinger', 'clr' choose to tranform your
#'        counts data defaults to FALSE
#' @param col.palette 'Spectral', 'Set2', 'Paired' use any of the RColorbrewer palette depending on number of groups #'                     you have in VaraibleA
#' @param samsize Number of reads to rarefy your data and save the rarefied phyobject.
#' @author Contact: SUdarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @return A \code{\link{ggplot}} plot object
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data(DynamicIBD)
#'     ps0 <- DynamicIBD
#'     run_microbiome_pipeline(biom = "myBiom.biom",
#'                       mapping = "myMapping.csv",
#'                      taxonomy = NULL,
#'                      treefilename = "myTree.tre",
#'                      type = "biom",
#'                      out_dir = "F:/PATH/to/my/output/directory",
#'                      VariableA = "MY_MainVariable",
#'                      VariableB = "MY_SecondaryVariable",
#'                      Distmethod = "bray",
#'                      OrdinatMethod = "MDS",
#'                      QC = TRUE,
#'                      filterCount = 4,
#'                      filterPrev = 0.01,
#'                      transformation = FALSE,
#'                      col.palette = "Set2",
#'                      filterpseq = TRUE, samsize=0)
#' @keywords utilities


run_microbiome_pipeline <- function(biom, mapping, 
                                    taxonomy, treefilename,
                                    type, out_dir, 
                                    VariableA, VariableB, 
                                    Distmethod, OrdinatMethod, 
                                    QC, filterCount, 
                                    filterPrev, transformation, 
                                    col.palette, filterpseq, 
                                    samsize) {

  ps0 <- qc_plot <- ps1 <- ps2 <- Variance.plot.a <- Variance.plot.b <- otu_tab_ps2 <- phy_tree_ps2 <- df.pd <- metadf <- metatab <- tab <- alpha_div <- alpha.plot <- ord <- NULL
  oripath <- getwd()
  setwd(out_dir)
  dir.create("QC")
  dir.create("AlphaDiversity")
  dir.create("BetaDiversity")
  dir.create("Others")
  dir.create("PhyloseqObjects")

  message("Building phyloseq object, may take time depending on size of the file")
  if (type == "biom") {
    ps0 <- read_phyloseq(otu.file = biom, metadata.file = mapping,
                         taxonomy.file = taxonomy, type = "biom")
    message("using biom as input")
    print(ps0)
  } else if (type == "simple") {
    ps0 <- read_phyloseq(otu.file = biom, metadata.file = mapping,
                         taxonomy.file = taxonomy, type = "simple")
    message("using simple csv format as input")
    print(ps0)
  } else if (type == "mothur") {
    ps0 <- read_phyloseq(otu.file = biom, metadata.file = mapping,
                         taxonomy.file = taxonomy, type = "mothur")
    message("using mothur format as input")
    print(ps0)
  } else {
    stop("Unrecognized type in read_phyloseq input. Exiting.")
  }

  print(ps0)

  # Add tree file is supplied by user

  if (!is.null(treefilename)){

    treefile_p0 <- read.tree(treefilename)

    ps0 <- merge_phyloseq(ps0, treefile_p0)

    saveRDS(ps0, "./PhyloseqObjects/ps_raw.rds")

    message("Below is the content of raw phyloseqobject stored as ps_raw.rds")
    print(ps0)
  } else {

    message("No tree file supplied hence phyloseq phy_tree is NULL");
    print(ps0)

  }


  qc_plot1 <- plot_ReadDistribution(ps0, groups = VariableA, plot.type = "histogram")
  ggsave("./QC/ReadDistribution.pdf")

  # Investigate library sizes
  message("Investigating library sizes")
  pdf("./QC/LibrarySizePerSample.pdf")
  barplot(sort(sample_sums(ps0)), horiz = TRUE, las = 2, sub = "Check if there is large differeence in library sizes")
  hist(sort(sample_sums(ps0)), las = 2)
  dev.off()

  barplot(sort(sample_sums(ps0)), horiz = TRUE, las = 2, sub = "Check if there is large differeence in library sizes")
  hist(sort(sample_sums(ps0)), las = 2)

  message("Investigating OTU counts distribution")
  pdf("./QC/Distribution_OTU_Counts.pdf")
  barplot(sort(taxa_sums(ps0)), horiz = TRUE, las = 2, sub = "Check if there is large differeence in otu counts")
  hist(taxa_sums(ps0), las = 2, main = "raw")
  dev.off()

  barplot(sort(taxa_sums(ps0)), horiz = TRUE, las = 2, sub = "Check if there is large differeence in otu counts")
  hist(taxa_sums(ps0), las = 2, main = "raw")

  ps1 <- prune_taxa(taxa_sums(ps0) > 0, ps0)

  message("Checking for variance in OTU counts")
  # variance in data

  Variance.plot.a <- qplot(log10(apply(otu_table(ps1), 1, var)), xlab = "log10(variance)",
                           main = "Variance in OTUs") + ggtitle("before filtering") + theme_minimal()
  print(Variance.plot.a)
  ggsave("./QC/Variance before filtering.pdf")

  if (filterpseq == TRUE) {
    message("Filtering OTUs with less than 10 counts in at least 5 % of the dataset")
    ps2 = filter_taxa(ps1, function(x) sum(x > filterCount) > (filterPrev *
                                                                 length(x)), TRUE)
    message("Saving the filtered phyobject"
            "Saving the transformed phyloseq object as ps_filtered.rds")

    saveRDS(ps2, "./PhyloseqObjects/ps_filtered.rds")

    message("Below is the content of filtered phyloseqobject (based on filterCount and filterPrev) stored as ps_filtered.rds")

    print(ps2)

  } else {
    message("filterpseq was false. Did not filter and hence, will not save the filtered phyloseq")

    ps2 <- ps1

  }

  # Filtering to remove spurious counts
  Variance.plot.b <- qplot(log10(apply(otu_table(ps1), 1, var)), xlab = "log10(variance)",
                           main = "Variance in OTUs") + ggtitle("after filtering") + theme_minimal()

  print(Variance.plot.b)

  ggsave("./QC/Variance After filtering.pdf")

  # Alpha diversities

  message("Calculating diversity indices")
  # otu_tab_ps2 <- abundances(ps2) phy_tree_ps2 <- ps2@phy_tree

  # df.pd <- pd(t(otu_tab_ps2), otu_tab_ps2, include.root = F)

  metadf <- meta(ps2)

  # metatab <- cbind(metadf, df.pd) #,df.pd)
  # write.csv('./AlphaDiversity/Alpha_diversity_metadata.csv')

  # alpha.plot <- ggboxplot(metatab, VariableA, 'PD', fill = VariableB,
  # palette = 'rickandmorty') + ggtitle('Phylogenetic diversity')
  # print(alpha.plot)
  # ggsave('./Figures/Phylogenetic_alpha_diversity.pdf')

  # All
  alpha_div <- plot_richness(ps2, color = VariableA, shape = VariableB,
                             measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
  alpha_div <- alpha_div + geom_boxplot()

  print(alpha_div)

  ggsave("./AlphaDiversity/Non-phylogenetic_alpha_diversity.pdf", height = 6,
         width = 18)

  # Beta diversity
  message("Analysing beta diversitiy")

  if (transformation != FALSE) {
    message("#Transformation was selected and will transform your phyloseq object"
            "#Saving the transformed phyloseq object as ps_transformed.rds")

    ps3 <- microbiome::transform(ps2, transform = transformation)

    saveRDS(ps3, "./PhyloseqObjects/ps_transformed.rds")
    message("Below is the content of transfored phyloseqobject based on transformation = TRUE stored as ps_transformed.rds")

    print(ps3)

  } else {
    ps3 <- ps2
    message("Below is the content of non-transfored phyloseqobject based on transformation = FALSE")

    print(ps3)
  }

  ord <- ordinate(ps3, distance = Distmethod, method = OrdinatMethod)

  ordi_plot <- plot_ordination(ps3, ord, color = VariableA, shape = VariableB)
  ordi_plot <- ordi_plot + scale_color_brewer(palette = col.palette) +
    theme_bw() + geom_point(size = 2) + ggtitle(paste("BetaDiversity",
                                                      OrdinatMethod, " based on ", Distmethod))
  print(ordi_plot)

  ggsave(filename = paste("./BetaDiversity/", OrdinatMethod, " based on ",
                          OrdinatMethod, ".pdf", sep = ""), plot = ordi_plot, height = 8,
         width = 12)

  message("Plotting composition barplot")

  taxic <- as.data.frame(ps2@tax_table)
  colourCount = length(unique(taxic$Phylum))  #define number of variable colors based on number of Family (change the level accordingly to phylum/class/order)
  getPalette = colorRampPalette(brewer.pal(8, col.palette))  # change the palette as well as the number of colors will change according to palette.
  guide_italics <- guides(fill = guide_legend(label.theme = element_text(size = 15,
                                                                         face = "italic", colour = "Black", angle = 0)))

  comp.plot <- plot_taxa_composition(ps2, taxonomic.level = "Phylum", transform = "compositional")

  comp.plot <- comp.plot + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Relative abundance") + guide_italics  #+
  # scale_fill_brewer(palette = getPalette(colorCount))

  if (nrow(metadf) > 30) {
    ggsave("./Others/compositionbarplot.pdf", plot = comp.plot, height = 8,
           width = 28)
  } else {
    ggsave("./Others/compositionbarplot.pdf", plot = comp.plot, height = 8,
           width = 18)
  }
  message("setting back to this", setwd(oripath))
  setwd(oripath)

  # Check for kurtosis

  message("Using the raw phyloseq to check for kurtosis in library size")
  require(data.table)

  df <- data.table(NumberReads = sample_sums(ps0), SampleID = sample_names(ps0))

  require(moments)

  n <- kurtosis(df$NumberReads)

  if (samsize != 0) {
    ps.1 <- rarefy_even_depth(ps, sample.size = samsize)
    saveRDS(ps.1, "./phyloseqObjects/ps_rarefyied.rds")
  }

  # check for kurtosis

  if (n > 3) {
    message("Your library size is heavily tailed, considering normalising them for further analysis")

  } else {
    message("Not rarefying since the variation in seq depth is not high")
    ps1 <- ps
  }
  print("Done with running the command!")
  message("Check the otuputs")
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
