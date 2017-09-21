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
