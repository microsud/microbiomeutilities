#' @title Convert phyloseq object to long data format
#' @description Faster alternative to psmelt function from \code{\link{phyloseq-class}} object. It is important that the sample name in metadata should have "SampleID" and this should match the sample names in OTU table.
#' @param x \code{\link{phyloseq-class}} object
#
#' @return A data frame in long format.
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data("dietswap")
#'     pseq <- subset_samples(dietswap, group == "DI" & nationality == "AFR")
#'     pseq_df <- phy_to_ldf(pseq)
#'           }
#' @keywords utilities
phy_to_ldf = function(x)
{
  message("sample_data must have a column named SampleID, \n consisting of namesm corresponding to samples names in otu_table")
  meta_df = data.frame(sample_data(x))
  tax_df = data.frame(tax_table(x)) %>% 
    rownames_to_column("OTU")
  otu_df = data.frame(otu_table(x), 
                      check.names = FALSE) %>% rownames_to_column("OTU")
  suppressWarnings(otu_df %>% left_join(tax_df) %>% gather_("SampleID", 
                                           "Abundance", setdiff(colnames(otu_df), 
                                                                "OTU")) %>% left_join(meta_df))
}