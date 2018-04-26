#' @title List of available datasets
#' @description Data are used from Duvallet et al 2017 https://www.nature.com/articles/s41467-017-01973-8.pdf?origin=ppub.
#' @details Data for practice, also an example for importing mothur files from Baxtrer et al 2016. The source file for these data is the microbiomedatarepo https://github.com/microsud/microbiomedatarepo.
#' @param printtab Print in console or not, defaut is TRUE and will print output.
#' @export
#' @examples
#'     library(microbiomeutilities)
#'
#'     df <- list_microbiome_data(printtab = FALSE)
#'
#' @keywords utilities
#'

list_microbiome_data <- function(printtab = TRUE){
  
  # TODO add all the dataset from Duvallet et al 201 as phyloseq objects
  Study <-
    c(
      "Son2015_ASD",
      "Kang2013_ASD",
      "Schubert2014_CDI",
      "Youngster2014_CDI",
      "Baxter2016_CRC",
      "Zackular2014_CRC",
      "Zeller2014_CRC",
      "Singh2015_EDD",
      "NogueraJulian2016_HIV",
      "Dinh2015_HIV",
      "Lozupone2013_HIV",
      "Gevers2014_IBD",
      "Zhang2013_LIV",
      "Wong2013_NASH",
      "Ross2015_OB",
      "Zupancic2012_OB",
      "Scher2013_PAR",
      "Alkanani2015_T1D",
      "Scheperjans2015_PAR",
      "Alkanani2015_T1D"
    )
  
  Disease <-
    c(
      "ASD",
      "ASD",
      "CDI",
      "CDI",
      "CRC",
      "CRC",
      "CRC",
      "EDD",
      "HIV",
      "HIV",
      "HIV",
      "IBD",
      "LIV",
      "NASH",
      "OB",
      "OB",
      "PAR",
      "T1D",
      "PAR",
      "T1D"
    )
  
  microbiomeDB_pseq <- as.data.frame(cbind(Study, Disease))
  if(printtab == TRUE){
    print(microbiomeDB_pseq)
  } else {
    
    return(microbiomeDB_pseq)
  }
  
  
}
