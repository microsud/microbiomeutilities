#' @title Download test microbiome data
#' @description Test microbiome data in phyloseq format.
#' @details You can download few example datasets in phyloseq format from Duvallet et al 2017 https://www.nature.com/articles/s41467-017-01973-8.pdf?origin=ppub. The source file for these data is the microbiomedatarepo https://github.com/microsud/microbiomedatarepo.
#' @param disease Disease of interest as shown in list_microbiome_data()
#' @param study Name of the study as shown in list_microbiome_data()
#' @return \code{\link{phyloseq-class}} object.
#' @import microbiome
#' @import phyloseq
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     library(microbiomeUtilities)
#'     list_microbiome_data()
#'     ps1 <- get_microbiome_data(disease = "CDI", "Schubert2014_CDI")
#'     print_ps(ps1)
#'           }
#' @keywords utilities
get_microbiome_data <- function(disease, study) {
  
  
  if(grepl(disease, study) == FALSE){
    
    stop("Disease-Study combination does not exists, check list_microbiome()")
    
  }
  
  fileloc <-
    paste0(
      "https://github.com/microsud/microbiomedatarepo/blob/master/datasets/microbiomeDB/",
      disease,
      "/",
      study,
      ".rds?raw=true"
    )
  
  psX <- readRDS(url(fileloc))
  
  return(psX)
}



