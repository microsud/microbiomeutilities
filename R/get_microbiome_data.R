#' @title Download test microbiome data
#' @description Test microbiome data in phyloseq format.
#' @details You can download few example datasets in phyloseq format from Duvallet et al 2017 https://www.nature.com/articles/s41467-017-01973-8.pdf?origin=ppub. The source file for these data is the microbiomedatarepo https://github.com/microsud/microbiomedatarepo.
#' @param studynumber A numerical value corresponsing to the one given by list_microbiome_data()
#' @return \code{\link{phyloseq-class}} object.
#' @import microbiome
#' @import phyloseq
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     library(microbiomeUtilities)
#'     list_microbiome_data()
#'     ps1 <- get_microbiome_data(studynumber = 1)
#'     # this will download the first study listed in list_microbiome_data() function
#'     print(p)
#'           }
#' @keywords utilities
get_microbiome_data <- function(studynumber = 1) {

  if (studynumber == 1) {

    ps1 <- readRDS(url("https://github.com/microsud/microbiomedatarepo/blob/master/datasets/microbiomeDB/cdi_schubert.rds?raw=true"))

  } else if (studynumber == 2) {

    ps1 <- readRDS(url("https://github.com/microsud/microbiomedatarepo/blob/master/datasets/microbiomeDB/cdi_youngster.rds?raw=true"))

  } else if (studynumber == 3) {

    ps1 <- readRDS(url("https://github.com/microsud/microbiomedatarepo/blob/master/datasets/microbiomeDB/ob_goodrich.rds?raw=true"))

  } else if (studynumber == 4) {

    ps1 <- readRDS(url("https://github.com/microsud/microbiomedatarepo/blob/master/datasets/microbiomeDB/hiv_dinh.rds?raw=true"))

  } else if (studynumber == 5) {

    ps1 <- readRDS(url("https://github.com/microsud/microbiomedatarepo/blob/master/datasets/microbiomeDB/crc_zeller.rds?raw=true"))

  } else if (studynumber == 6) {

    print("Use the following link and use the read_phyloseq function from microbiome package")
    message("consensus taxonomy file link https://raw.githubusercontent.com/microsud/microbiomedatarepo/master/datasets/Baxter_FITs_Microbiome_2016/baxter_con.taxonomy")
    message("shared otu file link https://raw.githubusercontent.com/microsud/microbiomedatarepo/master/datasets/Baxter_FITs_Microbiome_2016/baxter_shared.shared")
    message("metadata file link https://raw.githubusercontent.com/microsud/microbiomedatarepo/master/datasets/Baxter_FITs_Microbiome_2016/baxter_metadat.csv")

  }

}
