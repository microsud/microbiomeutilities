#' @title List of available datasets
#' @description Data are used from Duvallet et al 2017 https://www.nature.com/articles/s41467-017-01973-8.pdf?origin=ppub.
#' @details Data for practice, also an example for importing mothur files from Baxtrer et al 2016. The source file for these data is the microbiomedatarepo https://github.com/microsud/microbiomedatarepo.
#' @export
#' @examples
#'     library(microbiomeutilities)
#'
#'     list_microbiome_data()
#'
#' @keywords utilities
#'

list_microbiome_data <- function(){

  # TODO add all the dataset from Duvallet et al 201 as phyloseq objects

  print("data from Duvallet et al 2017 https://www.nature.com/articles/s41467-017-01973-8.pdf?origin=ppub")
  message("1] Schubert et al., CDI 336 samples [cdi_schubert.rds]")
  message("2] Youngster et al., CDI 23 samples [cdi_youngster.rds]")
  message("3] Goodrich et al., OB 613 samples [ob_goodrich.rds]")
  message("4] Dinh et al., HIV 36 samples [hiv_dinh.rds]")
  message("5] Zeller et al., CRC 116 samples [crc_zeller.rds]")
  message("6] example of mothur file from Baxter et al., 2016 \n i) baxter_shared.shared \n ii) baxter_con.taxonomy \n iii) baxter_metadat.csv")

}
