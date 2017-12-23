#' @title Microbiome analysis pipeline
#' @param 
#' @param otufile biom object, otu_table in csv format or mothur shared files.   
#' @param mapping Metadata variable to check for groups based sequencing effort csv format.              
#' @param taxonomy NULL or csv fomatted file.
#' @param treefilename for phylogenetic based diversity analysis "*.tre" file 
#' @param type  "biom", "mothur", "simple" simple is for *.csv file.
#' @param work_dir Working directory where input files are stored.
#' @param out_dir Output directory where all outputs are to be stored.
#' @param VariableA Main variable of interest.
#' @param VariableB Secondary variable of interest.           
#' @param UnConstOrd If unconstrained ordination to be plotted TRUE or FALSE.
#' @param heatmap Heatmap or not option is TRUE or FALSE.
#' @param filterCount Filter OTUs below this count number. 
#' @param filterPrev Filter OTUs not detected in less than this percent of samples.
#' @param col.palette = 'Spectral', 'Set2', 'Paired' use any of the RColorbrewer palette depending on number of groups 
#'                      you have in VaraibleA.
#' @param filterpseq TRUE or FALSE.
#' @param samsize Number of reads to rarefy your data and save the rarefied phyobject.
#' @param projectname Name of the project.
#' @param author Name of the author/investigator.
#' @author Contact: Sudarshan Shetty \email{sudarshanshetty9@@gmail.com}
#' @return A HTML report with graphs and data stats.
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     data(DynamicIBD)
#'     ps0 <- DynamicIBD
#'     library(microbiomeutilities)
#'     microbiome_pipeline(otufile = "NGTaxMerged_conv.biom",  
#'                    mapping = "mappingMerged_edit.csv",              
#'                    taxonomy = NULL,                    
#'                    treefilename = "combinedTree.tre",            
#'                    type = "biom",
#'                    work_dir = "F:/R_studio_pacakes/MICROBIOMEUTILITIES/microbiomeutilities-master/testdata",
#'                    out_dir = "F:/R_studio_pacakes/MICROBIOMEUTILITIES/microbiomeutilities-master/testdata/myoutput", 
#'                    VariableA = "MC_type1",               
#'                    VariableB = "Region",           
#'                    UnConstOrd = TRUE,                           
#'                    heatmap = TRUE,                          
#'                    filterCount = 4,                           
#'                    filterPrev = 0.01,                          
#'                    col.palette = "Paired",                    
#'                    filterpseq = TRUE,                          
#'                    samsize = NA,
#'                    projectname = "Mock_communities",
#'                    author = "Sudarshan")
#' @keywords utilities

microbiome_pipeline <- function(otufile,
                                mapping,
                                taxonomy,
                                treefilename,
                                type,
                                work_dir,
                                out_dir,
                                VariableA,
                                VariableB,
                                UnConstOrd,
                                heatmap,
                                filterCount,
                                filterPrev,
                                col.palette,
                                filterpseq,
                                samsize,
                                projectname,
                                author){
  {
    rmarkdown::render(input=system.file("microutility_template.rmd", package="microbiomeutilities"),
                      output_file=paste0(projectname, "_report.html"),
                      output_dir=work_dir,
                      knit_root_dir=work_dir,
                      run_pandoc=TRUE,
                      quiet=FALSE,
                      clean=TRUE,
                      params = list(
                        otufile,
                        mapping,
                        taxonomy,
                        treefilename,
                        type,
                        work_dir,
                        out_dir,
                        VariableA,
                        VariableB,
                        UnConstOrd,
                        heatmap,
                        filterCount,
                        filterPrev,
                        col.palette,
                        filterpseq,
                        samsize,
                        projectname,
                        author))
    cat("HTML report created in output directory\n")
  }
  
}                               
                                    

