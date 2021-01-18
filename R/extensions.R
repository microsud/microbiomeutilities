#' @title Convert Phyloseq Slots to Tibbles
#' @description Utility to convert phyloseq slots to tibbles.
#' @param x \code{\link{phyloseq-class}} object
#' @param slot Must be one of c("otu_table", "sam_data", "tax_table").
#'        Default= "otu_table"   
#' @param column_id Provide name for the column which will hold the rownames of slot.
#' @return A tibble
#' @examples
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' otu_tibble <- get_tibble(p0,slot="otu_table",column_id="OTUID")
#' head(otu_tibble)
#' @export
#' @keywords utilities

get_tibble <- function(x, 
                       slot="otu_table",
                       column_id="column_id") {
  
  if (class(x)!="phyloseq"){
    stop("Input is not an object of phyloseq class")
  }
  
  if(isFALSE(any(slot==getslots.phyloseq(x)))) {
   stop("slot must be one of 'otu_table', 'sam_data', 'tax_table'")
  }
  
  if(slot=="otu_table"){
    
    tib_dat <- abundances(x) %>% 
      as.data.frame(stringsAsFactors=FALSE) %>% 
      rownames_to_column(column_id) %>% 
      as_tibble()
    
    return(tib_dat)
    
  } else if(slot=="sam_data"){
    tib_dat <- meta(x) %>% 
      as.data.frame(stringsAsFactors=FALSE) %>% 
      rownames_to_column(column_id) %>% 
      as_tibble()
    
    return(tib_dat)
    
  } else if(slot=="tax_table"){
    tib_dat <- tax_table(x) %>% 
      as("matrix") %>% 
      as.data.frame(stringsAsFactors=FALSE) %>% 
      rownames_to_column(column_id) %>% 
      as_tibble()
    
    return(tib_dat)
    
  }
  
  
}


#' @title Join otu_table and tax_table to Tibble
#' @description Utility to join otu_table and tax_table to tibble.
#' @param x \code{\link{phyloseq-class}} object
#' @param column_id Provide name for the column which will hold the rownames of slot.
#' @return A tibble
#' @examples
#' library(microbiomeutilities)
#' data("zackular2014")
#' p0 <- zackular2014
#' otu_tax <- join_otu_tax(p0,column_id = "OTUID")
#' head(otu_tax)
#' @export
#' @keywords utilities

join_otu_tax <- function(x, column_id = "OTUID"){
  
  otu_tb <- get_tibble(x, slot="otu_table", column_id) 
  tax_tb <- get_tibble(x, slot="tax_table", column_id)
  otu_tb <- tax_tb %>% 
    left_join(otu_tb, by=column_id)
  
}


#' @title Add refseq slot for dada2 based phyloseq object
#' @description Utility to add refseq slot for dada2 based phyloseq object. Here, the
#' rownames which are unique sequences, are stored in refseq slot of phyloseq. Sequence 
#' ids are converted to ids using tag option.  
#' @param x \code{\link{phyloseq-class}} object with seqs as rownames.
#' @param tag Provide name for Ids, Default="ASV".
#' @return \code{\link{phyloseq-class}} object 
#' @examples
#' 
#' # ps <- add_refseq(p0,tag="ASV")
#' # ps
#' 
#' @export
#' @keywords utilities
#' @importFrom Biostrings DNAStringSet
add_refseq <- function(x, tag="ASV"){
  
  if (class(x)!="phyloseq"){
    stop("Input is not an object of phyloseq class")
  }
  
  nucl <- Biostrings::DNAStringSet(taxa_names(x))
  names(nucl) <- taxa_names(x)
  x <- merge_phyloseq(x, nucl)
  
  rm(nucl)
  
  if(is.na(tag) || is.null(tag)){
    taxa_names(x) <- paste0("taxa", seq(ntaxa(x)))
    return(x)
  } else{
    taxa_names(x) <- paste0(tag, seq(ntaxa(x)))
    return(x)
  }
  
}



