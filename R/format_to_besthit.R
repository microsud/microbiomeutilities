#' @title Formatting the Phyloseq Object advanced
#' @description Format the phyloseq object to add the best taxonomy in phyloseq object (tax_table and otu_table).
#' @details Most commonly it is observed that the taxonomy file has classification until a given taxonomic level.
#'          row.names for both tax_table and otu_table have best hit, until maximun genus level (species classification with short amplicons is a myth)is made available. This code is a 
#'          slight modification of the code from  \pkg{ampvis} \code{\link{phyloseq-class}}. 
#'          Here, we directly take the phyloseq object as input and make the necessary formatting.
#' @param x \code{\link{phyloseq-class}} object
#' @return  \code{\link{phyloseq-class}} object.
#' @export
#' @examples \dontrun{
#'   # Example data
#'     library(microbiome)
#'     library(microbiomeUtilities)
#'     data(DynamicsIBD)
#'     p0 <- DynamicsIBD
#'     p0.f <- format_to_besthit(p0)
#'           }
#' @keywords utilities

format_to_besthit <- function(pobj){
  
  Domain <- Phylum <- Class <- Order <- Family <- Genus <- Species <- x <- y <- NULL
  
  x <- pobj
  
  if (ncol(tax_table(x)) == 6){
    colnames(tax_table(x)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus" )
  } else if (ncol(tax_table(x)) == 7){
    colnames(tax_table(x)) <- c("Domain", "Phylum", "Class", "Order",  "Family", "Genus", "Species" )
  } else {
    stop("Taxonomic levels should be either 6 (untill genus) or 7 (until species) level")
  }
  
  
  tax_table(x)[,1][is.na(tax_table(x)[,1])] <- paste0(tolower(substring("kingdom", 1, 1)), "__")
  
  tax_table(x)[,2][is.na(tax_table(x)[,2])] <- paste0(tolower(substring("Phylum", 1, 1)), "__")
  
  tax_table(x)[,3][is.na(tax_table(x)[,3])] <- paste0(tolower(substring("Class", 1, 1)), "__")
  
  tax_table(x)[,4][is.na(tax_table(x)[,4])] <- paste0(tolower(substring("Order", 1, 1)), "__")
  
  tax_table(x)[,5][is.na(tax_table(x)[,5])] <- paste0(tolower(substring("Family", 1, 1)), "__")
  
  tax_table(x)[,6][is.na(tax_table(x)[,6])] <- paste0(tolower(substring("Genus", 1, 1)), "__")
  
  if (ncol(tax_table(x)) == 7){
    tax_table(x)[,7][is.na(tax_table(x)[,7])] <- paste0(tolower(substring("Species", 1, 1)), "__")
  } 
  
  
  y <- as.data.frame.matrix(tax_table(x))
  
  
  y$Domain <- gsub("k__", "", y$Domain)
  y$Phylum <- gsub("p__", "", y$Phylum)
  y$Class <- gsub("c__", "", y$Class)
  y$Order <- gsub("o__", "", y$Order)
  y$Family <- gsub("f__", "", y$Family)
  y$Genus <- gsub("g__", "", y$Family)
  
  if (ncol(tax_table(y)) == 7){
    y$Species <- gsub("s__", "", y$Species)
  }
  
  if (ncol(tax_table(x)) == 6){
    tax <- mutate(y, Domain, Domain = ifelse(Domain == "", "Unclassified", Domain)) %>% mutate(Phylum,
                                                                                               Phylum = ifelse(Phylum == "", paste("k__", Domain, "", sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__",
                                                                                          Phylum, "", sep = "")), Class)) %>% mutate(Order, Order = ifelse(Order ==
                                                                                                                                                                           "", ifelse(grepl("__", Class), Class, paste("c__", Class, "", sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__",
                                                                                           Order, "", sep = "")), Family)) %>% mutate(Genus, Genus = ifelse(Genus ==
                                                                                                                                                                            "", ifelse(grepl("__", Family), Family, paste("f__", Family, "", sep = "")),
                                                                                                                                                                          Genus))
    
  } else if (ncol(tax_table(x)) == 7){
    
    tax <- mutate(y, Domain, Domain = ifelse(Domain == "", "Unclassified", Domain)) %>% mutate(Phylum,
                                                                                               Phylum = ifelse(Phylum == "", paste("k__", Domain, "", sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__",
                                                                                          Phylum, "", sep = "")), Class)) %>% mutate(Order, Order = ifelse(Order ==
                                                                                                                                                                           "", ifelse(grepl("__", Class), Class, paste("c__", Class, "", sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__",
                                                                                           Order, "", sep = "")), Family)) %>% mutate(Genus, Genus = ifelse(Genus ==
                                                                                                                                                                            "", ifelse(grepl("__", Family), Family, paste("f__", Family, "", sep = "")),
                                                                                                                                                                          Genus)) %>% mutate(Species, Species = ifelse(Species == "", ifelse(grepl("__", Genus), Genus,
                                                                                                                                                                                                                                             paste("g__", Genus, "", sep = "")), Species))
    
  } 
  
  me <- as.data.frame(x@tax_table)
  me$Domain <- tax$Domain
  me$Phylum <- tax$Phylum
  me$Class <- tax$Class
  me$Order <- tax$Order
  me$Family <- tax$Family
  me$Genus <- tax$Genus
  
  if (ncol(y) == 7){
    me$Species <- tax$Species
  }
  
  taxmat <- as.matrix(me)
  new.tax <- tax_table(taxmat)
  tax_table(x) <- new.tax
  
  ### Now add the best hist as rownames and final column in tax_table
  taxa_names(x) <- paste("OTU-", taxa_names(x), sep="")
  tax.g <- as.data.frame.matrix(tax_table(x))
  tax.g$col1 <- tax.g$Genus
  tax.g$col2 <- rownames(tax.g)
  tax.unit <- tidyr::unite(tax.g, best_hit,c("col2", "col1"), sep = ":", remove = TRUE)
  rownames(tax.unit) <- tax.unit$best_hit
  otu.1 <- as.data.frame.matrix(abundances(x))
  rownames(otu.1) <- tax.unit$best_hit
  
  tax.new <- as.data.frame.matrix(tax_table(x))
  tax.new$top_hit <- tax.unit$best_hit
  rownames(tax.new) <- tax.new$top_hit
  
  OTU = otu_table(as.matrix(otu.1), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(tax.new))
  sampledata <- sample_data(meta(x))
  p.new <- merge_phyloseq(OTU, TAX, sampledata)
  
  return(p.new)
}


