

format_phylo_tax <- function(pobj) {
  Domain <- Phylum <- Class <- Order <- Family <- Genus <- Species <- x <- y <- NULL

  x <- pobj

  if (ncol(tax_table(x)) == 6) {
    colnames(tax_table(x)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  } else if (ncol(tax_table(x)) == 7) {
    colnames(tax_table(x)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  } else {
    stop("Taxonomic levels should be either 6 (untill genus) or 7 (until species) level")
  }


  tax_table(x)[, 1][is.na(tax_table(x)[, 1])] <- paste0(tolower(substring("kingdom", 1, 1)), "__")

  tax_table(x)[, 2][is.na(tax_table(x)[, 2])] <- paste0(tolower(substring("Phylum", 1, 1)), "__")

  tax_table(x)[, 3][is.na(tax_table(x)[, 3])] <- paste0(tolower(substring("Class", 1, 1)), "__")

  tax_table(x)[, 4][is.na(tax_table(x)[, 4])] <- paste0(tolower(substring("Order", 1, 1)), "__")

  tax_table(x)[, 5][is.na(tax_table(x)[, 5])] <- paste0(tolower(substring("Family", 1, 1)), "__")

  tax_table(x)[, 6][is.na(tax_table(x)[, 6])] <- paste0(tolower(substring("Genus", 1, 1)), "__")

  if (ncol(tax_table(x)) == 7) {
    tax_table(x)[, 7][is.na(tax_table(x)[, 7])] <- paste0(tolower(substring("Species", 1, 1)), "__")
  }


  y <- as.data.frame.matrix(tax_table(x))


  y[, 1] <- gsub("k__", "", y[, 1])
  y[, 2] <- gsub("p__", "", y[, 2])
  y[, 3] <- gsub("c__", "", y[, 3])
  y[, 4] <- gsub("o__", "", y[, 4])
  y[, 5] <- gsub("f__", "", y[, 5])
  y[, 6] <- gsub("g__", "", y[, 6])

  if (ncol(tax_table(y)) == 7) {
    y[, 7] <- gsub("s__", "", y[, 7])
  }

  if (ncol(tax_table(x)) == 6) {
    tax <- mutate(y, Domain, Domain = ifelse(Domain == "", "Unclassified", Domain)) %>%
      mutate(Phylum,
        Phylum = ifelse(Phylum == "", paste("k__", Domain, "_", rownames(y), sep = ""), Phylum)
      ) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__",
        Phylum, "_", rownames(y),
        sep = ""
      )), Class)) %>%
      mutate(Order, Order = ifelse(Order ==
        "", ifelse(grepl("__", Class), Class, paste("c__", Class, "_", rownames(y), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__",
        Order, "_", rownames(y),
        sep = ""
      )), Family)) %>%
      mutate(Genus, Genus = ifelse(Genus ==
        "", ifelse(grepl("__", Family), Family, paste("f__", Family, "_", rownames(y), sep = "")),
      Genus
      ))
  } else if (ncol(tax_table(x)) == 7) {
    tax <- mutate(y, Domain, Domain = ifelse(Domain == "", "Unclassified", Domain)) %>%
      mutate(Phylum,
        Phylum = ifelse(Phylum == "", paste("k__", Domain, "_", rownames(y), sep = ""), Phylum)
      ) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste("c__",
        Phylum, "_", rownames(y),
        sep = ""
      )), Class)) %>%
      mutate(Order, Order = ifelse(Order ==
        "", ifelse(grepl("__", Class), Class, paste("c__", Class, "_", rownames(y), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste("o__",
        Order, "_", rownames(y),
        sep = ""
      )), Family)) %>%
      mutate(Genus, Genus = ifelse(Genus ==
        "", ifelse(grepl("__", Family), Family, paste("f__", Family, "_", rownames(y), sep = "")),
      Genus
      )) %>%
      mutate(Species, Species = ifelse(Species == "", ifelse(grepl("__", Genus), Genus,
        paste("g__", Genus, "_", rownames(y), sep = "")
      ), Species))
  }

  me <- as.data.frame(x@tax_table)
  me$domain <- tax$Domain
  me$Phylum <- tax$Phylum
  me$Class <- tax$Class
  me$Order <- tax$Order
  me$Family <- tax$Family
  me$Genus <- tax$Genus

  if (ncol(y) == 7) {
    me$Species <- tax$Species
  }

  taxmat <- as.matrix(me)
  new.tax <- tax_table(taxmat)
  tax_table(x) <- new.tax
  return(x)
}
