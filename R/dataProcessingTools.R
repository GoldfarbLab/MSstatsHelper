#' Makes column names tidy
#'
#' @param df tibble or data frame
#' @return object of same class with tidy column names
#'
#' @export
#'
#' @import dplyr
formatColNames <- function(df) {
  cols <- colnames(df)
  cols <- str_replace_all(cols, ' ', '.')
  cols <- str_replace_all(cols, '\\.\\.', '\\.')
  cols <- str_replace_all(cols, '\\(', '')
  cols <- str_replace_all(cols, '\\)', '')
  colnames(df) <- cols

  return(df)
}


#' Given the protein name column, returns corresponding genes listed with a ';'
#'
#' @param x a column containing uniprot protein names listed by ';'
#' @param type a string specifying if the genes should be unique or not. If unique, specify "unique". If nonunique, do not include.
#' @return column of genes corresponding to input proteins. Listed with a ';'
#' @export
#'
#' @import dplyr
mapToGene <- function(x, type = ""){
  map <- invisible(map)

  hash <- hash(keys = map$Entry, values = map$gene.name)

  pboptions(type = "timer", char = "=", style = 3)
  mapped_genes <- pbsapply(str_split(x, ";"), function(protein_ids) ## pbsapply instead of sapply so that there is a progress bar
  { ### creates list of protein names

    genes <- sapply(str_match(protein_ids, "^[^-]*"), function(protein_id)
    { # individual protein level

      mapped_protein <- hash[[str_match(protein_id,"^[^_]*")]] # searches hash table
      mapped_protein <- mapped_protein[mapped_protein != "" & mapped_protein != " " &
                                         mapped_protein != "character(0)" & mapped_protein != "NULL"]

      return(as.character(mapped_protein)) # returns gene
    })

    if(type == "unique")
    { # for use in the Gene.names column, not for appending PTMs to name
      genes <- unique(genes)
    }

    genes <- genes[genes != "" & genes != " " & genes != "character(0)" & genes != "NULL"] # cleaning out non-matches

    return(str_c(genes, collapse = ";", sep = NULL)) # returns list of genes for one row of proteins
  })

  mapped_genes <- mapped_genes %>% replace_na("")
  return(trimws(mapped_genes))
}

#' Given the protein name column, returns corresponding organisms listed with a ';'
#'
#' @param x a column containing uniprot protein names listed by ';'
#' @param type a string specifying if the genes should be unique or not. If unique, specify "unique". If nonunique, do not include.
#' @return column of organisms corresponding to input proteins. Listed with a ';'
#' @export
#'
#' @import dplyr
mapToSpecies <- function(x, type = ""){
  map <- invisible(map)

  hash <- hash(keys = map$Entry, values = map$Organism)

  pboptions(type = "timer", char = "=", style = 3)
  mapped_species <- pbsapply(str_split(x, ";"), function(protein_ids) ## pbsapply instead of sapply so that there is a progress bar
  { ### creates list of protein names

    species <- sapply(str_match(protein_ids, "^[^-]*"), function(protein_id)
    { # individual protein level

      mapped_species <- hash[[str_match(protein_id,"^[^_]*")]] # searches hash table

      mapped_species <- mapped_species[mapped_species != "" & mapped_species != " " &
                                         mapped_species != "character(0)" & mapped_species != "NULL"]

      return(as.character(mapped_species)) # returns gene

    })

    if(type == "unique")
    { # for use in the Gene.names column, not for appending PTMs to name
      species <- unique(species)
    }

    species <- species[species != "" & species != " " & species != "character(0)" & species != "NULL"] # cleaning out non-matches

    return(str_c(species, collapse = ";", sep = NULL)) # returns list of genes for one row of proteins
  })

  mapped_species <- mapped_species %>% replace_na("")
  return(trimws(mapped_species))
}
