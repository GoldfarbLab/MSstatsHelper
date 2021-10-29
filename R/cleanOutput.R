#' Changes protein names in the final output of MSstats to a readable format
#'
#' @param df the output of the comparison function and normalization level used -- pg or nothing
#'
#' @export
#' @examples
#' \dontrun{
#' sites <- read_tsv("Phospho (STY)Sites.txt", sep = '\t')
#' evidence <- read.delim("evidence.txt", sep = '\t')
#'
#' new_sites <- sitesToEvidence(sites, evidence)
#' }
#'
#' @import dplyr

cleanOutput <- function(comparison, protein_groups, names_col = 'Protein', level = 'p', normalization = ''){
  ### cleaning output for novel normalization
  if(normalization == "" & level == 'p'){
    df <- comparison %>%
      mutate(identifiers = sapply(str_split(!!sym(names_col), ';'), function(x){
        x <- sapply(x, function(y) str_split(y, '_')[[1]][1])
        return(str_c(x, collapse = ';'))
      })) %>%
      mutate(Site = sapply(str_split(!!sym(names_col), ';'), function(x){
        x <- sapply(x, function(y) str_split(y, '_')[[1]][length(str_split(y, '_')[[1]])])
        return(str_c(x, collapse = ';'))
      }))

    df <- df %>%
      mutate("Genes_in_protein" = mapToGene(!!sym(names_col), "nonunique"))

    df <- df %>%
      mutate(Genes_in_protein = ifelse(Genes_in_protein == "", identifiers, Genes_in_protein)) %>%
      mutate(Protein = ifelse(Genes_in_protein == identifiers, # meaning we normalized to gene not protein
                                   sapply(str_split(!!sym(names_col), ';'), function(x){ #then we need to retrieve protein name
                                     x <- sapply(x, function(y) str_split(y, '_')[[1]][2])
                                     return(str_c(x, collapse = ';'))
                                     # " ",
                                   }),
                                   identifiers)) %>% # otherwise we normalized to protein and the identifier is the protein name
      select(-identifiers)

    df <- df %>%
      filter(Site != Protein) %>%
      mutate("Genes_in_protein" = mapToGene(Protein, "Genes_in_protein"))

    return(df)
  }

  else if(tolower(normalization) %in% c("pg", 'protein group', 'proteingroup') & level == 'p'){
    df <- comparison %>%
      mutate(proteinGroup = sapply(str_split(!!sym(names_col), ';'), function(x){
        x <- sapply(x, function(y) str_split(y, '_')[[1]][1])
        x <- str_extract(x, pattern = '[0-9]*')
        return(as.numeric(str_c(unique(x), collapse = ';')))
      })) %>%

      mutate(Protein = sapply(str_split(!!sym(names_col), ';'), function(x){
        x <- sapply(x, function(y) str_split(y, '_')[[1]][2])
        return(str_c(x, collapse = ';'))
      })) %>%

      mutate(Site = sapply(str_split(!!sym(names_col), ';'), function(x){
        x <- sapply(x, function(y) str_split(y, '_')[[1]][3])
        return(str_c(x, collapse = ';'))
      })) %>%

      mutate(Genes_in_protein = mapToGene(Protein, 'unique'))

    protein_groups <- formatColNames(protein_groups)
    razor.peptide.ids <- protein_groups %>%
      filter(!Protein.IDs %in% c("", " ", "NULL")) %>%
      select(`Peptide.is.razor`, `Peptide.IDs`, `id`, `Protein.IDs`) %>%
      separate_rows(c(`Peptide.is.razor`, `Peptide.IDs`), sep=";", convert=T) %>%
      filter(`Peptide.is.razor` == "True") %>%
      select(`id`, `Protein.IDs`) %>%
      mutate(Genes_in_pg = mapToGene(Protein.IDs, "unique"),
             Proteins_in_pg = str_replace_all(Protein.IDs, ';', '-')) %>%
      select(-Protein.IDs) %>%
      unique()

    df <- df %>%
      left_join(razor.peptide.ids, by = c("proteinGroup" = "id"))

    return(df)
  }

  else if(level == 'wc'){
    df <- df %>%
      mutate("Gene.names" = mapToGene(Protein, "unique"))

    return(df)
  }

  else{
    message("\n** Unknown normalization method.
            \n*** For phospho level data: Set normalization to ProteinGroup or nothing and p for level.
            \n*** For whole cell level data: Set level to wc")
    return()
  }
}
