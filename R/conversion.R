#' Sorts the Proteins and Gene names columns in alphanumeric order
#'
#' @param df tibble or data frame with Proteins and Gene names columns, e.g. evidence or sites tables from MaxQuant
#'
#' @export
#' @examples
#' \dontrun{
#' evidence <- read_tsv("evidence.txt")
#' evidence <- sortIdentifiers(evidence)
#'
#' sites <- read_tsv("Phospho (STY)Sites.txt")
#' sites <- sortIdentifiers(sites)
#' }
#'
#' @import dplyr
#'
sortIdentifiers <- function(df) {
  df <- df %>%
    mutate("Sorted Proteins" = sortMultiColumn(Proteins)) %>%
    mutate("Sorted Gene names" = sortMultiColumn(Gene.names))

  return(df)
}

sortMultiColumn <- function(x) {
  # pboptions(type = "timer", char = "=", style = 3)
  sapply(str_split(x, ";"), function(p) str_c(sort(p), collapse=";"))
}

sortByPattern <- function(x, pattern) {
  # matched by gene pattern: "(^[^_]*)_([^_]*)_([^_]*)"
  # match by protein pattern: "(^[^_]*)_([^_]*)"

  # pboptions(type = "timer", char = "=", style = 3)
  y <- sapply(str_split(x, ';'), function(x){ # take in one row of annotated names
    # print(x)
    df <- as.data.frame(str_match(x, pattern)) # break a list of PTMs into gene, protein, and position
    df <- df[order(df[,2]),] # sort by gene name only
    return(str_c(df[,1], collapse = ';', sep = NULL))
  })
  return(y)
}

mapToGene <- function(x, type = ""){
  map <- invisible(mapped)
  setkey(map, Entry)

  pboptions(type = "timer", char = "=", style = 3)
  mapped_genes <- pbsapply(str_split(x, ";"), function(protein_ids) ## pbsapply instead of sapply so that there is a progress bar
  { ### creates list of protein names

    genes <- sapply(protein_ids, function(protein_id)
    { # individual protein level
      # mapped <- map %>% filter(Entry %in% str_match(protein_id, "^[^-]*")) %>% select(`gene.name`) # find gene name, step to be chagned to a hash sort
      protein_id <- str_match(protein_id, "^[^-]*") # retrieves canonical protein name
      mapped_protein <- map[.(protein_id)] # maybe a hash table that uses the Entry (protein name) as the key
      mapped_protein <- mapped_protein %>% filter(gene.name != "" & gene.name != " " & # cleaning the non-matches
                                    gene.name != "NULL" & gene.name != "character(0)" &
                                    !is.null(gene.name) &
                                    !is.na(gene.name))
      return(str_c(as.character(mapped_protein$gene.name), collapse = ";", sep = NULL)) # returns list of genes for one protein

    })

    if(type == "unique")
    { # for use in the Gene.names column, not for apending PTMs to name
      genes <- unique(genes)
    }

    genes <- genes[genes != "" & genes != " " & genes != "character(0)" & genes != "NULL"] # cleaning out non-matches

    if(length(genes) == 0)
    {
      return("")
    }
    # print(genes)
    return(str_c(genes, collapse = ";", sep = NULL)) # returns list of genes for one row of proteins
  })

  return(mapped_genes)
}

update.gene.names <- function(df) {
  df <- df %>%
    mutate("Gene.names" = mapToGene(Proteins, "unique"))

  return(df)
}

addSitesToNames <- function(ids, positions, amino_acid) {
  split_ids = str_split(ids, ";")
  split_positions = str_split(positions, ";")

  ids_with_sites <- mapply(function(s_ids, s_pos, aa) {
    str_c(s_ids, "_", aa, s_pos, collapse=";")
  }, s_ids=split_ids, s_pos=split_positions, aa=amino_acid)

  return(ids_with_sites)
}

addSitesToGeneNames <- function(ids, positions, amino_acid) {
  split_ids <- str_split(ids, ";")
  split_genes <- str_split(mapToGene(ids), ";") # returns redundant list of genes but it's unsorted

  # message("\n== Getting gene names for gene based identifiers")
  genes_and_ids <- mapply(function(protein, gene) {
    str_c(gene, "_", protein, collapse = ';')
  }, protein = split_ids, gene = split_genes)

  # returning only id column
  return(addSitesToNames(genes_and_ids, positions, amino_acid))
}

prepareForMSstats <- function(phosphosites, global_evidence, min_match=2, min_global_nonMissing=2) {
  # updating gene names
  message("\n== Updating gene names for phosphosites file")
  phosphosites <- update.gene.names(phosphosites)
  message("\n== Updating gene names for global evidence file")
  global_evidence <- update.gene.names(global_evidence)

  # order protein IDs and gene names
  message("\n== Sorting gene and protein names for the global evidence file")
  global_evidence <- sortIdentifiers(global_evidence)
  message("\n== Sorting gene and protein names for the phosphosites file")
  phosphosites <- sortIdentifiers(phosphosites)

  # global_evidence <- global_evidence %>% mutate(`Sorted Gene names` = Gene.names, `Sorted Proteins` = Proteins)
  # phosphosites <- phosphosites %>% mutate(`Sorted Gene names` = Gene.names, `Sorted Proteins` = Proteins)

  message("\n** Removing rows with too many missing values from evidence file")
  # filter out rows with too many missing values
  global_evidence <- global_evidence %>% filter(rowSums(across(matches("Reporter.intensity.corrected.\\d+"), function(x) x > 0)) >= min_global_nonMissing)
  # find phosphosites with enough distinct peptides in the global data that have exact protein matches

  message("\n** Classifying identifer type")
   matching_phospho_protein <- phosphosites %>%
    inner_join(global_evidence, by = c("Sorted Proteins")) %>%
    filter(`Sorted Proteins` != "" & "Sorted Proteins" != " " & !is.na("Sorted Proteins")) %>%
    group_by(id.x, `Sorted Proteins`) %>%
    summarize(num_matches = n()) %>%
    filter(num_matches >= min_match)

  # find phosphosites that didn't pass the previous filter, but do have gene-level matches
  matching_phospho_gene <- phosphosites %>%
    filter(!id %in% matching_phospho_protein$id.x) %>%
    inner_join(global_evidence, by = c("Sorted Gene names")) %>%
    filter(`Sorted Gene names` != "" & "Sorted Gene names" != " " & !is.na("Sorted Gene names")) %>%
    group_by(id.x, `Sorted Gene names`) %>%
    summarize(num_matches = n())

  not_matching_phospho <- phosphosites %>%
    filter(!id %in% matching_phospho_protein$id.x) %>%
    filter(!id %in% matching_phospho_gene$id.x) %>%
    select(id, Proteins, Gene.names) %>%
    rename("id.x" = id,
           "Sorted Proteins" = Proteins,
           "Sorted Gene names" = "Gene.names")

  message("\n** Creating new mapping identifiers")
  # create new identifiers for global for MSstats
  # evidence that matches proteins exactly
  global_evidence_protein <- global_evidence %>%
    filter(`Sorted Proteins` %in% matching_phospho_protein$`Sorted Proteins`) %>%
    mutate("matched_ids" = `Sorted Proteins`) %>%
    mutate(Modified.sequence = str_c(Modified.sequence, "protein"))
  # evidence that matches genes exactly (will double count some rows with the previous)
  global_evidence_gene <- global_evidence %>%
    filter( `Sorted Gene names` %in% matching_phospho_gene$`Sorted Gene names`) %>%
    mutate("matched_ids" =  `Sorted Gene names`) %>%
    mutate(Modified.sequence = str_c(Modified.sequence, "gene"))
  # evidence that didn't match anything (still keep it for normalization procedures)
  global_evidence_no_match <- global_evidence %>%
    filter(!`Sorted Proteins` %in% matching_phospho_protein$`Sorted Proteins` &
           !`Sorted Gene names` %in% matching_phospho_gene$`Sorted Gene names`) %>%
    mutate("matched_ids" =  `Sorted Proteins`) %>%
    mutate(Modified.sequence = str_c(Modified.sequence, "protein"))
  # merge all (doesn't include unquantified evidence)
  global_evidence <- bind_rows(global_evidence_gene, global_evidence_protein, global_evidence_no_match)

  # create new identifier columns for phospho for MSstats
  message("\n** Updating identifiers in phosphosites table")
  phosphosites_protein <- phosphosites %>%
    filter(id %in% matching_phospho_protein$id.x) %>%
    mutate("matched_ids" = addSitesToNames(Proteins, Positions.within.proteins, Amino.acid)) %>%
    mutate("matched_ids" = sortByPattern(matched_ids, pattern = "(^[^_]*)_([^_]*)"))

  phosphosites_genes <- phosphosites %>%
    filter(id %in% matching_phospho_gene$id.x) %>%
    mutate("matched_ids" = addSitesToGeneNames(Proteins, Positions.within.proteins, Amino.acid)) %>%
    mutate("matched_ids" = sortByPattern(matched_ids, pattern = "(^[^_]*)_([^_]*)_([^_]*)"))

  phosphosites_no_match <- phosphosites %>%
    filter(id %in% not_matching_phospho$id.x) %>%
    mutate("matched_ids" = addSitesToNames(Proteins, Positions.within.proteins, Amino.acid)) %>%
    mutate("matched_ids" =  sortByPattern(matched_ids, pattern = "(^[^_]*)_([^_]*)"))

  phosphosites <- bind_rows(phosphosites_genes, phosphosites_protein, phosphosites_no_match)

  global_evidence <- global_evidence %>% select(-Proteins) %>% rename(Proteins = matched_ids)
  phosphosites <- phosphosites %>% select(-Proteins) %>% rename(Proteins = matched_ids)
  return(list(global_evidence, phosphosites))

}
