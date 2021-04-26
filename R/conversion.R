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
    mutate("Sorted Gene names" = sortMultiColumn(`Gene names`))

  return(df)
}

sortMultiColumn <- function(x) {
  sapply(str_split(x, ";"), function(p) str_c(sort(p), collapse=";"))
}

prepareForMSstats <- function(phosphosites, global_evidence, min_match=2, min_global_nonMissing=2) {
  # filter out rows with too many missing values
  global_evidence <- global_evidence %>% filter(rowSums(across(matches("Reporter intensity corrected \\d+"), function(x) x > 0)) >= min_global_nonMissing)

  # find phosphosites with enough distinct peptides in the global data that have exact protein matches
  matching_phospho_protein <- phosphosites %>%
    inner_join(global_evidence, by = c("Sorted Proteins")) %>%
    group_by(id.x, `Sorted Proteins`) %>%
    summarize(num_matches = n()) %>%
    filter(num_matches >= min_match)

  # remove global data that was used in the previous so we don't double count them
  #global_evidence <- global_evidence %>%
  #  filter(!id %in% matching_phospho_protein$id.y)

  # find phosphosites that didn't pass the previous filter, but do have gene-level matches
  matching_phospho_gene <- phosphosites %>%
    filter(!id %in% matching_phospho_protein$id.x) %>%
    filter("Sorted Gene names" != "") %>%
    inner_join(global_evidence, by = c("Sorted Gene names")) %>%
    group_by(id.x, `Sorted Gene names`) %>%
    summarize(num_matches = n())

  # create new identifiers for global for MSstats
  # evidence that matches proteins exactly
  global_evidence_protein <- global_evidence %>%
    filter(`Sorted Proteins` %in% matching_phospho_protein$`Sorted Proteins`) %>%
    mutate("matched_ids" = `Sorted Proteins`)
  # evidence that matches genes exactly (will double count some rows with the previous)
  global_evidence_gene <- global_evidence %>%
    filter( `Sorted Gene names` %in% matching_phospho_gene$`Sorted Gene names`) %>%
    mutate("matched_ids" =  `Sorted Gene names`)
  # evidence that didn't match anything (still keep it for normalization procedures)
  global_evidence_no_match <- global_evidence %>%
    filter(!`Sorted Proteins` %in% matching_phospho_protein$`Sorted Proteins` &
           !`Sorted Gene names` %in% matching_phospho_gene$`Sorted Gene names`) %>%
    mutate("matched_ids" =  `Sorted Proteins`)
  # merge all (doesn't include unquantified evidence)
  global_evidence <- bind_rows(global_evidence_gene, global_evidence_protein, global_evidence_gene)

  # create new identifier columns for phospho for MSstats
  phosphosites_protein <- phosphosites %>%
    filter(id %in% matching_phospho_protein$id.x) %>%
    mutate("matched_ids" = addSitesToNames(Proteins, `Positions within proteins`, `Amino acid`)) %>%
    mutate("matched_ids" = sortMultiColumn(matched_ids))

  phosphosites_genes <- phosphosites %>%
    filter(id %in% matching_phospho_gene$id.x) %>%
    mutate("matched_ids" = addSitesToGeneNames(`Gene names`, Proteins, `Positions within proteins`, `Amino acid`)) %>%
    mutate("matched_ids" = sortMultiColumn(matched_ids))

  phosphosites <- bind_rows(phosphosites_genes, phosphosites_protein)

  return(list(global_evidence, phosphosites))

}

addSitesToNames <- function(ids, positions, amino_acid) {
  split_ids = str_split(ids, ";")
  split_positions = str_split(positions, ";")

  ids_with_sites <- mapply(function(s_ids, s_pos, aa) {
    str_c(s_ids, "_", aa, s_pos, collapse=";")
  }, s_ids=split_ids, s_pos=split_positions, aa=amino_acid)

  return(ids_with_sites)
}

addSitesToGeneNames <- function(genes, ids, positions, amino_acid) {
  split_genes = str_split(genes, ";")
  split_ids = str_split(ids, ";")
  split_positions = str_split(positions, ";")

  ids_with_sites <- mapply(function(s_genes, s_ids, s_pos, aa) {
    unique_ids <- unique(sapply(s_ids, function(x) str_replace(x, "-.*", "")))

    unique_indices = sapply(unique_ids, function(x) min(which(str_detect(s_ids, x))))

    #if (length(s_genes) != length(unique_indices)) {
    #  print(s_genes)
    #  print(s_ids)
    #  print(unique_ids)
    #  print(unique_indices)
    #}

    str_c(s_genes, "_", s_ids[unique_indices], "_", aa, s_pos[unique_indices], collapse=";")
  }, s_genes=split_genes, s_ids=split_ids, s_pos=split_positions, aa=amino_acid)

  return(ids_with_sites)
}




