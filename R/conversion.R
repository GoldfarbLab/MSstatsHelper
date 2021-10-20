#' Updates the Gene names and sorts the Proteins and Gene names columns in alphanumeric order
#'
#' @param df tibble or data frame with Proteins and Gene names columns, e.g. evidence or sites tables from MaxQuant
#'
#' @export
#' @examples
#' \dontrun{
#' evidence <- read.delim("evidence.txt", sep = '\t')
#' sites <- read.delim("Phospho (STY)Sites.txt", sep = '\t')
#'
#' formatted_data <- prepareForMSstats(sites, evidence)
#' new_evidence <- formatted_data[[1]]
#' new_sites <- formatted_data[[2]]
#' }
#'
#' @import dplyr
#'
#'

formatColNames <- function(df) {
  cols <- colnames(df)
  cols <- str_replace_all(cols, '.', ' ')
  colnames(df) <- cols

  return(df)
}
sortIdentifiers <- function(df) {
  df <- df %>%
    mutate("Sorted Proteins" = sortMultiColumn(Proteins)) %>%
    mutate("Sorted Gene names" = sortMultiColumn(Gene.names))

  return(df)
}

sortMultiColumn <- function(x) {
  sapply(str_split(x, ";"), function(p) str_c(sort(p), collapse=";"))
}

sortByPattern <- function(x, pattern) {
  # matched by gene pattern: "(^[^_]*)_([^_]*)_([^_]*)"
  # match by protein pattern: "(^[^_]*)_([^_]*)"

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
  hash <- hash(keys = map$Entry, values = map$gene.name)

  pboptions(type = "timer", char = "=", style = 3)
  mapped_genes <- pbsapply(str_split(x, ";"), function(protein_ids) ## pbsapply instead of sapply so that there is a progress bar
  { ### creates list of protein names

    genes <- sapply(str_match(protein_ids, "^[^-]*"), function(protein_id)
    { # individual protein level
      if(protein_id == " " | protein_id == ""){return("")}
      mapped_protein <- hash[[str_match(protein_id,"^[^_]*")]] # searches hash table

      mapped_protein <- mapped_protein[mapped_protein != "" & mapped_protein != " " &
                                         mapped_protein != "character(0)" & mapped_protein != "NULL"]

      return(as.character(mapped_protein)) # returns gene
    })

    # print(genes)

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

update.gene.names <- function(df) {
  df <- df %>%
    mutate("Gene.names" = mapToGene(Proteins, "unique"))

  return(df)
}


addPGtoNames <- function(Protein.group.IDs, Proteins) {
  split_proteins = str_split(Proteins, ";")

  ids_with_sites <- mapply(function(pro, id) {
    str_c('pg', id, "_", pro, collapse=";")
  }, pro=split_proteins, id=Protein.group.IDs)

  return(ids_with_sites)
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

prepareForMSstats <- function(phosphosites = 0, global_evidence, protein_groups, ptm_match = '', min_match=2, min_global_nonMissing=2) {
  # phosphosites <- formatColNames(phosphosites)
  # global_evidence <- formatColNames(global_evidence)

  if(!tolower(ptm_match) %in% c("", " ", "", 'proteingroup', 'pg', 'protein group')){
    message("ERROR: unknown ptm_match. Specify pg for normalization using MaxQuant\'s protein grouping or nothing to use protein and gene level normalization")
    return()
  }

  # filter out rows with too many missing values
  message("\n** Removing rows with too many missing values from evidence file")
  global_evidence <- global_evidence %>%
    filter(rowSums(across(matches("Reporter.intensity.corrected.\\d+"), function(x) x > 0)) >= min_global_nonMissing)


  ### If the output will use the 'novel' level of phospho normalization
  if(ptm_match == ''){
    # order protein IDs and gene names
    phosphosites <- phosphosites %>% filter(Localization.prob > 0.75)
    message("\n** Sorting gene and protein names for the global evidence file")
    global_evidence <- sortIdentifiers(global_evidence)
    message("\n** Sorting gene and protein names for the phosphosites file")
    phosphosites <- sortIdentifiers(phosphosites)

    # updating gene names
    message("\n== Updating gene names for phosphosites file")
    phosphosites <- update.gene.names(phosphosites)
    message("\n== Updating gene names for global evidence file")
    global_evidence <- update.gene.names(global_evidence)

    message("\n** Classifying identifer type for Protein/Gene normalization")
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

    global_evidence <- global_evidence %>%
      select(-Proteins) %>%
      rename(Proteins = matched_ids) %>%
      filter(!str_detect(Protein.group.IDs, ';'))
    phosphosites <- phosphosites %>%
      select(-Proteins) %>%
      rename(Proteins = matched_ids) %>%
      filter(!str_detect(Protein.group.IDs, ';'))

    return(list("GlobalEvidence" = global_evidence, "Sites" = phosphosites))
  }

  ### If the output will use the protein group level of phospho normalization
  else if(tolower(ptm_match) %in% c("pg", 'protein group', 'proteingroup')){
    message("\n** Adding Protein Groups to Protein names for normalization using MQ's protein grouping")

    pgmap <- protein_groups %>%
      select(Protein.IDs, id) %>%
      mutate(Protein.IDs = str_replace_all(Protein.IDs, ';', '-'))

    phosphosites <- phosphosites %>%
      filter(Proteins != "") %>%
      filter(!str_detect(Protein.group.IDs, ';')) %>%
      mutate(Protein.group.IDs = as.numeric(Protein.group.IDs)) %>%
      # inner_join(pgmap, by=c("Protein.group.IDs" = "id")) %>%
      mutate("Proteins" = addSitesToNames(Proteins, Positions.within.proteins, Amino.acid),
             "Proteins" = addPGtoNames(Protein.group.IDs, Proteins)) #%>%
      # mutate("Proteins" = sortByPattern(Proteins, pattern = "([^_]*)_([^_]*)_([^_]*)"))

    global_evidence <- global_evidence %>%
      filter(Proteins != "") %>%
      filter(!str_detect(Protein.group.IDs, ';')) %>%
      mutate(Protein.group.IDs = as.numeric(Protein.group.IDs)) %>%
      mutate(Proteins = str_c('pg', Protein.group.IDs)) #%>%
      # inner_join(pgmap, by=c("Protein.group.IDs" = "id")) %>%
      # select(-Proteins) %>%
      # rename("Proteins" = "Protein.group.IDs")

    return(list("GlobalEvidence" = global_evidence, "Sites" = phosphosites))
  }

}
