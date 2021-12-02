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
  cols <- str_replace_all(cols, ' ', '.')
  cols <- str_replace_all(cols, '\\.\\.', '\\.')
  cols <- str_replace_all(cols, '\\(', '')
  cols <- str_replace_all(cols, '\\)', '')
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
    # message(x)
    df <- as.data.frame(str_match(x, pattern)) # break a list of PTMs into gene, protein, and position
    df <- df[order(df[,2]),] # sort by gene name only
    return(str_c(df[,1], collapse = ';', sep = NULL))
  })
  return(y)
}

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

update.gene.names <- function(df) {
  df <- df %>%
    filter(!Proteins %in% c("", " ", "NULL")) %>%
    mutate(Proteins = str_replace_all(Proteins, " ", "")) %>%
    mutate("Gene.names" = mapToGene(Proteins, "unique"))

  return(df)
}


addPGtoNames <- function(Protein.group.IDs, Proteins) {
  split_proteins = str_split(Proteins, ";")

  ids_with_sites <- mapply(function(pro, id) {
    str_c(id, "_",  pro, collapse=";")
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

checkReferenceImputation <- function(df, type = 'ev', annotation, summary = ''){
  if(type == 'ev'){
    message("\n** Removing peptides missing in each reference channel from evidence file.")
    annotation <- annotation %>%
      filter(str_detect(tolower(Condition), 'norm')) %>%
      mutate(number = str_c("Is.missing.", str_match(Channel, "(.*)\\.([0-9]*)")[,3]))

    reference.channels <- unique(annotation$number)

    ### filter any peptide that is missing in every in every reference channel
    df <- df %>%
      mutate(sums = select(., all_of(reference.channels)) %>% rowSums()) %>%
      filter(sums < length(reference.channels)) %>%
      select(-sums)

    return(df)
  }

  else if(type == 'ph'){ ### only do this AFTER conveting sites file to evidence
    message("\n** Removing phosphosites missing in each reference channel from sites file.")

    annotation <- annotation %>%
      filter(str_detect(tolower(Condition), 'norm')) %>%
      mutate(number = str_c("Is.missing.", str_match(Channel, "(.*)\\.([0-9]*)")[,3]))

    reference.channels <- unique(annotation$number)

    ### filter any peptide that is missing in every in every reference channel
    df <- df %>%
      mutate(sums = select(., all_of(reference.channels)) %>% rowSums()) %>%
      filter(sums < length(reference.channels)*3) %>%
      select(-sums)

    return(df)
  }

  else{
    message("\nERROR: Incorrect file type in checkReferenceImputation")
    return()
  }

}

removeFewMeas <- function(df, min_measure = 3){
  is.missing.cols <- colnames(df)[grep("Is\\.missing\\.", colnames(df))]

  df <- df %>%
    rowwise %>%
    mutate(measurements = sum(c_across(all_of(is.missing.cols)) == FALSE)) %>%
    filter(measurements >= min_measure) %>%
    select(-measurements)

  return(df)

}

prepareForMSstats <- function (phosphosites = '',
                               phospho_annotation = "",
                               global_annotation,
                               global_evidence,
                               phospho_evidence = '',
                               protein_groups,
                               ptm_match = '',
                               evidence_only = FALSE,
                               min_match=2,
                               min_measure=3,
                               min_global_nonMissing=2) {
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

  global_evidence <- formatColNames(global_evidence)
  protein_groups <- formatColNames(protein_groups)

  message("\n** Removing peptides with few measurements")
  global_evidence <- removeFewMeas(global_evidence, min_measure=min_measure)
  # if(max(str_detect(tolower(global_annotation$Condition), "norm"))){
  #   message("** Reference channel found in whole cell data.")
  #   global_evidence <- checkReferenceImputation(global_evidence, annotation = global_annotation, type = 'ev')
  # }

  if(evidence_only){
    message("\n== Updating Gene names in Protein Groups file")
    razor.peptide.ids <- protein_groups %>%
      filter(!Protein.IDs %in% c("", " ", "NULL")) %>%
      select(`Peptide.is.razor`, `Peptide.IDs`, `id`, `Protein.IDs`) %>%
      separate_rows(c(`Peptide.is.razor`, `Peptide.IDs`), sep=";", convert=T) %>%
      filter(`Peptide.is.razor` == "True") %>%
      select(`Peptide.IDs`, `id`, `Protein.IDs`) %>%
      mutate(Genes = mapToGene(Protein.IDs, "unique"),
             Protein.IDs = str_replace_all(Protein.IDs, ';', '--'))

    message("\n** Changing protein names to protein group names.")
    global_evidence <- global_evidence %>%
      separate_rows(`Protein.group.IDs`, sep=";") %>%
      mutate(`Protein.group.IDs` = as.numeric(`Protein.group.IDs`)) %>%
      inner_join(razor.peptide.ids, by=c("Peptide.ID" = "Peptide.IDs", "Protein.group.IDs" = "id")) %>%
      mutate(Proteins = str_c('pg', Protein.group.IDs, "_", Protein.IDs))

    return(list("GlobalEvidence" = global_evidence))
  }

  else if (ptm_match == ''){

    phosphosites <- formatColNames(phosphosites)
    ### If the output will use the 'novel' level of phospho normalization

    # updating gene names
    message("\n== Updating gene names for global evidence file")
    global_evidence <- update.gene.names(global_evidence)
    message("\n== Updating gene names for phosphosites file")
    phosphosites <- update.gene.names(phosphosites)

    # order protein IDs and gene names
    message("\n** Sorting gene and protein names for the global evidence file")
    global_evidence <- sortIdentifiers(global_evidence)
    message("\n** Sorting gene and protein names for the phosphosites file")
    phosphosites <- sortIdentifiers(phosphosites)

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

    message("\n** Converting Sites to Evidence format")
    phospho.as.evidence <- sitesToEvidence(phosphosites, evidence = phospho_evidence, proteinGroups = protein_groups)
    phosphosites <- phospho.as.evidence$Sites

    if(max(str_detect(tolower(phospho_annotation$Condition), "norm"))){
      message("** Reference channel found.")
      phosphosites <- checkReferenceImputation(phosphosites, annotation = phospho_annotation, type = 'ph')
      global_evidence <- checkReferenceImputation(global_evidence, annotation = global_annotation, type = 'ev')
    }

    phosphostes <- phosphosites %>% filter(Localization.prob > 0.75)
    message("\n** Removing phosphosites with few measurements")
    phosphosites <- removeFewMeas(phosphosites, min_measure=min_measure)

    return(list("GlobalEvidence" = global_evidence, "Sites" = phosphosites, "ProteinGroups" = phospho.as.evidence$proteinGroups))
  }

  ### If the output will use the protein group level of phospho normalization
  else if(tolower(ptm_match) %in% c("pg", 'protein group', 'proteingroup')){
    phosphosites <- formatColNames(phosphosites)
    message("\n** Adding Protein Groups to Protein names for normalization using MQ's protein grouping")

    message("\n== Updating Gene names in Protein Groups file.")
    razor.peptide.ids <- protein_groups %>%
      filter(!Protein.IDs %in% c("", " ", "NULL")) %>%
      select(`Peptide.is.razor`, `Peptide.IDs`, `id`, `Protein.IDs`) %>%
      separate_rows(c(`Peptide.is.razor`, `Peptide.IDs`), sep=";", convert=T) %>%
      filter(`Peptide.is.razor` == "True") %>%
      select(`Peptide.IDs`, `id`, `Protein.IDs`) %>%
      mutate(Genes = mapToGene(Protein.IDs, "unique"),
             Protein.IDs = str_replace_all(Protein.IDs, ';', '--'))

    message("\n** Changing protein names to protein group names.")
    global_evidence <- global_evidence %>%
      separate_rows(`Protein.group.IDs`, sep=";", convert=T) %>%
      separate_rows(`Peptide.ID`, sep=";", convert=T) %>%
      inner_join(razor.peptide.ids, by=c("Peptide.ID" = "Peptide.IDs", "Protein.group.IDs" = "id")) %>%
      # mutate(Proteins = Protein.IDs)
      # mutate(Proteins = str_c('pg', Protein.group.IDs, "-", Protein.IDs))
      mutate(Proteins = Protein.group.IDs)

    phosphosites <- phosphosites %>%
      filter(Proteins != "") %>%
      separate_rows(`Protein.group.IDs`, sep=";", convert=T) %>%
      separate_rows(`Peptide.IDs`, sep=";", convert=T) %>% ### creates multiple rows for one id.... should sty's have one to many peptide Id mapping?
      # mutate(`Peptide.IDs` = as.numeric(`Peptide.IDs`)) %>%
      inner_join(razor.peptide.ids, by=c("Peptide.IDs" = "Peptide.IDs", "Protein.group.IDs" = "id")) %>%
      mutate("Proteins" = addSitesToNames(Proteins, Positions.within.proteins, Amino.acid),
             # "ProteinGroupIDs" = str_c('pg', Protein.group.IDs, "-", Protein.IDs),
             # "Proteins" = addPGtoNames(Protein.IDs, Proteins)) %>%
             "Proteins" = addPGtoNames(Protein.group.IDs, Proteins)) %>%
      group_by(id) %>%
      slice(1)

    message("\n** Converting Sites to Evidence format")
    phospho.as.evidence <- sitesToEvidence(phosphosites, evidence = phospho_evidence, proteinGroups = protein_groups)
    phosphosites <- phospho.as.evidence$Sites

    if(max(str_detect(tolower(phospho_annotation$Condition), "norm"))){
      message("** Reference channel found in phosphosites data.")
      phosphosites <- checkReferenceImputation(phosphosites, annotation = phospho_annotation, type = 'ph')
    }

    phosphosites <- phosphosites %>% filter(Localization.prob > 0.75)

    message("\n** Removing phosphosites with few measurements")
    phosphosites <- removeFewMeas(phosphosites, min_measure=min_measure)

    return(list("GlobalEvidence" = global_evidence, "Sites" = phosphosites, "ProteinGroups" = phospho.as.evidence$proteinGroups))
  }
}
