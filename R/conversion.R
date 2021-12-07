#' Sorts 'Protein" and "Gene.names" column
#' @param df tibble or data frame containing columns "Proteins" and "Gene.names"
#' @return tibble or data frame with sorted "Proteins" and "Gene.names" columns
#' @export
#' @examples
#'
#' @import dplyr
sortIdentifiers <- function(df) {
  df <- df %>%
    mutate("Sorted Proteins" = sortMultiColumn(Proteins)) %>%
    mutate("Sorted Gene names" = sortMultiColumn(Gene.names))

  return(df)
}

#' Updates the Gene names and sorts the Proteins and Gene names columns in alphanumeric order
#'
#' @param x character column with words separated by ';'
#' @return column sorted by words
#'
#' @export
#' @examples
#'
#' @import dplyr
sortMultiColumn <- function(x) {
  sapply(str_split(x, ";"), function(p) str_c(sort(p), collapse=";"))
}

#' Sorts a column based on some component of a pattern
#'
#' @param x character column with words separated by ';'
#' @param pattern regular expression with capture groups representing each component of the word to be sorted
#' @param n_sort specifies which capture group to sort the words by
#' @return column sorted by words
#'
#' @export
#' @examples
#'
#' @import dplyr
sortByPattern <- function(x, pattern, n_sort = 1) {
  # matched by gene pattern: "(^[^_]*)_([^_]*)_([^_]*)"
  # match by protein pattern: "(^[^_]*)_([^_]*)"

  y <- sapply(str_split(x, ';'), function(x){ # take in one row of annotated names
    # message(x)
    df <- as.data.frame(str_match(x, pattern)) # break a list of PTMs into gene, protein, and position
    df <- df[order(df[,(n_sort+1)]),] # sort by gene name only
    return(str_c(df[,1], collapse = ';', sep = NULL))
  })
  return(y)
}

#' Adds column of genes mapped to protein columns
#'
#' @param df data frame or tibble contianing the column "Proteins"
#' @return data frame or tibble with included genes columns
#'
#' @export
#' @examples
#'
#' @import dplyr
update.gene.names <- function(df) {
  df <- df %>%
    filter(!Proteins %in% c("", " ", "NULL")) %>%
    mutate(Proteins = str_replace_all(Proteins, " ", "")) %>%
    mutate("Gene.names" = mapToGene(Proteins, "unique"))

  return(df)
}

#' appends protein group to protein names
#'
#' @param Protein.group.IDs column of protein group identifiers
#' @param Proteins column containing protein names
#' @return column with Protein.group.IDs_Proteins
#'
#' @export
#' @examples
#'
#' @import dplyr
addPGtoNames <- function(Protein.group.IDs, Proteins) {
  split_proteins = str_split(Proteins, ";")

  ids_with_sites <- mapply(function(pro, id) {
    str_c(id, "_",  pro, collapse=";")
  }, pro=split_proteins, id=Protein.group.IDs)

  return(ids_with_sites)
}

#' appends site and amino acid to protein name
#'
#' @param ids column containing protein name or identifier
#' @param positions column containing the position of phosphorylated amino acid in protein
#' @param amino_acid column containing letter of phosphorylated amino acid
#' @return column with id_aa_position;id_aa_position...
#'
#' @export
#' @examples
#'
#' @import dplyr
addSitesToNames <- function(ids, positions, amino_acid) {
  split_ids = str_split(ids, ";")
  split_positions = str_split(positions, ";")

  ids_with_sites <- mapply(function(s_ids, s_pos, aa) {
    str_c(s_ids, "_", aa, s_pos, collapse=";")
  }, s_ids=split_ids, s_pos=split_positions, aa=amino_acid)

  return(ids_with_sites)
}

#' appends site to geneName_proteinName
#'
#' @param ids column containing protein name or identifier
#' @param positions column containing the position of phosphorylated amino acid in protein
#' @param amino_acid column containing letter of phosphorylated amino acid
#' @return column with geneName_id_aa_position;geneName_id_aa_position...
#'
#' @export
#' @examples
#'
#' @import dplyr
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

#' if annotation file contains a reference channel, the input data is filtered so that each row in evidence file contains a non-missing value for at least one reference channel
#'
#' @param df tibble or data frame from MQ evidence or proteinGroups file
#' @param type string. If file is an evidence or proteinGroup file specify 'ev' if it's a phosphosite file specify 'ph'
#' @param annotation tibble or data frame of annotation file corresponding with input file
#' @param profiel boolean. if data was made in profile mode or not
#'
#' @export
#' @examples
#'
#' @import dplyr
checkReferenceImputation <- function(df, type = 'ev', annotation, profile){
  if(type == 'ev'){
    message("\n** Removing peptides missing in each reference channel from evidence file.")
    if(profile){

    }

    else{
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

  }

  else if(type == 'ph'){ ### only do this AFTER conveting sites file to evidence
    message("\n** Removing phosphosites missing in each reference channel from sites file.")

    if(profile){

    }

    else{
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
  }

  else{
    message("\nERROR: Incorrect file type in checkReferenceImputation")
    return()
  }

}

#' Removes rows with less nonmissing values than the specified number
#'
#' @param df tibble or data frame containing "Is missing" columns
#' @param min_measure default 3, the fewest measurements tolerated
#' @param profile boolean, if data was made in profile mode or not
#' @return data frame or tibble filtered or rows with min_measure or more
#'
#' @export
#' @examples
#'
#' @import dplyr
removeFewMeas <- function(df, min_measure = 3, profile){
  if(profile){
    df <- df %>%
      filter(rowSums(across(matches("Reporter.intensity.corrected.\\d+"), function(x) x > 0)) >= min_measure)
  }

  else{
    is.missing.cols <- colnames(df)[grep("Is\\.missing\\.", colnames(df))]

    df <- df %>%
      rowwise %>%
      mutate(measurements = sum(c_across(all_of(is.missing.cols)) == FALSE)) %>%
      filter(measurements >= min_measure) %>%
      select(-measurements)
  }
  return(df)

}

#'
#'
#' @param phosphosites tibble or data frame of phosphosites file output by TMTPurityCorrection package
#' @param phospho_annotation tibble or data frame of phospho data's annotation file
#' @param global_annotation tibble or data frame of global data's annotation file
#' @param global_evidence tibble or data frame of evidence file output by TMTPuritycorrection package
#' @param phospho_evidence tibble or data frame of evidence only corresponding to phosphosites
#' @param protein_groups tibble or data frame of protein groups file output by TMTPurityCorrection package
#' @param ptm_match string, if normalizing phospho by protein groups, specify "pg", otherwise do not specify
#' @param evidence_only TRUE/FALSE, if not using phospho file, specify TRUE
#' @param min_match numeric. for novel normalization, fewest number of matches to normalize by isoform or gene
#' @param min_measure numeric. fewest tolerated number of nonmissing values in row to be
#' @param profile_mode boolean if data was made in profile mode
#' @return If evidence only, returns file filtered as specified with the protein name as pgID_proteinname--proteinname2--...
#' If including phosphosite file and normalizing protein group returns phosphosite file that has been converted to evidence format with protein names as pgID_prot_aa_site;...
#' also returns the evidence file with names of format pgID and a protein group file with updated id's. Both files should be used for phosphosite normalization.
#' If including phosphosite file and using novel normalization, returns all as stated above but identifiers are different.
#' @export
#' @examples
#'
#' @import dplyr
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
                               profile_mode = FALSE) {
  # phosphosites <- formatColNames(phosphosites)
  # global_evidence <- formatColNames(global_evidence)

  if(!tolower(ptm_match) %in% c("", " ", "", 'proteingroup', 'pg', 'protein group')){
    message("ERROR: unknown ptm_match. Specify pg for normalization using MaxQuant\'s protein grouping or nothing to use protein and gene level normalization")
    return()
  }

  if(profile_mode){
    message("Data was run in profile mode. Filtering by missingness will not be performed.")
  }

  # filter out rows with too many missing values
  # message("\n** Removing rows with too many missing values from evidence file")
  global_evidence <- formatColNames(global_evidence)
  protein_groups <- formatColNames(protein_groups)

  if(!profile_mode){
    message("\n** Removing peptides with few measurements")
    global_evidence <- removeFewMeas(global_evidence, min_measure=min_measure, profile_mode)
  }

  if(!profile_mode & max(str_detect(tolower(global_annotation$Condition), "norm"))){
    message("** Reference channel found.")
    global_evidence <- checkReferenceImputation(global_evidence, annotation = global_annotation, type = 'ev', profile_mode)
  }

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
    phospho.as.evidence <- sitesToEvidence(phosphosites, evidence = phospho_evidence, proteinGroups = protein_groups, profile_mode)
    phosphosites <- phospho.as.evidence$Sites

    if(!profile_mode & max(str_detect(tolower(phospho_annotation$Condition), "norm"))){
      message("** Reference channel found.")
      phosphosites <- checkReferenceImputation(phosphosites, annotation = phospho_annotation, type = 'ph', profile_mode)
    }

    phosphosites <- phosphosites %>% filter(Localization.prob > 0.75)


    if(!profile_mode){
      message("\n** Removing phosphosites with few measurements")
      phosphosites <- removeFewMeas(phosphosites, min_measure=min_measure, profile_mode)
    }


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
      mutate(Proteins = str_c('pg', Protein.group.IDs)) %>%
      mutate(Proteins = Protein.group.IDs)

    phosphosites <- phosphosites %>%
      filter(Proteins != "") %>%
      separate_rows(`Protein.group.IDs`, sep=";", convert=T) %>%
      separate_rows(`Peptide.IDs`, sep=";", convert=T) %>%
      # mutate(`Peptide.IDs` = as.numeric(`Peptide.IDs`)) %>%
      inner_join(razor.peptide.ids, by=c("Peptide.IDs" = "Peptide.IDs", "Protein.group.IDs" = "id")) %>%
      mutate("Proteins" = addSitesToNames(Proteins, Positions.within.proteins, Amino.acid),
             "ProteinGroupIDs" = str_c('pg', Protein.group.IDs),
             # "Proteins" = addPGtoNames(Protein.IDs, Proteins)) %>%
             "Proteins" = addPGtoNames(Protein.group.IDs, Proteins)) %>%
      group_by(id) %>%
      slice(1)

    message("\n** Converting Sites to Evidence format")
    phospho.as.evidence <- sitesToEvidence(phosphosites, evidence = phospho_evidence, proteinGroups = protein_groups, profile_mode)

    phosphosites <- phospho.as.evidence$Sites
    if(!profile_mode & max(str_detect(tolower(phospho_annotation$Condition), "norm"))){
      message("** Reference channel found in phosphosites data.")
      phosphosites <- checkReferenceImputation(phosphosites, annotation = phospho_annotation, type = 'ph', profile_mode)
    }

    if(!profile_mode){
      message("\n** Removing phosphosites with few measurements")
      phosphosites <- removeFewMeas(phosphosites, min_measure=min_measure, profile_mode)
    }
    phosphosites <- phosphosites %>% filter(Localization.prob > 0.75)

    return(list("GlobalEvidence" = global_evidence, "Sites" = phosphosites, "ProteinGroups" = phospho.as.evidence$proteinGroups))
  }
}
