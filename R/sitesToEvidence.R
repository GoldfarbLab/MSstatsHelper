#' Changes MaxQuant's phosphosites (STY) file so that it can be used in place of evidence.txt in the MSstats Pipeline
#'
#' @param df tibble or data frame with Proteins and Gene names columns, e.g. evidence or sites tables from MaxQuant
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

sitesToEvidence <- function(sites, evidence, proteinGroups) {

  sites <- formatColNames(sites)
  contains_sn <- max(str_detect(colnames(sites), "SN"))
  evidence <- formatColNames(evidence)
  proteinGroups <- formatColNames(proteinGroups)

  ##############################################################################
  ### sum together phosphosites into one intensity column
  ##############################################################################
  message("\n** Merging data from site level to peptide level")
  map <- sites %>%
    select(id, Evidence.IDs, Proteins, Protein.group.IDs)

  sn.cols <- colnames(sites)[grep("SN\\.", colnames(sites))]
  is.missing.cols <- colnames(sites)[grep("Is\\.missing\\.", colnames(sites))]
  r.i.c.cols <- colnames(sites)[grep("Reporter\\.intensity\\.corrected\\.", colnames(sites))]
  r.i.cols <- colnames(sites)[grep("Reporter\\.intensity\\.", colnames(sites))]
  loc.cols <- colnames(sites)[grep("Localization\\.prob\\.", colnames(sites))]
  keep <- c(sn.cols, is.missing.cols, r.i.c.cols, r.i.cols, loc.cols)

  # reformat data
  cols.desc <- unique(tibble(cols = colnames(sites))) %>%
    mutate(Experiment = str_match(cols, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,4],
           Experiment = ifelse(is.na(Experiment), str_match(cols, pattern = "(Localization\\.prob\\.)(.*)")[,3], Experiment),
           Channel = str_match(cols, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,3],
           Column.type = str_c(str_match(cols, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,2]),
           Column.type = ifelse(is.na(Column.type), 'Localization.prob', Column.type),
           Phosphosite = str_match(cols, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,5],
           Phosphosite = ifelse(is.na(Phosphosite), 1, Phosphosite)) %>%
    # filter(!is.na(Experiment)) %>%
    mutate(new.col.name = ifelse(is.na(Channel),
                                 Column.type,
                                 str_c(Column.type, ".", Channel))) %>%
    select(-Channel, -Column.type) %>%
    filter(!is.na(Experiment))

  sites <-  sites %>%
    select(id, all_of(keep)) %>%
    pivot_longer(cols = all_of(keep),
                 values_to = "Value",
                 names_to = 'old.name') %>%
    inner_join(cols.desc, by = c('old.name' = 'cols')) %>%
    filter(!is.na(Experiment)) %>%
    select(-old.name) %>%
    pivot_wider(names_from = new.col.name, values_from = Value)

  filter.col <- colnames(sites)[grep("Is\\.missing\\.", colnames(sites))]
  localization.col <- colnames(sites)[grep("Localization\\.prob", colnames(sites))]
  reporter.corr.col <- colnames(sites)[grep("Reporter\\.intensity\\.corrected\\.", colnames(sites))]
  reporter.col <- colnames(sites)[grep("Reporter\\.intensity\\.[0-9]+", colnames(sites))]
  SN.col <- colnames(sites)[grep("SN\\.", colnames(sites))]

  if(contains_sn){
     summed <- sites %>%
      group_by(Experiment, id) %>%
      summarise_at(c(all_of(SN.col), all_of(reporter.col),
                   all_of(reporter.corr.col), all_of(filter.col)), funs(sum)) %>%
       mutate_at(all_of(filter.col), funs(ifelse(. < 3, FALSE, TRUE)))
                # Localization.prob = max(Localization.prob)) ### because each peptide has only one localization prob

     maxed <- sites %>%
       group_by(Experiment, id) %>%
       summarise_at(all_of(localization.col), funs(max), na.rm=TRUE) %>%
       filter(!Localization.prob %in% c(-Inf, Inf))
  }

  else {
    summed <- sites %>%
      group_by(Experiment, id) %>%
      summarise_at(c(all_of(reporter.col),
                     all_of(reporter.corr.col), all_of(filter.col)), funs(sum)) %>%
      mutate_at(all_of(filter.col), funs(ifelse(. < 3, 0, 1))) %>%
      by_row(~ sum(is.na(.)), .collate = 'cols')
    # Localization.prob = max(Localization.prob)) ### because each peptide has only one localization prob

    maxed <- sites %>%
      group_by(Experiment, id) %>%
      summarise_at(all_of(localization.col), funs(max), na.rm=TRUE) %>%
      filter(!Localization.prob %in% c(-Inf, Inf))
  }

  sites <- inner_join(summed, maxed, by=c('id', 'Experiment')) %>%
    inner_join(map, by = 'id')

  ##############################################################################
  ### merge sites and evidence
  ##############################################################################
  sites <- sites %>%
    rename('Phospho.STY.site.IDs' = 'id') %>%
    separate_rows(Evidence.IDs, sep = ';', convert=TRUE)

  remove <- c(filter.col, reporter.col, reporter.corr.col, SN.col)
  evidence <- evidence %>%
    select(-all_of(remove), -Proteins, -Protein.group.IDs) %>%
    # select('Modified.sequence', 'Raw.file', 'Phospho..STY..site.IDs', 'id') %>%
    rename('Evidence.IDs' = 'id') %>%
    separate_rows(Phospho.STY.site.IDs, sep = ';', convert = TRUE)


  input <- sites %>% left_join(evidence, by = c("Phospho.STY.site.IDs", "Evidence.IDs", "Experiment"))
  input <- input %>%
    group_by(Phospho.STY.site.IDs, Experiment) %>%
    arrange(Localization.prob, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup

  input <- input %>% ## redo IDs for table -- multiple experiments makes phosphoID redundant
    mutate(id = row_number()) %>%
    select(-Evidence.IDs) %>%
    ungroup()

  ##############################################################################
  ### update Evidence ID in the proteinGroups table
  ##############################################################################
  map <- input %>%
    # filter(!str_detect(Protein.group.IDs, ';')) %>%
    select(Protein.group.IDs, id) %>%
    separate_rows(Protein.group.IDs, sep=';', convert=TRUE) %>%
    unique() %>%
    group_by(Protein.group.IDs) %>%
    summarize(Evidence.IDs = str_c(id, collapse = ';'))

  proteinGroups <- proteinGroups %>% select(-Evidence.IDs)
  proteinGroups <- proteinGroups %>%
    left_join(map, by = c('id' = 'Protein.group.IDs'))

  return(list("Sites" = input, "proteinGroups" = proteinGroups))
}
