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
  evidence <- formatColNames(evidence)
  proteinGroups <- formatColNames(proteinGroups)

  ##############################################################################
  ### sum together phosphosites into one intensity column
  ##############################################################################
  ## this all worked for singular TMT plexes but doesn't work for multiple/multiple experiments.
  # I need to transform the separate columns for experiment abundances into rows with the experiment labeled
  # for intensity, intensity corrected, SN, Is.missing

  r.i.cols <- colnames(sites)[grep("Reporter\\.intensity\\.([0-9]*)\\.", colnames(sites))]
  experiments <- unique(str_match(r.i.cols, pattern = "(Reporter\\.intensity\\.([0-9]*)\\.)(.*)(___([0-9]*))")[,4])
  channels <- unique(str_match(r.i.cols, pattern = "(Reporter\\.intensity\\.([0-9]*)\\.)(.*)(___([0-9]*))")[,3])

  new.intensities <- list()

  pb = txtProgressBar(max=length(channels), style=3)
  for(i in 1:length(channels)){

    # Grab column names for this particular channel
    sn.cols <- colnames(sites)[grep(str_c("SN\\.", channels[i], "\\."), colnames(sites))]
    is.missing.cols <- colnames(sites)[grep(str_c("Is\\.missing\\.", channels[i], '\\.'), colnames(sites))]
    r.i.c.cols <- colnames(sites)[grep(str_c("Reporter\\.intensity\\.corrected\\.", channels[i], "\\."), colnames(sites))]
    r.i.cols <- colnames(sites)[grep(str_c("Reporter\\.intensity\\.", channels[i],"\\."), colnames(sites))]
    loc.cols <- colnames(sites)[grep(str_c("Localization\\.prob\\."), colnames(sites))]

    SN.col <- sym(str_c("SN.", channels[i]))
    filter.col <- sym(str_c("Is.missing.", channels[i]))
    reporter.corr.col <- sym(str_c("Reporter.intensity.corrected.", channels[i]))
    reporter.col <- sym(str_c("Reporter.intensity.", channels[i]))

    new.intensities[[i]] <-  sites %>%
      select(id, all_of(sn.cols), all_of(is.missing.cols), all_of(r.i.cols), all_of(r.i.c.cols), all_of(loc.cols)) %>%
      pivot_longer(cols = c(all_of(sn.cols), all_of(is.missing.cols), all_of(r.i.cols), all_of(r.i.c.cols), all_of(loc.cols)),
                   values_to = "Value",
                   names_to = 'old.name') %>%
      mutate(Experiment = str_match(old.name, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,4],
             Experiment = ifelse(is.na(Experiment), str_match(old.name, pattern = "(Localization\\.prob\\.)(.*)")[,3], Experiment),
             Column.type = str_c(str_match(old.name, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,2], ".", channels[i]),
             Column.type = ifelse(is.na(Column.type), 'Localization.prob', Column.type),
             Phosphosite = str_match(old.name, pattern = "(.*)\\.([0-9]*)\\.(.*)___([0-9]*)")[,5]) %>%
      replace_na(list(Phosphosite=1)) %>%
      select(-old.name) %>%
      pivot_wider(names_from = Column.type, values_from = Value) %>%
      replace_na(list(Localization.prob=0)) %>%
      group_by(Experiment, id) %>%
      summarise(!!SN.col := sum(!!SN.col),
                !!reporter.col := sum(!!reporter.col),
                !!reporter.corr.col := sum(!!reporter.corr.col),
                !!filter.col := sum(!!filter.col),
                Localization.prob = max(Localization.prob)) %>% ### because each peptide has only one localization prob
      ungroup() %>%
      filter(!!rlang::sym(filter.col) < 3)

    setTxtProgressBar(pb, i)
  }
  close(pb)

  updated <- new.intensities[[1]]
  for(i in 2:(length(channels))){
    updated <- updated %>%
      full_join(new.intensities[[i]], by = c('id', 'Experiment'))
  }

  keep <- c(colnames(updated)[grep("Reporter.intensity.", colnames(updated))],
            colnames(updated)[grep("Reporter.intensity.corrected,", colnames(updated))],
            colnames(updated)[grep("SN.", colnames(updated))],
            colnames(updated)[grep("Is.missing.", colnames(updated))])

  sites <- sites %>%
    left_join(updated, by = c('id')) %>%
    select(all_of(keep), Evidence.IDs, id, Proteins, Localization.prob, Protein.group.IDs, Experiment)

  ##############################################################################
  ### merge sites and evidence
  ##############################################################################
  sites <- sites %>%
    rename('Phospho.STY.site.IDs' = 'id') %>%
    separate_rows(Evidence.IDs, sep = ';', convert=TRUE)

  remove <- colnames(evidence)[grep("Reporter.intensity.", colnames(evidence))]
  evidence <- evidence %>%
    select(-all_of(remove), -all_of(keep), -Proteins, -Protein.group.IDs) %>%
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
