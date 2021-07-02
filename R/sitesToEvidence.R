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

  ##############################################################################
  ### sum together phosphosites into one intensity column
  ##############################################################################
  cols <- as_tibble(str_match(colnames(sites), 'Reporter.intensity.corrected.([0-9]*)___([0-9]*)')) %>%
    mutate(output_column_name = str_c("Reporter.intensity.corrected.",V2)) %>%
    drop_na()

  for(i in 1:length(unique(cols$V2))){
    small_cols <- cols %>% filter(V2 == as.character(i))
    sites <- sites %>% mutate(
      !!as.character(small_cols[1,4]) :=
        !!sym(as.character(small_cols[1,1])) +
        !!sym(as.character(small_cols[2,1])) +
        !!sym(as.character(small_cols[3,1]))
    )
  }

  remove <- colnames(sites)[grep("___", colnames(sites))]
  sites <- sites %>% select(-all_of(remove))

  keep <- colnames(sites)[grep("Reporter.intensity.corrected.", colnames(sites))]
  sites <- sites %>% select(all_of(keep), Evidence.IDs, id, Proteins)

  ##############################################################################
  ### merge sites and evidence
  ##############################################################################
  sites <- sites %>%
    rename('Phospho..STY..site.IDs' = 'id') %>%
    mutate(Evidence.IDs = str_split(Evidence.IDs, ';')) %>%
    unnest(Evidence.IDs) %>%
    mutate(Evidence.IDs = as.integer(Evidence.IDs))

  remove <- colnames(evidence)[grep("Reporter.intensity.", colnames(evidence))]
  evidence <- evidence %>%
    select(-all_of(remove), -Proteins) %>%
    # select('Modified.sequence', 'Raw.file', 'Phospho..STY..site.IDs', 'id') %>%
    rename('Evidence.IDs' = 'id')

  evidence <- evidence %>%
    mutate(Phospho..STY..site.IDs = str_split(Phospho..STY..site.IDs, ';')) %>%
    unnest(Phospho..STY..site.IDs) %>%
    mutate(Phospho..STY..site.IDs = as.integer(Phospho..STY..site.IDs))

  input <- sites %>% left_join(evidence, by = c("Phospho..STY..site.IDs", "Evidence.IDs"))
  input <- input %>%
    group_by(Phospho..STY..site.IDs) %>%
    arrange(PEP, .by_group = TRUE) %>%
    slice(1)

  input <- input %>%
    rename(id = Phospho..STY..site.IDs) %>%
    select(-Evidence.IDs) %>%
    unique()

  ##############################################################################
  ### update Evidence ID in the proteinGroups table
  ##############################################################################
  map <- input %>%
    select(Protein.group.IDs, id) %>%
    mutate(Protein.group.IDs = str_split(Protein.group.IDs, ';')) %>%
    unnest(Protein.group.IDs) %>%
    mutate(Protein.group.IDs = as.numeric(Protein.group.IDs)) %>%
    unique()

  map <- map %>%
    group_by(Protein.group.IDs) %>%
    summarize(Evidence.IDs = str_c(id, collapse = ';'))

  proteinGroups <- proteinGroups %>% select(-Evidence.IDs)
  proteinGroups <- proteinGroups %>%
    left_join(map, by = c('id' = 'Protein.group.IDs'))

  return(list("Sites" = input, "proteinGroups" = proteinGroups))
}
