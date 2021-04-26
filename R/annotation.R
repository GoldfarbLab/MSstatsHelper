#' Writes first annotation file to fill out mixture and tech rep mixture.
#' Output is based on the header of the evidence file and the full summary file
#'
#' @param evidence_path input path to MaxQuant evidence file
#' @param summary_path input path to MaxQuant summary file
#' @param out_path output path to mixture annotation file
#'
#' @export
#' @examples
#' \dontrun{
#' MQ_TMT_writeMixtureAnnotationFile("evidence.txt", "summary.txt", "mixtureAnnotation.csv")
#' }
#'
#' @import dplyr
#' @import readr
#'
MQ_TMT_writeMixtureAnnotationFile <- function(evidence_path, summary_path, out_path) {
  evidence <- read_tsv(evidence_path, n_max=1)
  summary <- read_tsv(summary_path)

  mixtureAnnotation <- summary %>%
    filter(`Raw file` != "Total") %>%
    select(Run = `Raw file`, Fraction) %>%
    mutate(Mixture = "",
           TechRepMixture = "")

  write_csv(mixtureAnnotation, out_path)
}

#' Writes second annotation file to fill out channel, condition, bioReiplicate for each mixture.
#'
#' @param evidence_path input path to MaxQuant evidence file
#' @param mix_annot_path input path to mixture annotation file
#' @param out_path output path to channel annotation file
#'
#' @export
#' @examples
#' \dontrun{
#' MQ_TMT_writeChannelAnnotationFile("evidence.txt", "mixtureAnnotation.csv", "channelAnnotation.csv")
#' }
#'
#' @import dplyr
#' @import readr
#' @import stringr
#'
MQ_TMT_writeChannelAnnotationFile <- function(evidence_path, mix_annot_path, out_path) {
  evidence <- read_tsv(evidence_path, n_max=1)
  mixtureAnnotation <- read_csv(mix_annot_path)

  col_pattern <- "Reporter intensity [:digit:]+"

  num_channels <- sum(str_detect(col_pattern, colnames(evidence)))
  channels <- tibble(Channel = str_c("channel.", seq(0, num_channels-1)))

  channelAnnotation <- mixtureAnnotation %>%
    select(Mixture) %>%
    unique() %>%
    crossing(channels) %>%
    mutate(Condition="",
           BioReplicate="")

  write_csv(channelAnnotation, out_path)
}

#' Writes final annotation file that merges channel and run info
#'
#' @param mix_annot_path input path to mixture annotation file
#' @param chan_annot_path input path to channel annotation file
#' @param out_path output path to final annotation file
#'
#' @export
#' @examples
#' \dontrun{
#' MQ_TMT_writeFullAnnotationFile("mixtureAnnotation.csv", "channelAnnotation.csv", "annotation.csv")
#' }
#'
#' @import dplyr
#' @import readr
#'
MQ_TMT_writeFullAnnotationFile <- function(mix_annot_path, chan_annot_path, out_path) {
  channelAnnotation <- read_csv(chan_annot_path)
  mixtureAnnotation <- read_csv(mix_annot_path)

  fullAnnotation <- channelAnnotation %>%
    full_join(mixtureAnnotation, by=c("Mixture")) %>%
    relocate(Run, Fraction, Mixture, TechRepMixture, Channel, Condition, BioReplicate)

  write_csv(fullAnnotation, out_path)
}
