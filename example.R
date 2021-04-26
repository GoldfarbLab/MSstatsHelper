library(tidyverse)
library(MSstatsHelper)

in_path <- "~/Box/CellBio-GoldfarbLab/Projects/Current/Blumer - Uveal Melanoma/txt-mp41_global/"
phos_path <- "~/Box/CellBio-GoldfarbLab/Projects/Current/Blumer - Uveal Melanoma/txt-mp41_IMAC/"
out_path <- "~/Documents/annotations/"

dir.create(file.path(out_path), showWarnings = FALSE)

evidence_path <- str_c(in_path, "/evidence.txt")
summary_path <- str_c(in_path, "/summary.txt")
mix_annot_path <- str_c(out_path, "/mixtureAnnotation.csv")
chan_annot_path <- str_c(out_path, "/channelAnnotation.csv")
annot_path <- str_c(out_path, "/annotation.csv")

# Initial file to fill out mixture and tech rep mixture
MQ_TMT_writeMixtureAnnotationFile(evidence_path, summary_path, mix_annot_path)
# Second file to fill out channel, condition, bioReiplicate for each mixture
MQ_TMT_writeChannelAnnotationFile(evidence_path, mix_annot_path, chan_annot_path)
# Merge channel and run info
MQ_TMT_writeFullAnnotationFile(mix_annot_path, chan_annot_path, annot_path)

evidence <- read_tsv(evidence_path, guess_max = 100000)
proteinGroups <- read_tsv(str_c(in_path, "/proteinGroups.txt"), guess_max = 100000)
sites <- read_tsv(str_c(phos_path, "/Phospho (STY)Sites.txt"), guess_max = 100000)

# order protein IDs and gene names
evidence <- sortIdentifiers(evidence)
sites <- sortIdentifiers(sites)

# how many Protein IDs match between phosphosites and global evidence?
localized.sites <- sites %>% filter(`Localization prob` >= 0.75)

# formatting
formatted_data <- prepareForMSstats(localized.sites, evidence)
new_evidence <- global_evidence[[1]]
new_sites <- global_evidence[[2]]



matching.localized.sites <- new_sites %>% filter(`Sorted Proteins` %in% new_evidence$`Sorted Proteins`)
matching.localized.sites.gene <- new_sites %>% filter(`Sorted Gene names` %in% new_evidence$`Sorted Gene names`)

not.matching.localized.sites.gene <- new_sites %>% filter(!(`Gene names` %in% new_evidence$`Gene names`))




#mapping <- read_tsv("~/Box/CellBio-GoldfarbLab/Data/Annotations/uniprot_mapping.tsv.zip")
