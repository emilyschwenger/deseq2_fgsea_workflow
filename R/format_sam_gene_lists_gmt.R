## Format Sam's CHIP and CutnTag gene lists into gmt file for GSEA ##

# Load libraries
library(tidyverse)
library(here)

all <- read_tsv(here("data/msigdb_genesets/custom/chip_chip_cutntag_all.gmt"), col_names = F)
xl <- readxl::read_xlsx(here("data/msigdb_genesets/custom/URE summary lists Chip for RNA seq anal.xlsx"), sheet = "LOST gene lists")

xl <- xl %>%
  distinct %>%
  mutate(CHIP_PU1_LOST_URE = toupper(`All Lost peak genes`),
         CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE = toupper(`->Promoter/Intron/Exon/TTS`),
         CHIP_PU1_PROMOTER_LOST_URE = toupper(`-> PROMOTER only`)) %>%
  dplyr::select(CHIP_PU1_LOST_URE, CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE, CHIP_PU1_PROMOTER_LOST_URE)

xl <- rbind(c("> All Lost peak genes", "> Promoter Intron Exon TTS", "> Promoter only"), xl)
xl <- rbind(colnames(xl), xl)

v <- xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE[!is.na(xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE)]
v <- v[v==""]

v <- xl$CHIP_PU1_PROMOTER_LOST_URE[!is.na(xl$CHIP_PU1_PROMOTER_LOST_URE)]
v <- v[v==""]



v <- data.frame(xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE[!is.na(xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE)])

write(paste(as.character(v), collapse="\t"),
      here("data/msigdb_genesets/custom/chip_cutntag_promoter_exon_intron_tts_ure.gmt"), append="TRUE")

# Sanity checks
xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE[xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE==""]

all <- rbind(all, c(xl$CHIP_PU1_LOST_URE[!is.na(xl$CHIP_PU1_LOST_URE)],
                    xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE[!is.na(xl$CHIP_PU1_PROMOTER_EXON_INTRON_TTS_LOST_URE)],
                    xl$CHIP_PU1_PROMOTER_LOST_URE[!is.na(xl$CHIP_PU1_PROMOTER_LOST_URE)]))
dim(all)
dim(all %>%
  na.omit() )

write_tsv(all, here("data/msigdb_genesets/custom/chip_chip_cutntag_all_ure3.gmt"), na = "", col_names = F)
write_tsv(all[11,], here("data/msigdb_genesets/custom/chip_cutntag_ure_only.gmt"), na = "", col_names = F)
