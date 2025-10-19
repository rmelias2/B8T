library(here)
library(tidyverse)
library(consensusMIBC)
library(SummarizedExperiment)
source(here('code/function_stackedbar.R'))

here::i_am("code/supfig8.R")
output.path <- here("output/supfig8")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

#Compare B8T by Concensus and TCGA subgroups ####

se <- readRDS(here("output/00.0_data_prep/IMvigor010_tpm_se.rds"))
mat <- assays(se)[["tpm"]]

log2tpm <- log2(mat + 1)

consensus <- getConsensusClass(log2tpm, minCor = .2, gene_id = c("hgnc_symbol"))

data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))
data$consensus <- consensus$consensusClass

#Concensus Groups

cols = c("Hi/Hi" = "#d1d43f",
         "Hi/Lo" = "#2b628c",
         "Lo/Hi" = "#3fd485",
         "Lo/Lo" = "#523070")

g <- stacked_bar(df = data, x_variable = "consensus", group_variable = "B8T", stat = "chisq.test",
                 group.colors = cols,annotate_n = TRUE, x = "Concensus Classifier", 
                 y = "Frequency")

ggsave(plot = g, filename = "concensus_b8t.tiff", path = here(output.path), device = "tiff", width = 3.5, height = 3.5, unit = "in", dpi = 320, bg = "white")


#TCGA Groups
g <- stacked_bar(df = data, x_variable = "TCGA_subtype", stat = "chisq.test",
                 group_variable = "B8T", group.colors = cols,annotate_n = TRUE,
                 x = "TCGA Subtype", y = "Frequency")

ggsave(plot = g, filename = "tcga_b8t.tiff", path = here(output.path), device = "tiff", width = 3.5, height = 3.5, unit = "in", dpi = 320, bg = "white")

