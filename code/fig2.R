library(here)
library(tidyverse)
library(ggsurvfit)
library(consensusMIBC)

source(here('code/function_kmcurve.R'))
source(here('code/function_stackedbar.R'))

here::i_am("code/fig2.R")
output.path <- here("output/fig2")


#Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

#Reformat variables
data <- data %>% mutate(B8T_hihivselse = if_else(B8T == "Hi/Hi", "Hi/Hi", "non-Hi/Hi"),
                        ARM = if_else(ARM == "Observation", "Obs", "Atezo"),
                        IFN_group = if_else(IFN_gamma > median(IFN_gamma), "high", "low"))



# Plot B8T by ctDNA status and Treatment Arm ####
colors = c("Hi/Hi" = "#d1d43f", "Hi/Lo"= "#2b628c", "Lo/Hi" = "#3fd485", "Lo/Lo" = "#523070")

#ctDNA positive, Atezo 
df <- data %>% filter(ctDNA_call_C1D1 == "Positive", ARM == "Atezo")

km_curve(Data = df, groups = "B8T", group_label = "B8T", subset_label = "Atezo, ctDNA Positive", plot_colors = colors, width = 3.5)

df <- data %>% filter(ctDNA_call_C1D1 == "Positive", ARM == "Obs")

km_curve(Data = df, groups = "B8T", group_label = "B8T", subset_label = "Obs, ctDNA Positive", plot_colors = colors, width = 3.5)


# Plot B8T Hi/Hi by Treatment Arm
colors = c("Atezo" = "red", "Obs"= "blue")
df <- data %>% filter(ctDNA_call_C1D1 == "Positive", B8T == "Hi/Hi")

km_curve(Data = df, groups = "ARM", group_label = "Arm", subset_label = "Hi/Hi, ctDNA Positive", plot_colors = colors)

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


sessionInfo()




