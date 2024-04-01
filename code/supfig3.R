library(here)
library(tidyverse)
library(ggsurvfit)


here::i_am("code/supfig3.R")
output.path <- here("output/supfig3")


#Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

#Reformat variables
data <- data %>% mutate(B8T_hihivselse = if_else(B8T == "Hi/Hi", "Hi/Hi", "non-Hi/Hi"),
                        ARM = if_else(ARM == "Observation", "Obs", "Atezo"),
                        IFN_group = if_else(IFN_gamma > median(IFN_gamma), "high", "low"))



# Plot B8T by ctDNA status and Treatment Arm ####
colors = c("Hi/Hi" = "#d1d43f", "Hi/Lo"= "#2b628c", "Lo/Hi" = "#3fd485", "Lo/Lo" = "#523070")

#ctDNA positive, Atezo 
df <- data %>% filter(ctDNA_call_C1D1 == "Negative", ARM == "Atezo")

km_curve(Data = df, groups = "B8T", group_label = "B8T", subset_label = "Atezo, ctDNA Negative", plot_colors = colors, width = 3.5)

df <- data %>% filter(ctDNA_call_C1D1 == "Negative", ARM == "Obs")

km_curve(Data = df, groups = "B8T", group_label = "B8T", subset_label = "Obs, ctDNA Negative", plot_colors = colors, width = 3.5)

