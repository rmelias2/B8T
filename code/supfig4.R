library(here)
library(tidyverse)
library(ggsurvfit)


here::i_am("code/supfig4.R")
output.path <- here("output/supfig4")


#Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

#Reformat variables
data <- data %>% mutate(B8T_hihivselse = if_else(B8T == "Hi/Hi", "Hi/Hi", "non-Hi/Hi"),
                        ARM = if_else(ARM == "Observation", "Obs", "Atezo"),
                        IFN_group = if_else(IFN_gamma > median(IFN_gamma), "high", "low"))


# Plot B8T Hi/Hi by Treatment Arm
colors = c("Atezo" = "red", "Obs"= "blue")
df <- data %>% filter(B8T == "Hi/Hi",ctDNA_call_C1D1 == "Negative")

km_curve(Data = df, groups = "ARM", group_label = "Arm", subset_label = "Hi/Hi, ctDNA Negative", plot_colors = colors)

df <- data %>% filter(B8T == "Hi/Hi")

km_curve(Data = df, groups = "ARM", group_label = "Arm", subset_label = "Hi/Hi", plot_colors = colors)