library(here)
library(tidyverse)
library(ggsurvfit)


here::i_am("code/fig1.R")
output.path <- here("output/fig1")


#Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

#Reformat variables
data <- data %>% mutate(B8T_hihivselse = if_else(B8T == "Hi/Hi", "Hi/Hi", "non-Hi/Hi"),
                        ARM = if_else(ARM == "Observation", "Observation", "Atezo"),
                        IFN_group = if_else(IFN_gamma > median(IFN_gamma), "high", "low"))
                    

source(here('code/function_kmcurve.R'))


# Plot BCGS Hi vs Low by ctDNA status and Treatment Arm ####
colors = c("High"= "#d1d43f", "Low" = "#523070")

#ctDNA positive, Atezo 
df <- data %>% filter(ctDNA_call_C1D1 == "Positive", ARM == "Atezo")

km_curve(Data = df, groups = "b_cell_bin", group_label = "BCGS", subset_label = "Atezo, ctDNA Positive", plot_colors = colors)

#ctDNA positive, Observation 
df <- data %>% filter(ctDNA_call_C1D1 == "Positive", ARM == "Observation")
km_curve(Data = df, groups = "b_cell_bin", group_label = "BCGS", subset_label = "Obs, ctDNA Positive", plot_colors = colors)

#ctDNA negative, Observation 
df <- data %>% filter(ctDNA_call_C1D1 == "Negative", ARM == "Observation")
km_curve(Data = df, groups = "b_cell_bin", group_label = "BCGS", subset_label = "Obs, ctDNA Negative", plot_colors = colors)

#ctDNA negative, Atezo 
df <- data %>% filter(ctDNA_call_C1D1 == "Negative", ARM == "Atezo")
km_curve(Data = df, groups = "b_cell_bin", group_label = "BCGS", subset_label = "Atezo, ctDNA Negative", plot_colors = colors)

