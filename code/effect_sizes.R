library(here)
library(tidyverse)
library(ggsurvfit)


here::i_am("code/effect_sizes.R")
output.path <- here("output/effect_sizes")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


#Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

#Reformat variables
data <- data %>% mutate(B8T_hihivselse = if_else(B8T == "Hi/Hi", "Hi/Hi", "non-Hi/Hi"),
                        ARM = if_else(ARM == "Observation", "Observation", "Atezo"),
                        IFN_group = if_else(IFN_gamma > median(IFN_gamma), "high", "low"))

data$OS_followup <- abs(data$OS_numeric - 1)


## Effect sizes for Median Follow Up
surv.object <- survfit2(Surv(df[["OS_months"]], df[["OS_followup"]]) ~ 1, data = df) 
surv.object



atezo <- data %>% filter(ARM == "Atezo", ctDNA_call_C1D1 == "Positive")
surv.object <- survfit2(Surv(atezo[["OS_months"]], atezo[["OS_numeric"]]) ~ B8T, data = atezo) 
surv.object

obs <- data %>% filter(ARM == "Observation", ctDNA_call_C1D1 == "Positive")
surv.object <- survfit2(Surv(obs[["OS_months"]], obs[["OS_numeric"]]) ~ B8T, data = obs) 
surv.object


##B8T Hi/Lo survival difference: 
Hilo <- data %>% filter(B8T == "Hi/Lo")
surv.object <- survfit2(Surv(Hilo[["OS_months"]], Hilo[["OS_numeric"]]) ~ ARM, data = Hilo) 
surv.object

## Effect Size for B8T 

km_curve(Data = df, groups = "b_cell_bin", group_label = "BGS", subset_label = "Atezo, ctDNA Positive", plot_colors = colors)