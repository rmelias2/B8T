library(here)
library(tidyverse)
library(Gmisc)

here::i_am("code/tab1.R")
output.path <- here("output/tab1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


#' ## Load Data

data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))


df <- data %>% 
  mutate(tumor_stage = if_else(tumor_stage == "PT3/PT4", "PT3/PT4", "PT1/PT2"),
         ARM = if_else(ARM == "Observation", "Observation", "Atezo"), 
         Baseline_ECOG = as.character(Baseline_ECOG))

df <- df %>% 
  set_column_labels(
    tumor_stage = "T stage",
    nodal_status = "Nodal Status",
    Baseline_ECOG = "ECOG",
    prior_neoadjuvant_chemotherapy = "Prior NeoAdj Chemo",
    SEX = "Sex",
    ctDNA_call_C1D1 = "ctDNA C1D1",
    ctDNA_call_C3D1 = "ctDNA C3D1",
    ctDNA_clearance_C1D1_to_C3D1 = "ctDNA clearance",
    PDL1_status = "PDL1",
    TMB_status = "TMB",
    OS_event = "Survival Status at Last Follow Up")



### Descriptive Stats pooled cohort
t <- df %>% getDescriptionStatsBy(tumor_stage, nodal_status, Baseline_ECOG, prior_neoadjuvant_chemotherapy, SEX, ctDNA_call_C1D1, PDL1_status, TMB_status ,
                                  B8T, OS_event,
                                  by = ARM, add_total_col = TRUE,header_count = TRUE, 
                                  show_all_values = TRUE, statistics = list(proportion = getPvalFisher)) %>% htmlTable()

htmltools::save_html(htmlTable(t), file = here(output.path, paste0("IMVIGOR010_cohort_overview", ".html")))




sessionInfo()


