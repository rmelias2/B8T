
library(here)
library(tidyverse)
library(ggsurvfit)
library(survival)


here::i_am("code/table1.R")
output.path <- here("output/table1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))


## Clean Variables ----

data$survival_obj <- Surv(data$OS_months, data$OS_numeric)
data$Baseline_ECOG <- factor(data$Baseline_ECOG)
data <- data %>% mutate(arm_label = if_else(ARM == "Observation", "Obs", "Atezo"))
data$arm_label <- factor(data$arm_label, levels = c("Obs", "Atezo"))
data$TMB_status <- factor(data$TMB_status, levels = c("TMB_low", "TMB_high"))
data <- data %>% mutate(ctDNA_call_C1D1 = if_else(is.na(ctDNA_call_C1D1), "Missing", ctDNA_call_C1D1))
data$B8T <- factor(data$B8T, levels = c("Lo/Lo", "Lo/Hi", "Hi/Lo", "Hi/Hi"))
data$ctDNA_call_C1D1 <- factor(data$ctDNA_call_C1D1, levels = c("Negative", "Missing", "Positive"))


#### Subset by B8T cohort ----

B8T_groups = unique(data$B8T)

for(i in seq_along(1:length(B8T_groups))){
  label <- str_replace(B8T_groups[i], "/", "_")
  df <- data %>% filter(B8T == B8T_groups[i])
  
  vars_for_mva <- c('tumor_stage', 'nodal_status', 'PDL1_status','arm_label','ctDNA_call_C1D1', 'TMB_status')
  
  multivariate_formula <- as.formula(paste('survival_obj ~', paste(vars_for_mva, collapse = ' + ')))
  multivariate_cox_model <- coxph(multivariate_formula, data = df)
  multivariate_summary <- summary(multivariate_cox_model)
  
  t1 <- multivariate_summary$coefficients
  t2 <- multivariate_summary$conf.int
  
  t3 <- bind_cols(t1,t2)
  t3$var = rownames(t2)
  t3 <- t3 %>% select(var, `exp(coef)...2`, `lower .95`, `upper .95`,`Pr(>|z|)`)
  t3 %>% write_csv(here(output.path, paste0(label,".csv")))
  
}
