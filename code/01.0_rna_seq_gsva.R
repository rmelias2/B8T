library(matrixStats)
library(tidyverse)
library(here)
library(GSVA)
library(SummarizedExperiment)
library(readxl)

here::i_am("code/01.0_rna_seq_gsva.R")
output.path <- here("output/01.0_rna_seq_gsva")
if (!dir.exists(here(output.path))) dir.create(here(output.path))



#Read in SE Object with TPM

se <- readRDS(here("output/00.0_data_prep/IMvigor010_tpm_se.rds"))
mat <- assays(se)[["tpm"]]

## Read in Genesets
gene_sigs <- read_excel(here('data/gene_sets/b_t_cell_gene_sets.xlsx'))
gene_sets <-
  (gene_sigs %>% group_by(Signature) %>% summarise(genes = list(unique(Gene))) %>% as.list(.))[[2]]
names(gene_sets) <-
  (gene_sigs %>% group_by(Signature) %>% summarise(genes = list(unique(Gene))) %>% as.list(.))[[1]]


## Log Transform and perform GSVA 
mat <- log2(mat+1)
results <- gsva(mat, gene_sets)

## Adding results to metadata file
data <- as.data.frame(t(results))
se$b_cell <- data$b_cell
se$t_cell <- data$t_cell
se$IFN_gamma <- data$IFN_gamma

meta <- as_tibble(colData(se)) %>% 
  mutate(
    t_cell_bin = if_else(t_cell > 0, "High", "Low"),
    b_cell_bin = if_else(b_cell > 0, "High", "Low"),
    B8T = str_c(b_cell_bin, "/", t_cell_bin),
    OS_numeric = if_else(OS_event == "Alive", 0, 1),
    DFS_numeric = if_else(DFS_event == "UC Free", 0, 1)
  ) 

meta <- meta %>% mutate(B8T = case_when(
  B8T == "High/High" ~ "Hi/Hi",
  B8T == "High/Low" ~ "Hi/Lo",
  B8T == "Low/High" ~ "Lo/Hi",
  B8T == "Low/Low" ~ "Lo/Lo"
))
meta <- meta %>% mutate(ctDNA_call_C1D1 = if_else(is.na(ctDNA_call_C1D1), "NA", ctDNA_call_C1D1))

write_tsv(meta, here(output.path, "gsva_results_log2.txt"))


#ABACUS ####
log2.tpm <- readRDS(here("output/00.0_data_prep/abacus.log2.tpm.clean.RDS"))
meta <- read_csv(here("output/00.0_data_prep/abacus.meta.clean.csv"))


### Run GSVA
gene_sigs <- read_excel(here('data/gene_sets/b_t_cell_gene_sets.xlsx'))
gene_sets <-
  (gene_sigs %>% group_by(Signature) %>% summarise(genes = list(unique(Gene))) %>% as.list(.))[[2]]
names(gene_sets) <-
  (gene_sigs %>% group_by(Signature) %>% summarise(genes = list(unique(Gene))) %>% as.list(.))[[1]]

results <- gsva(mat, gene_sets)


## Generate output: 
data <- as.data.frame(t(results))
data$RNAID <- rownames(t(results))
data <- left_join(data,meta)
data <- as_tibble(data)

data <- data %>% mutate(
  t_cell_bin = if_else(t_cell > 0, "Hi", "Lo"),
  b_cell_bin = if_else(b_cell > 0, "Hi", "Lo"),
  B8T = str_c(b_cell_bin, "/", t_cell_bin))

write_tsv(data, here(output.path, "ABACUS_gsva_results.txt"))


sessionInfo()

