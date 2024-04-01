library(here)
library(tidyverse)


here::i_am("code/02.0_cibersortx.R")
output.path <- here("output/02.0_cibersortx")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


# CIBERSORT performed using : https://cibersortx.stanford.edu/ on tpm data with default settings and LM22 signature matrix

cibersortx <- read_tsv(here('data/cibersort_output/CIBERSORTx_output.txt'))
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

cibersortx <- cibersortx %>% select(Mixture:Neutrophils)
cibersortx <- cibersortx %>% dplyr::rename(RNA_ID = Mixture)
annotated_output <- left_join(data, cibersortx)

write_tsv(annotated_output, here(output.path,"cibersortx.tsv"))
