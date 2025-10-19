library(matrixStats)
library(tidyverse)
library(here)
library(clusterProfiler)
library(SummarizedExperiment)

here::i_am("code/00.0_data_prep.R")
output.path <- here("output/00.0_data_prep")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


## TPM data and clinical annotation for IMmVigor010 were accessed from EGA, EGA Study ID: EGAS00001004997 ----
tpm <- read_csv(here("data/IMvigor010/tpm.csv"))
clinical1 <- read_csv(here("data/IMvigor010/wo29636_clinical_20210301.csv"))
clinical2 <- read_csv(here("data/IMvigor010/wo29636_clinical_20210423.csv"))
id_key <- read_tsv(here("data/IMvigor010/Run_Sample_meta_info.map"), col_names = FALSE)


## Link SampleID to PatientID
id_key <- id_key %>% dplyr::rename(RNA_ID = X1, clin = X8) %>% select(RNA_ID, clin)
id_key <- id_key %>% tidyr::separate(clin, into = "Pt_ID", sep = ";")
id_key <- id_key %>% tidyr::separate(Pt_ID, into = c("del", "PTID"), sep = "=")
id_key <- id_key %>% dplyr::select(-del)


## Join Data
clinical1 <- clinical1 %>% tidyr::separate(UNI_ID, into = "pt_id", sep = 12)
clinical1 <- clinical1 %>% mutate(PTID = str_c("PAT-", pt_id)) %>% select(PTID, everything())
clinical2 <- clinical2 %>% dplyr::rename(PTID = PATIENT_ID)
meta <- left_join(id_key, clinical1)
meta <- left_join(meta, clinical2)
meta <- meta %>% select(-pt_id)

## Convert Gene IDs to Symbols in TPM File
tpm <- tpm %>% dplyr::rename(GeneID = `Unnamed: 0`)
gen_ids <- tibble(GeneID = tpm$GeneID) %>% 
  tidyr::separate(GeneID, into = c("discard", "ENTREZID"), 
           sep = ":", remove = FALSE)

eg = bitr(gen_ids$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
gen_ids <- left_join(gen_ids, eg) %>% select(-discard)

dup_list <- (gen_ids %>% count(SYMBOL) %>% filter(n>1))$SYMBOL
dup_GeneID <- gen_ids %>% filter(SYMBOL %in% dup_list) %>% select(GeneID)
tpm <- tpm %>% filter(!GeneID %in% dup_GeneID$GeneID)
gen_ids <- gen_ids %>% filter(GeneID %in% tpm$GeneID)
which(!gen_ids$GeneID == tpm$GeneID)

mat <- as.matrix(
  tpm %>% select(2:ncol(tpm))
)
rownames(mat) <- gen_ids$SYMBOL


meta$RNA_ID <- factor(meta$RNA_ID, levels = c(colnames(mat)))
meta <- arrange(meta, RNA_ID)

# Save Summarized Experiment Object 
se <- SummarizedExperiment(assays = list("tpm" = mat), rowData = gen_ids, colData = meta)
saveRDS(se, here(output.path, "IMvigor010_tpm_se.rds"))

# Save TSV 
tpm_df <- as.data.frame(assays(se)[["tpm"]])
tpm_df$Gene_ID <- rownames(se)
tpm_df <- tpm_df %>% select(Gene_ID, everything())
write_tsv(tpm_df, here(output.path, "IMvigor010_tpm_dataframe.txt"))

# ABACUS DATA ####
# Accessed from the EGA, Dataset ID: EGAD00001006199
# Refers to the following paper (https://www.nature.com/articles/s41591-019-0628-7)

log2.tpm <- read_csv(here("data/ABACUS/abacus.tpm.log2.csv"))
geneIds <-  read_csv(here("data/ABACUS/abacus.RNAseq.genes.csv"))

#Convert geneIDs to Symbols:

## Remove genes without a Symbol or with duplicate symbols
geneIds <- geneIds %>% filter(!is.na(symbol)) %>% 
  filter(!duplicated(symbol)  & !duplicated(symbol, fromLast = TRUE))

## Subset TPM matrix for GeneIDs contained in the above dataframe
log2.tpm <- log2.tpm %>% filter(Gene %in% geneIds$full_gene_id)

## Verify order of GeneIDs and Symbols are Equivallent in the two dataframes
which(!geneIds$full_gene_id == log2.tpm$Gene) 

## Concert to matrix
mat <- as.matrix(
  log2.tpm %>% select(2:ncol(log2.tpm))
)
rownames(mat) <- geneIds$symbol

## Save
saveRDS(mat, here(output.path, "abacus.log2.tpm.clean.RDS"))

#Clinical Annotations: 
abacus_clinical <- read_csv(here("data/ABACUS/anon_abacus.RNAseq.pdata.public.csv"))

#Clean ID variable names
abacus_clinical <- abacus_clinical %>% dplyr::rename(RNAID = alias, ptID = subjectId)
write_csv(abacus_clinical, here(output.path,"abacus.meta.clean.csv"))


sessionInfo()


