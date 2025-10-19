library(here)
library(tidyverse)
library(fgsea)
library(readxl)
library(SummarizedExperiment)
library(ggsurvfit)
library(survival)

here::i_am("code/suptab1.R")
output.path <- here("output/suptab1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


# Load Gene Sets ####
gene_sigs <- read_excel(here('data/gene_sets/b_t_cell_gene_sets.xlsx'))
bcgs.genes <- gene_sigs %>% filter(Signature == "b_cell")
bcgs <- list("BCGA" = bcgs.genes$Gene)


gene_sets <-
  (gene_sigs %>% group_by(Signature) %>% summarise(genes = list(unique(Gene))) %>% as.list(.))[[2]]
names(gene_sets) <-
  (gene_sigs %>% group_by(Signature) %>% summarise(genes = list(unique(Gene))) %>% as.list(.))[[1]]


C2_KEGG <- gmtPathways(here("data/gene_sets/msigdb_v2022.1.Hs_GMTs/c2.cp.kegg.v2022.1.Hs.symbols.gmt"))
C2_CGP <- gmtPathways(here("data/gene_sets/msigdb_v2022.1.Hs_GMTs/c2.cgp.v2022.1.Hs.symbols.gmt"))
C7_Immunesig <- gmtPathways(here("data/gene_sets/msigdb_v2022.1.Hs_GMTs/c7.immunesigdb.v2022.1.Hs.symbols.gmt"))
C8_CellType <- gmtPathways(here("data/gene_sets/msigdb_v2022.1.Hs_GMTs/c8.all.v2023.2.Hs.symbols.gmt"))


# Identify B_CELL related gene sets ####
identify_gene_sets <- function(x, pattern){
  names(x)[str_detect(names(x), pattern = pattern)]
  
}


#c2cgp_bgs <- identify_gene_sets(C2_CGP, "B_LYMPHOCYTE")
#c7_bgs <- identify_gene_sets(C7_Immunesig, "BCELL")
c8_bgs <- identify_gene_sets(C8_CellType, "_B_CELL")


#C2_CGP <- C2_CGP[c2cgp_bgs]
#C7_Immunesig <- C7_Immunesig[c7_bgs]
C8_CellType <- C8_CellType[c8_bgs]

#gene_sets <- c(gene_sets, C2_CGP, C7_Immunesig, C8_CellType)
gene_sets <- c(gene_sets,C8_CellType)

# Identify Overlaps ####
overlaps <- tibble::tibble(gene_set = rep(NA, length(gene_sets)), overlap = rep(NA, length(gene_sets)), size = rep(NA, length(gene_sets)))

for(i in seq_along(gene_sets)){
  gs = names(gene_sets[i])
  ol = length(BiocGenerics::intersect(bcgs[[1]], gene_sets[[i]]))
  size = length(gene_sets[[i]])
  overlaps$gene_set[i] <- gs
  overlaps$overlap[i] <- ol
  overlaps$size[i] <- size
}

head(overlaps)
ol = intersect(bcgs[[1]], c("BLK", "CD19"))


#Identify GSVA Correlation ####
library(GSVA)

se <- readRDS(here("output/00.0_data_prep/IMvigor010_tpm_se.rds"))
mat <- assays(se)[["tpm"]]

mat <- log2(mat+1)

gsvapar <- gsvaParam(mat, gene_sets, maxDiff=TRUE, kcdf="Gaussian")
results <- gsva(param = gsvapar)
cors <- cor(t(results))

gsva.cors.df <- data.frame(gsva.cor = cors[,"b_cell"], gene_set = names(cors[,"b_cell"]))
tables <- left_join(overlaps, gsva.cors.df)


# Add Hazard Ratio and Pvals ####
meta <- as_tibble(colData(se)) 
meta <- meta %>% mutate(OS_numeric = if_else(OS_event == "Alive", 0, 1))


tables$HR <- NA
tables$p.val <- NA

for(i in seq_along(1:nrow(results))){
  gs_val <- results[i,]
  gs_name <- rownames(results)[i]
  
  meta$gs_test = gs_val
  cox_model <- coxph(Surv(meta$OS_months, meta$OS_numeric) ~ meta$gs_test, data = meta)
  surv.stats <- summary(cox_model)  
  HR <- surv.stats$coefficients[2]
  p.val <- surv.stats$coefficients[5]
  
  tables$HR[i] <- HR
  tables$p.val[i] <- p.val
}
tables$fdr <- p.adjust(tables$p.val, method = "fdr")
tables %>% write_csv(here(output.path, "suptab1.csv"))

