library(here)
library(tidyverse)

here::i_am("code/supfig1.R")
output.path <- here("output/supfig1")
if (!dir.exists(here(output.path))) dir.create(here(output.path))


## Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))


plot_and_save_density <- function(data, xvar, outputname, xlabname){
  g <- ggplot(data, aes(x = .data[[xvar]])) + geom_density(color="darkblue", fill="darkblue") + ylab("Density") + xlab(xlabname)
  g <- g + 
    theme_bw() +
    theme(text = element_text(size = 12)) 
  
  ggsave(filename = paste0(outputname,".tiff"), path = output.path, plot = g, device = "tiff", width = 2.5, height = 2.5, units = "in", dpi = 320)
}


#Sup Figure 1
# All Patients ####
#BCGS 
plot_and_save_density(data = data, x = "b_cell", outputname = "b_cell_densityplot", xlabname = "BSGS")

#CD8GS 
plot_and_save_density(data = data, x = "t_cell", outputname = "t_cell_densityplot", xlabname = "CD8TGS")

#ctDNA Positive Subset ####
df <- data %>% filter(ctDNA_call_C1D1 == "Positive")

#BCGS
plot_and_save_density(data = df, x = "b_cell", outputname = "ctdnaPos_b_cell_densityplot", xlabname = "BSGS")

#CD8GS
plot_and_save_density(data = df, x = "t_cell", outputname = "ctdnaPos_t_cell_densityplot", xlabname = "CD8TGS")

#ctDNA Negative Subset ####
df <- data %>% filter(ctDNA_call_C1D1 == "Negative")

#BCGS
plot_and_save_density(data = df, x = "b_cell", outputname = "ctdnaNegative_b_cell_densityplot", xlabname = "BSGS")

#CD8GS
plot_and_save_density(data = df, x = "t_cell", outputname = "ctdnaNegative_t_cell_densityplot", xlabname = "CD8TGS")

sessionInfo()


