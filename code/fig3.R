library(here)
library(tidyverse)
library(ggsurvfit)


here::i_am("code/fig3.R")
output.path <- here("output/fig3")


#Load Data
data <- read_tsv(here("output/01.0_rna_seq_gsva/gsva_results_log2.txt"))

#Reformat variables
data <- data %>% mutate(B8T_hilovselse = if_else(B8T == "Hi/Lo", "Hi/Lo", "non-Hi/Lo"),
                        ARM = if_else(ARM == "Observation", "Obs", "Atezo"),
                        IFN_group = if_else(IFN_gamma > median(IFN_gamma), "high", "low"))


# Plot B8T Hi/Hi by Treatment Arm
colors = c("Atezo" = "red", "Obs"= "blue")
df <- data %>% filter(B8T == "Hi/Lo")

km_curve(Data = df, groups = "ARM", group_label = "Arm", subset_label = "Hi/Lo", plot_colors = colors, pval.x = 10)

df <- data %>% filter(B8T == "Hi/Lo",ctDNA_call_C1D1 == "Negative")
km_curve(Data = df, groups = "ARM", group_label = "Arm", subset_label = "Hi/Lo, ctDNA Negative", plot_colors = colors)

df <- data %>% filter(B8T == "Hi/Lo",ctDNA_call_C1D1 == "Positive")
km_curve(Data = df, groups = "ARM", group_label = "Arm", subset_label = "Hi/Lo, ctDNA Positive", plot_colors = colors)


# Plot CIBERSORT B Cells

source(here('code/function_boxplot.R'))
cibersortx <- read_tsv(here('output/02.0_cibersortx/cibersortx.tsv'))

#Plot theme
theme <- theme_minimal() + 
  theme(
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_blank(),
    text = element_text(family = "sans", color = "black"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))

colors = c("Hi/Hi" = "#d1d43f", "Hi/Lo"= "#2b628c", "Lo/Hi" = "#3fd485", "Lo/Lo" = "#523070")
    
comparisons <- list( c("Hi/Hi", "Hi/Lo"), c("Hi/Lo", "Lo/Hi"), c("Lo/Hi", "Lo/Lo"), c("Hi/Hi", "Lo/Hi"),
                        c("Hi/Lo", "Lo/Lo"), c("Hi/Hi", "Lo/Lo"))

g <- ggplot_samplelvl_boxplot(data = cibersortx, x_variable = "B8T", y_variable = "B cells memory", 
                         y.lab = "Estimated Proportion", group.colors = colors, default.theme = theme, comparisons = comparisons)
g <- g + ylim(0,0.5)

ggsave(filename = paste0("mem_bcells",".tiff"), path = output.path, plot = g, device = "tiff", width = 3, height = 3, units = "in", dpi = 320)

g <- ggplot_samplelvl_boxplot(data = cibersortx, x_variable = "B8T", y_variable = "B cells naive", 
                              y.lab = "Estimated Proportion", group.colors = colors, default.theme = theme, comparisons = comparisons)
g <- g + ylim(0,0.5)
ggsave(filename = paste0("naive_bcells",".tiff"), path = output.path, plot = g, device = "tiff", width = 3, height = 3, units = "in", dpi = 320)

