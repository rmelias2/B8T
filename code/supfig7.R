library(here)
library(tidyverse)
library(ggsurvfit)


here::i_am("code/supfig7.R")
output.path <- here("output/supfig7")



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

of_interest <- c("Macrophages M0", "Macrophages M1","Macrophages M2", "Monocytes", "NK cells activated", "NK cells resting",
                 "Plasma cells")

for(i in of_interest){
  g <- ggplot_samplelvl_boxplot(data = cibersortx, x_variable = "B8T", y_variable =  i, 
                                y.lab = "Estimated Proportion", group.colors = colors, default.theme = theme, comparisons = comparisons)
  
  ggsave(filename = paste0(i,".tiff"), path = output.path, plot = g, device = "tiff", width = 3, height = 3, units = "in", dpi = 320)
  
}

