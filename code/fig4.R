library(here)
library(tidyverse)


here::i_am("code/fig4.R")
output.path <- here("output/fig4")

#Load Data:
results <- read_tsv(here("output/01.0_rna_seq_gsva/ABACUS_gsva_results.txt"))

#Clean variables
results <- results %>% mutate(response = case_when(PCR == "Yes" ~ "PCR", 
                                                   MPR == "Yes" ~"MPR", .default = "NR"))
results <- results %>% mutate(response_bin = case_when(RELAPSE == "Yes" ~ "Relapse",
                                                       response %in% c("PCR", "MPR") ~ "Responder",
                                                       .default = "SD"))

#Pre/Post Impact on CD8TCGS and BCGS (B/C) ####
#Identify patients with two measurements (Pre/Post)
no_follow_up <- (results %>% dplyr::count(ptID) %>% filter(n == 1))$ptID
pre.post_analysis <- results %>% filter(!ptID %in% no_follow_up)
pre.post_analysis$VISIT <- factor(pre.post_analysis$VISIT, levels = c("PRE", "POST"))

source(here('code/function_boxplot.R'))

g <- ggplot_samplelvl_boxplot(data = pre.post_analysis, x_variable = "VISIT", y_variable = "t_cell", default.theme = "default",
                         x.lab = "Treatment Visit", y.lab = "CD8TCGS", group.colors = NULL, comparisons = NULL)

ggsave(filename = "CD8TGCS_ABACUS_Pre_Post.tiff", path = output.path, plot = g, device = "tiff", width = 2, height = 3, units = "in", dpi = 320)

g <- ggplot_samplelvl_boxplot(data = pre.post_analysis, x_variable = "VISIT", y_variable = "b_cell", default.theme = "default",
                              x.lab = "Treatment Visit", y.lab = "BCGS", group.colors = NULL, comparisons = NULL)

ggsave(filename = "BCGS_ABACUS_Pre_Post.tiff", path = output.path, plot = g, device = "tiff", width = 2, height = 3, units = "in", dpi = 320)

#Response by B8T group (D/E/F) ####
source(here('code/function_stackedbar.R'))

#Subset for initial biopsy 
plot <- results %>% filter(VISIT == "PRE") 
plot$response <- factor(plot$response, levels = c("NR", "MPR", "PCR"))
cols <- c("PCR" = "darkgreen", "MPR" = "orange", "NR" = "red")

g <- stacked_bar(df = plot, x_variable = "B8T", group_variable = "response", 
                 annotate_n = FALSE, group.colors = cols,default.theme = "default", 
                 x = "B8T", y = "Frequency", fill = "Response")
                   
ggsave(filename = "B8T_ABACUS_Pre_Response.tiff", path = output.path, plot = g, device = "tiff", width = 3, height = 2, units = "in", dpi = 320, bg = "white")

cols = c("No" = "darkgreen", "Yes" = "orange")

g <- stacked_bar(df = plot, x_variable = "B8T", group_variable = "RELAPSE", annotate_n = FALSE, group.colors = cols,
                 x = "B8T", y = "Frequency", fill = "Relapse") 

ggsave(filename = "B8T_ABACUS_Pre_Relapse.tiff", path = output.path, plot = g, device = "tiff", width = 3, height = 2, units = "in", dpi = 320, bg = "white")

plot  <- results %>% filter(VISIT == "POST") 
g <- stacked_bar(df = plot, x_variable = "B8T", group_variable = "RELAPSE", annotate_n = FALSE, group.colors = cols,
                 x = "B8T", y = "Frequency", fill = "Relapse") 
ggsave(filename = "B8T_ABACUS_Post_Relapse.tiff", path = output.path, plot = g, device = "tiff", width = 3, height = 2, units = "in", dpi = 320, bg = "white")

# Sankey Diagram:
library(networkD3)

plot <- results %>% select(ptID, VISIT,B8T, response)

no_follow_up <- (plot %>% count(ptID) %>% filter(n == 1))$ptID

plot <- plot %>% filter(!ptID %in% no_follow_up)


nodes <- data.frame(name = c("Pre Hi/Hi", "Pre Hi/Lo", "Pre Lo/Hi", "Pre Lo/Lo", "Post Hi/Hi", "Post Hi/Lo", "Post Lo/Hi", "Post Lo/Lo"))

reshaped_data <- plot %>%
  pivot_wider(names_from = VISIT, values_from = B8T) %>%
  rename(PRE = `PRE`, POST = `POST`)
reshaped_data <- reshaped_data %>% mutate(PRE = str_c("Pre ", PRE),
                                          POST = str_c("Post ", POST))
links <- reshaped_data %>%
  group_by(PRE, POST) %>%
  summarise(Count = n(), .groups = 'drop')

nodes <- data.frame(name = unique(c(reshaped_data$PRE, reshaped_data$POST)))

links <- links %>%
  mutate(
    source = match(PRE, nodes$name) - 1,  # Get index of source node
    target = match(POST, nodes$name) - 1  # Get index of target node
  ) %>%
  select(source, target, Count)


node.cols = rep(c("Hi/Hi" = "#d1d43f",
                  "Hi/Lo" = "#2b628c",
                  "Lo/Hi" = "#3fd485",
                  "Lo/Lo" = "#523070"), 2)

nodes$color= node.cols
sankeyPlot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", Target = "target", Value = "Count", NodeID = "name", fontSize = 12, nodeWidth = 30)

print(sankeyPlot)
