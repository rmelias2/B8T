#' Function to generate box plot with statistic
#' 
#' @param data dataframe object with x and y variable to plot
#' @param x_variable x variable for box plot. Also becomes the fill
#' @param y_variable y variable for box plot.
#' @param subset_variable Variable on which to subset values in the dataframe
#' @param subset_values vector of values within the subset variable on which to subset the data
#' @param default.theme TRUE. Uses default custom theme. Can accept a ggplot object to replace
#' @param x.lab NULL - If specified will alter the x label 
#' @param y.lab NULL - If specified will alter the y label 
#' @param group.colors named vector of colors for groups


ggplot_samplelvl_boxplot <- function(data, x_variable, y_variable, interaction.term = FALSE,
                                     subset_variable = NULL, subset_values = NULL, default.theme = "default", x.lab = NULL,
                                     y.lab = NULL, group.colors = NULL, comparisons = NULL, y.limits = NULL, pval.position = "centre"){
  
  required_packages <- c("ggpubr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is required but not installed. Please install it using install.packages()."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  
  ## Subset Plot Data
  if(!is.null(subset_variable)){
    data <- data %>% filter(.data[[subset_variable]] %in% subset_values)
  }
  plot_data <- data
  
  ## Set Plot theme
  if(is.character(default.theme) && default.theme == "default"){
    default.theme <- theme_minimal() + 
      theme(
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_blank(),
        text = element_text(family = "sans", color = "black"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8, face = "bold"),
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        legend.position = "none") 
  } else{
    default.theme = default.theme
  }
  
  symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", ""))
  ## Generate Plot
  if(interaction.term != FALSE){
    g <- ggplot(plot_data, aes(x = .data[[x_variable]] , y = .data[[y_variable]], fill = .data[[interaction.term]], 
                               group = interaction(.data[[x_variable]], .data[[interaction.term]]))) + 
      geom_boxplot(outlier.alpha = 0) 
  } else {
    g <- ggplot(plot_data, aes(x = .data[[x_variable]] , y = .data[[y_variable]], fill = .data[[x_variable]])) +
      geom_boxplot(outlier.alpha = 0) + 
      geom_jitter(size = 0.5, width = 0.2, height = 0)
  }
  if(!is.null(group.colors)){
    g <- g + scale_fill_manual(values = group.colors)
  } else{
    g <- g+scale_fill_brewer(palette = "Set1")
  }
  if(!is.null(y.limits)){
    g <- g + ylim(y.limits[1], y.limits[2])
  }
  if(!is.null(comparisons)){
    g <- g + stat_compare_means(symnum.args = symnum.args, label.x.npc = pval.position, size = 3, comparisons = comparisons)
  } else{
    g <- g + stat_compare_means(label = "p.format", label.x.npc = pval.position, size = 3)
  }
  if(!is.null(x.lab)){
    g <- g + xlab(x.lab)
  }
  if(!is.null(y.lab)){
    g <- g +  ylab(y.lab)
  }
  g <- g + default.theme 
  
  return(g)
}
