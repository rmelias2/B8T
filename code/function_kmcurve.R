#' @param Comparison Groups for curves
#' @param Data Dataframe which contains survival_time and survival_event variables
#' @param survival_time Name of column which contains survival time
#' @param survival_event Name of column which contains survival event
#' @param subset_label Character describing how the datawas subset, goes above number at risk 
#' @param xaxis.label Label for x.axis 
#' @param plot_colors names vector of colors corresponding to groups in comparison 


km_curve <- function(Data, groups, group_label = NULL, survival_time = "OS_months", survival_event = "OS_numeric", subset_label,
                     xaxis.label = "Overall Survival (mo)", plot_colors, width = 3, height = 3.5, pval.x = 5, ...) {
  
  format_pvalue <- function(pvalue) {
  
    if (pvalue < 0.001) {
      return(format(pvalue, scientific = TRUE, digits = 3))
    } else {
      # Ensure there are always 3 digits after the decimal point
      return(sprintf("%.3f", pvalue))
    }
  }
  
  required_packages <- c("ggsurvfit","tidyverse","here")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Package", pkg, "is required but not installed. Please install it using install.packages()."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  
  surv.object <- survfit2(Surv(Data[[survival_time]], Data[[survival_event]]) ~ Data[[groups]], data = Data) 
  
  g <- surv.object %>% ggsurvfit(linewidth = 1) + 
    add_risktable(risktable_stats = c("n.risk"),  
                  stats_label = paste0(subset_label,  ", No. at risk"), 
                  size = 2.5, 
                  risktable_height = 0.2,
                  theme = list(theme_risktable_default(
                    axis.text.y.size = 8, 
                    plot.title.size = 8),
                    theme(plot.title = element_text(face = "bold")))) + 
    add_censor_mark(size = 2, alpha = 0.8) +
    add_pvalue("annotation", x = pval.x, y = 0.1, size = 3, pvalue_fun = format_pvalue) + 
    xlab(xaxis.label) +
    labs(color = group_label) + 
    theme_classic() + 
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent, 
      expand = c(0.01, 0)) +
    scale_color_manual(values = plot_colors) +
    theme(legend.position = "top") +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
  
  outputname = paste0(str_replace_all(subset_label, pattern = ",| |/", replacement = "_"), group_label)

  ggsave(filename = paste0(outputname,".tiff"), path = output.path, plot = g, device = "tiff", width = width, height = height, units = "in", dpi = 320)
  
}
