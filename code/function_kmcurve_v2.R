#' @param Comparison Groups for curves
#' @param Data Dataframe which contains survival_time and survival_event variables
#' @param survival_time Name of column which contains survival time
#' @param survival_event Name of column which contains survival event
#' @param subset_label Character describing how the datawas subset, goes above number at risk 
#' @param xaxis.label Label for x.axis 
#' @param plot_colors names vector of colors corresponding to groups in comparison 

km_curve <- function(
    Data,
    groups,
    group_label = NULL,
    survival_time  = "OS_months",
    survival_event = "OS_numeric",
    subset_label,
    xaxis.label    = "Months from randomization", 
    plot_colors,
    width = 3, height = 3.5,
    ...
) {
  
  # Load required packages
  required_packages <- c("ggsurvfit","tidyverse","here")
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) message("Please install: ", pkg)
    else library(pkg, character.only = TRUE)
  }))
  
  surv.object <- survfit2(
    Surv(Data[[survival_time]], Data[[survival_event]]) ~ Data[[groups]],
    data = Data
  )
  
  # Extract tidy survival data
  tidy_df <- ggsurvfit::tidy_survfit(surv.object)
  
  # compute the last time before any group drops below 5 at risk [Per Journal Guidelines]
  t_cut <- tidy_df %>%
    group_by(strata) %>%
    summarize(last_ok = max(time[n.risk >= 5], na.rm = TRUE), .groups = "drop") %>%
    summarize(cutoff_time = min(last_ok, na.rm = TRUE)) %>%
    pull(cutoff_time)
  
  # if nothing meets the threshold, default to max observed time
  if (is.infinite(t_cut) || is.na(t_cut)) {
    t_cut <- max(tidy_df$time, na.rm = TRUE)
  } else {
    max(tidy_df$time, na.rm = TRUE)
  }
  
  # ---- (C) build the figure ----
  g <- surv.object %>%
    ggsurvfit(linewidth = 1) +
    add_risktable(
      risktable_stats = c("n.risk"),
      stats_label     = paste0(subset_label, ", No. at risk"),
      size            = 2.5,
      risktable_height = 0.20,
      theme = list(
        theme_risktable_default(axis.text.y.size = 8, plot.title.size = 8),
        theme(plot.title = element_text(face = "bold"))
      )
    ) +
    
    coord_cartesian(xlim = c(0, t_cut)) +    
    xlab(xaxis.label) +                       
    ylab("Overall survival") +              
    labs(color = group_label) +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent,
      expand = c(0.01, 0)
    ) +
    scale_color_manual(values = plot_colors) +
    theme(
      legend.position = "top",
      text = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # ---- (D) save ----
  outputname <- paste0(stringr::str_replace_all(subset_label, pattern = ",| |/", replacement = "_"), group_label)
  ggsave(filename = paste0(outputname, ".tiff"),
         path = output.path, plot = g, device = "tiff",
         width = width, height = height, units = "in", dpi = 320)
  
  
  invisible(g)
}

km_pvalue <- function(Data, groups, survival_time = "OS_months", survival_event = "OS_numeric") {
  stopifnot(all(c(survival_time, survival_event) %in% names(Data)))
  f <- survival::survdiff(stats::as.formula(
    sprintf("Surv(%s, %s) ~ `%s`", survival_time, survival_event, groups)
  ), data = Data)
  p <- 1 - pchisq(f$chisq, length(f$n) - 1)
  p
}


#' @title Landmark Survival
#' @description Returns the survival probability at a fixed time (e.g., 24 months)
#' @param Data Data frame containing survival time, event, and group columns
#' @param groups Variable defining comparison groups
#' @param survival_time Name of survival time variable (in months)
#' @param survival_event Name of event indicator (1=event, 0=censored)
#' @param landmark Landmark time in months (default = 24 for 2-year survival)
#' @return Tibble with group and estimated survival probability at landmark
#' @examples
#' km_landmark(df, groups = "b_cell_bin")
km_landmark <- function(Data,
                        groups,
                        survival_time  = "OS_months",
                        survival_event = "OS_numeric",
                        landmark       = 24) {
  
  # Fit KM
  fit <- survfit2(
    stats::as.formula(sprintf("Surv(%s, %s) ~ `%s`",
                              survival_time, survival_event, groups)),
    data = Data
  )
  
  # Extract tidy survival data
  surv_df <- ggsurvfit::tidy_survfit(fit) %>%
    group_by(strata) %>%
    summarize(
      surv_landmark = {
        idx <- which.min(abs(time - landmark))
        estimate[idx]  # <- FIX: use 'estimate', not 'surv'
      },
      .groups = "drop"
    ) %>%
    mutate(landmark_months = landmark)
  
  return(surv_df)
}
