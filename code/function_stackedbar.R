#' @param data dataframe object with x and y variable to plot
#' @param x_variable x variable for box plot. Also becomes the fill
#' @param group_variable y variable for box plot.
#' @param group.colors named vector of colors for groups
#' @param annotate_n Default 'FALSE'. If true, will annotate the count (n) within each element of the stacked bar plot
#' @param default.theme TRUE. Uses default custom theme. Can accept a ggplot object to replace
#' @param ... Arguments past on to "labs" 

df = plot
x_variable = "t_cell"
group_variable = "response"
stacked_bar <- function(df, x_variable, group_variable, stat = c("fisher.test", "chisq.test")[1],
                        group.colors = NULL, annotate_n = FALSE, 
                        default.theme = "default", ...) {
  
  format_pvalue <- function(pvalue) {
    # If the p-value is exactly 1, return "1"
    if(pvalue == 1) {
      return("1")
    }
    # For very low p-values, report as "<0.001"
    if(pvalue < 0.001) {
      return("<0.001")
    }
    # If the p-value is in the range considered "close to 0.05" (0.01 to 0.2 inclusive),
    # report it with 2 significant figures.
    if(pvalue >= 0.01 && pvalue <= 0.2) {
      return(as.character(signif(pvalue, 2)))
    }
    # Otherwise, report with 1 significant figure.
    return(as.character(signif(pvalue, 1)))
  }
  
  #Calculate pvalue: 
  if(stat == "fisher.test"){
    p_value <- fisher.test(table(df[[group_variable]], df[[x_variable]]))
  }
  if(stat == "chisq.test"){
    p_value <- chisq.test(table(df[[group_variable]], df[[x_variable]]))
  }
  p_value = p_value$p.value
  formatted_p_value <- format_pvalue(p_value)
  
  
  ## Set Plot theme
  if(is.character(default.theme) && default.theme == "default"){
    default.theme <- theme_minimal() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"),
            text = element_text(size = 10, colour = "black"))
    
  } else{
    default.theme = default.theme
  }
  
  
  ## Generate Plot
  g <- ggplot(data = df, aes(x = .data[[x_variable]])) + geom_bar(position = "fill", aes(fill = .data[[group_variable]])) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), breaks = seq(from = 0, to = 1, by = 0.1)) +
    annotate("text",x = 2, y = 1.05, label = paste0("p = ", formatted_p_value), size = 3, color = "black", fontface = "bold")
  if(!is.null(group.colors)){
    g <- g + scale_fill_manual(values = group.colors)
    } 
  if(annotate_n){
    #Calculate frequency annotation and position
    df_counts <- df %>%
      group_by(.data[[x_variable]], .data[[group_variable]]) %>%
      summarise(n = n()) %>%
      mutate(freq = n / sum(n),  # the frequency within each subgroup
             cumfreq = cumsum(freq)) %>%  # cumulative sum to stack the bars
      ungroup()
    # Position
    df_counts$text_pos <- with(df_counts, 1- cumfreq + (0.5 * freq))
   g <- g + geom_text(data = df_counts, aes(x = .data[[x_variable]], label = n, y = text_pos), vjust = 0.5, size = 4)
  }
  
  g <- g + labs(...) + default.theme 
  
  return(g)
  
  
}

