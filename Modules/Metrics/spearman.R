# Spearman Rank Correlation Metric
#
#' This function calculates the Spearman Rank Correlation metric, a non-parametric measure of rank correlation.
#' It assesses how well the relationship between two variables can be described using a monotonic function.
#' The function takes a data frame of ranks and a vector of selection intensities and returns the Spearman Rank Correlation.
#'
#' @param df_rank A data frame containing the ranks.
#' @param s_seq A vector of selection intensities, defaulting to c(1.0). Represents the percentage of total pedigrees one intends to select. Should be between 0 and 1.
#' @param ... Additional arguments to pass to the rank_to_metric function.
#'
#' @return A data frame containing the Spearman Rank Correlation metric.
#' 
rank_to_spearman = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "spearman", metric_spearman, ...)
}

# Spearman Rank Correlation Metric Function
#
#' This function computes the Spearman Rank Correlation between predicted and reference ranks.
#' It assesses how well the relationship between two variables can be described using a monotonic function.
#'
#' @param selection_intensity A numeric value representing the selection intensity.
#' @param predict_rank A numeric vector representing the predicted ranks.
#' @param reference_rank A numeric vector representing the reference ranks.
#' @return A numeric value representing the Spearman Rank Correlation.
#' 
metric_spearman = function(selection_intensity, predict_rank, reference_rank){
  if(sd(predict_rank) < 1e-12 || sd(reference_rank) < 1e-12) return(NA)
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    filter(predict_rank <= selection_intensity*length(predict_rank) | reference_rank <= selection_intensity*length(reference_rank)) %>%
    summarise(spearman = cor(predict_rank, reference_rank, method = "spearman") ) %>%
    pull(spearman)
}
