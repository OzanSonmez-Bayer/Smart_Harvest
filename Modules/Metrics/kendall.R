# Kendall Rank Correlation Metric
#
#' This function calculates the Kendall Rank Correlation metric, a measure of ordinal association between two measured quantities.
#' It assesses the strength and direction of the relationship between two ranks.
#' The function takes a data frame of ranks and a vector of selection intensities and returns the Kendall Rank Correlation.
#'
#' @param df_rank A data frame containing the ranks.
#' @param s_seq A vector of selection intensities, defaulting to c(1.0). Represents the percentage of total pedigrees one intends to select. Should be between 0 and 1.
#' @param ... Additional arguments to pass to the rank_to_metric function.
#'
#' @return A data frame containing the Kendall Rank Correlation metric.
#' 
rank_to_kendall = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "kendall", metric_kendall, ...)
}

metric_kendall = function(selection_intensity, predict_rank, reference_rank){
  if(sd(predict_rank) < 1e-12 || sd(reference_rank) < 1e-12) return(NA)
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    filter(predict_rank <= selection_intensity*length(predict_rank) | reference_rank <= selection_intensity*length(reference_rank)) %>%
    summarise(kendall = cor(predict_rank, reference_rank, method = "kendall")) %>%
    pull(kendall)
}
