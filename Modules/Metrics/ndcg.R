# nDCG Metric
#
#' This function calculates the normalized Discounted Cumulative Gain (nDCG) metric, a measure used to evaluate ranking quality.
#' It takes into account the position of the rank and is particularly useful in recommendation systems.
#' The function takes a data frame of ranks and a vector of selection intensities and returns the nDCG.
#'
#' @param df_rank A data frame containing the ranks.
#' @param s_seq A vector of selection intensities, defaulting to c(1.0). Represents the percentage of total pedigrees one intends to select. Should be between 0 and 1.
#' @param ... Additional arguments to pass to the rank_to_metric function.
#'
#' @return A data frame containing the nDCG metric corresponding to selection intensities.
#'
rank_to_ndcg = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "ndcg", metric_ndcg, ...)
}

# nDCG Metric Function
#
#' This function computes the normalized Discounted Cumulative Gain (nDCG) between predicted and reference ranks.
#' It considers the position of the rank and is useful in evaluating ranking quality.
#'
#' @param selection_intensity A numeric value representing the selection intensity.
#' @param predict_rank A numeric vector representing the predicted ranks.
#' @param reference_rank A numeric vector representing the reference ranks.
#' @return A numeric value representing the nDCG.
#' 
metric_ndcg = function(selection_intensity, predict_rank, reference_rank){
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    mutate(scores = reference_rank <= selection_intensity * length(reference_rank)) %>%
    filter(predict_rank <= selection_intensity * length(predict_rank)) %>%
    arrange(predict_rank) %>%
    summarise(nDCG = sum(scores/log2(predict_rank + 1))/sum(1/log2(predict_rank + 1)) ) %>%
    pull(nDCG)
}

# Function to calculate DCG
calc_DCG <- function(scores) {
  DCG = sum(scores / log2(1:length(scores) + 1))
  return(DCG)
}
