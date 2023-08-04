#' Jaccard Similarity Metric
#'
#' @param df_rank Data frame containing the rank information.
#' @param s_seq A vector of selection intensities. Represents the percentage of total pedigrees one intends to select. Should be between 0 and 1.
#' @param benchmark_harvest_repetition An integer indicating which harvest is the benchmark (optional).
#' @return Jaccard Similarity
#' 
rank_to_jaccard_similarity = function(df_rank, s_seq, ...) {
  rank_to_metric(df_rank, s_seq, "jaccard_similarity", metric_jaccard_similarity, ...)
}

#' Jaccard Similarity Metric Function
#'
#' This function computes the Jaccard Similarity between predicted and reference ranks.
#' It calculates the intersection over union of selected items based on the selection intensity.
#'
#' @param selection_intensity A numeric value representing the selection intensity.
#' @param predict_rank A numeric vector representing the predicted ranks.
#' @param reference_rank A numeric vector representing the reference ranks.
#' @return A numeric value representing the Jaccard Similarity.
#' 
metric_jaccard_similarity = function(selection_intensity, predict_rank, reference_rank){
  # Implementation of Jaccard Similarity Metric
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
  mutate(predict_selected = predict_rank <= selection_intensity * length(predict_rank)) %>%
  mutate(reference_selected = reference_rank <= selection_intensity * length(reference_rank)) %>%
  summarise(js = sum(predict_selected & reference_selected, na.rm = T)/max(sum(predict_selected | reference_selected, na.rm = T), 1e-10)) %>%
    pull(js)
}
