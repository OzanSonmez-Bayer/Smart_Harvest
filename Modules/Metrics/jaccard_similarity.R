# Jaccard Similarity Metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of precision.
#' 
#'
rank_to_jaccard_similarity = function(df_rank, s_seq, ...) {
  rank_to_metric(df_rank, s_seq, "jaccard_similarity", metric_jaccard_similarity, ...)
}

metric_jaccard_similarity = function(selection_intensity, predict_rank, reference_rank){
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
  mutate(predict_selected = predict_rank <= selection_intensity * length(predict_rank)) %>%
  mutate(reference_selected = reference_rank <= selection_intensity * length(reference_rank)) %>%
  summarise(js = sum(predict_selected & reference_selected, na.rm = T)/max(sum(predict_selected | reference_selected, na.rm = T), 1e-10)) %>%
    pull(js)
}
