# Jaccard Similarity Metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of precision.
#' 
#'
rank_to_jaccard_similarity = function(df_rank, s_seq) {
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                jaccard_similarity = map(data,
                                function(d) summarise(d,
                                                      across(-PEDIGREE_NAME, ~ metric_jaccard_similarity(s, ., d[[as.character(full_harvest_repetition)]]))
                                )
                ), .groups="drop") %>%
      unnest(cols = jaccard_similarity)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, selection_intensity), names_to = "current_harvest_repetition", values_to = "jaccard_similarity") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}

metric_jaccard_similarity = function(selection_intensity, predict_rank, reference_rank){
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
  mutate(predict_selected = predict_rank <= selection_intensity * length(predict_rank)) %>%
  mutate(reference_selected = reference_rank <= selection_intensity * length(reference_rank)) %>%
  summarise(js = sum(predict_selected & reference_selected)/max(sum(predict_selected | reference_selected), 1e-10)) %>%
    pull(js)
}

metric_jaccard_similarity_rank_df = function(selection_intensity, predict_rank, reference_rank){
  predict_rank %>% mutate(selected = rank <= selection_intensity * length(rank)) %>%
    left_join(reference_rank %>% mutate(selected = rank <= selection_intensity * length(rank)), by = "PEDIGREE_NAME") %>%
    summarise(js = sum(selected.x & selected.y)/max(sum(selected.x | selected.y), 1e-10)) %>%
    pull(js)
}