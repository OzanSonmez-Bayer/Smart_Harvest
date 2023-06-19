# Kendall Rank Correlation metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Kendall Rank Correlation.
rank_to_kendall = function(df_rank, s_seq) {
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                kendall = map(data,
                                function(d) summarise(d,
                                                      across(-PEDIGREE_NAME, ~ metric_kendall(s, ., d[[as.character(full_harvest_repetition)]]))
                                )
                ), .groups="drop") %>%
      unnest(cols = kendall)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, selection_intensity), names_to = "current_harvest_repetition", values_to = "kendall") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}

metric_kendall = function(selection_intensity, predict_rank, reference_rank){
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    filter(predict_rank <= selection_intensity*length(predict_rank) | reference_rank <= selection_intensity*length(reference_rank)) %>%
    summarise(kendall = cor(predict_rank, reference_rank, method = "kendall")) %>%
    pull(kendall)
}

metric_kendall_rank_df = function(selection_intensity, predict_rank_df, reference_rank_df){
  predict_rank_df %>%
    left_join(reference_rank_df, by = "PEDIGREE_NAME") %>%
    filter(rank.x <= selection_intensity*length(rank.x) | rank.y <= selection_intensity*length(rank.y)) %>%
    summarise(kendall = cor(rank.x, rank.y, method = "kendall") ) %>%
    pull(kendall)
}
