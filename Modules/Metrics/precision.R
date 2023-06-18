# Precision Metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1
#' @param predict_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of precision.
#' 
#' 
#' 
#' 
rank_to_precision = function(df_rank, s_seq) {
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                precision = map(data,
                                function(d) summarise(d,
                                                      across(-PEDIGREE_NAME, ~ metric_precision(s, ., d[[as.character(full_harvest_repetition)]]))
                                )
                ), .groups="drop") %>%
      unnest(cols = precision)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, selection_intensity), names_to = "current_harvest_repetition", values_to = "precision") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}

metric_precision = function(selection_intensity, predict_rank, reference_rank){
  predict_selected = predict_rank <= selection_intensity*length(predict_rank)
  reference_selected = reference_rank <= selection_intensity*length(reference_rank)
  return(sum(predict_selected * reference_selected)/sum(predict_selected))
}

metric_precision_legacy = function(selection_intensity, df_model){
  df_rank = df_model %>%
    group_by(PEDIGREE_NAME, current_max_r) %>%
    summarise(mean_re = mean(fitted_value_re), .groups = "drop") %>%
    group_by(current_max_r) %>%
    mutate(rank = rank(mean_re)) %>% ungroup() %>%
    arrange(current_max_r) %>%
    select(PEDIGREE_NAME, current_max_r, rank)
  
  metric_precision_rank_df(selection_intensity,
                        df_rank %>% filter(current_max_r == 5),
                        df_rank %>% filter(current_max_r == 52))
  
}

metric_precision_rank_df = function(selection_intensity, predict_rank_df, reference_rank_df){
  predict_rank %>% mutate(selected = rank <= selection_intensity * length(rank)) %>%
    left_join(reference_rank %>% mutate(selected = rank <= selection_intensity * length(rank)), by = "PEDIGREE_NAME") %>%
    summarise(precision = sum(selected.x * selected.y)/sum(selected.x)) %>%
    pull(precision)
}





