# nDCG metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of nDCG.
#' 
#' 
rank_to_ndcg = function(df_rank, s_seq) {
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                ndcg = map(data,
                                function(d) summarise(d,
                                                      across(-PEDIGREE_NAME, ~ metric_ndcg(s, ., d[[as.character(full_harvest_repetition)]]))
                                )
                ), .groups="drop") %>%
      unnest(cols = ndcg)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, selection_intensity), names_to = "current_harvest_repetition", values_to = "ndcg") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}

metric_ndcg = function(selection_intensity, predict_rank, reference_rank){
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    mutate(scores = reference_rank <= selection_intensity * length(reference_rank)) %>%
    filter(predict_rank <= selection_intensity * length(predict_rank)) %>%
    arrange(predict_rank) %>%
    summarise(nDCG = sum(scores/log2(predict_rank + 1))/sum(1/log2(predict_rank + 1)) ) %>%
    pull(nDCG)
}


metric_ndcg_rank_df = function(selection_intensity, predict_rank_df, reference_rank_df){
  predict_rank_df %>% filter(rank <= selection_intensity * length(rank)) %>%
    arrange(rank) %>%
    left_join(reference_rank_df %>% mutate(scores = rank <= selection_intensity * length(rank)) %>% select(PEDIGREE_NAME, scores), by = "PEDIGREE_NAME") %>%
    summarise(nDCG = sum(scores/log2(rank + 1))/sum(1/log2(rank + 1)) ) %>%
    pull(nDCG)
}

# Function to calculate DCG
calc_DCG <- function(scores) {
  DCG = sum(scores / log2(1:length(scores) + 1))
  return(DCG)
}
