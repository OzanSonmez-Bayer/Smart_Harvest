# Genetic Gain metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Genetic Gain metric.

rank_to_genetic_gain = function(df_rank, s_seq, benchmark_harvest_repetition = NULL) {
  if(is.null(benchmark_harvest_repetition)) {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = full_harvest_repetition)
  } else {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = benchmark_harvest_repetition)
  }
  
  Reduce(rbind, lapply(s_seq, function(selection_intensity) {
    df_rank %>% left_join(df_rank %>%
                            filter(current_harvest_repetition == benchmark_harvest_repetition) %>%
                            rename(benchmark_harvest_rank = rank, benchmark_harvest_mean_value = mean_value) %>%
                            select(PEDIGREE_NAME, full_harvest_repetition, skip, benchmark_harvest_repetition, benchmark_harvest_rank, benchmark_harvest_mean_value),
                          by = c("PEDIGREE_NAME", "full_harvest_repetition", "skip", "benchmark_harvest_repetition")) %>%
      group_by(full_harvest_repetition, skip, benchmark_harvest_repetition, current_harvest_repetition) %>%
      mutate(selected = (rank <= selection_intensity*length(rank)),
             benchmark_harvest_selected = (benchmark_harvest_rank <= selection_intensity*length(benchmark_harvest_rank)) ) %>%
      summarise(selection_intensity = selection_intensity,
                genetic_gain = sum(benchmark_harvest_mean_value * selected)/sum(selected) - sum(benchmark_harvest_mean_value * benchmark_harvest_selected)/sum(benchmark_harvest_selected),
                .groups = "drop")
  })) %>%
    return()
}


