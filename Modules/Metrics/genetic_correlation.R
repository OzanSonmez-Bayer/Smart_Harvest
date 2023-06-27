# Genetic Correlation metric
#
#' @param df_rank
#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' 
#' @return the numeric of Genetic Correlation.
rank_to_genetic_correlation = function(df_rank, s_seq = c(1.0), benchmark_harvest_repetition = NULL) {
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
      filter( (rank <= selection_intensity*length(rank)) | (benchmark_harvest_rank <= selection_intensity*length(benchmark_harvest_rank)) ) %>%
      summarise(selection_intensity = selection_intensity,
                genetic_correlation = cor(mean_value, benchmark_harvest_mean_value), .groups = "drop")
  })) %>%
    return()
}

