# Genetic Gain metric
#
#' @param df_rank
#' @param s_seq a vector of selection intensities, what % of total pedigrees one inteds to select. Each entry in the vector should be a number bewteen 0 and 1.
#' @param benchmark_harvest_repetition a integer indicating which harvest is the benchmark
#'
#' @return a data frame recording the numeric of Genetic Gain corresponding to the selection intensities.
#' 
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


