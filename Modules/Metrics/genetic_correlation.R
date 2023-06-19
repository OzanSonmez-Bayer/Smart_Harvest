# Genetic Correlation metric
#
#' @param df_rank
#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' 
#' @return the numeric of Genetic Correlation.
fitted_to_genetic_correlation = function(df, s_seq) {
  df_rank = df %>%
    group_by(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip) %>%
    summarise(mean_value = mean(fitted_value_re), .groups = "drop") %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    mutate(rank = rank(mean_value)) %>% ungroup() %>%
    arrange(current_harvest_repetition) %>%
    select(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip, mean_value, rank)
  
  Reduce(rbind, lapply(s_seq, function(selection_intensity) {
    df_rank %>% left_join(df_rank %>%
                            filter(current_harvest_repetition == full_harvest_repetition) %>%
                            rename(full_harvest_rank = rank, full_harvest_mean_value = mean_value) %>%
                            select(PEDIGREE_NAME, full_harvest_repetition, skip, full_harvest_rank, full_harvest_mean_value),
                          by = c("PEDIGREE_NAME", "full_harvest_repetition", "skip")) %>%
      group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
      filter( (rank <= selection_intensity*length(rank)) | (full_harvest_rank <= selection_intensity*length(full_harvest_rank)) ) %>%
      summarise(selection_intensity = selection_intensity, cor = cor(mean_value, full_harvest_mean_value), .groups = "drop")
  })) %>%
    return()
}



