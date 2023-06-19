# Genetic Gain metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Genetic Gain metric.

fitted_to_genetic_gain = function(df, s_seq) {
  Reduce(rbind, lapply(s_seq, function(selection_intensity) {
    df %>%
      fitted_to_rank(fitted_value_re) %>%
      group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
      filter(rank <= selection_intensity*length(rank)) %>%
      left_join(df %>%
                  filter(current_harvest_repetition == full_harvest_repetition) %>%
                  select(PEDIGREE_NAME, full_harvest_repetition, skip, fitted_value_re) %>%
                  group_by(PEDIGREE_NAME, full_harvest_repetition, skip) %>%
                  summarise(mean_re_full_harvest = mean(fitted_value_re), .groups = "drop") %>%
                  group_by(full_harvest_repetition, skip) %>%
                  mutate(single_mean_re_full_harvest=mean(mean_re_full_harvest)),
                by = c("PEDIGREE_NAME", "full_harvest_repetition", "skip")) %>%
      group_by(full_harvest_repetition, current_harvest_repetition) %>%
      summarise(
        selection_intensity = selection_intensity,
        genetic_gain = mean(mean_re_full_harvest - single_mean_re_full_harvest), .groups = "drop") %>%
      mutate(normalized_genetic_gain = genetic_gain/max(abs(genetic_gain)))
  }))
}


