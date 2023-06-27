# Genetic Correlation metric
#
#' @param df_rank
#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' 
#' @return the numeric of Genetic Correlation.
fitted_to_heritability = function(df) {
  df %>%
    select(full_harvest_repetition, current_harvest_repetition, skip, sigma_u, sigma_e) %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    summarise(h2=mean(sigma_u^2/(sigma_u^2+sigma_e^2)), .groups="drop") %>%
    return()
}


