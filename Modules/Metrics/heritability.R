# Heritability metric
#
#' @param df the data frame of fitted mixed effect model with sigma_u (sd of random effect) and sigma_e (sd of noise) 
#'
#' @return a data frame recording the numeric of Heritability corresponding to different current harvest repetition
#' 
fitted_to_heritability = function(df) {
  df %>%
    select(full_harvest_repetition, current_harvest_repetition, skip, sigma_u, sigma_e) %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    summarise(h2=mean(sigma_u^2/(sigma_u^2+sigma_e^2)), .groups="drop") %>%
    return()
}
