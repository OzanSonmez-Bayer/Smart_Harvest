fitted_to_rank = function(df, fitted_col) {
  df %>%
    group_by(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip) %>%
    summarise(mean_value = mean({{fitted_col}}), .groups = "drop") %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    mutate(rank = rank(mean_value)) %>% ungroup() %>%
    arrange(current_harvest_repetition) %>%
    select(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip, rank) %>%
    return()
}

rank_random_effect = function(df, rank_sign = 1) {
  # If rank_sign = 1 then larger values has larger rank (e.g. (-5.4, 2.6, 7.2) => (1, 2, 3))
  # If rank_sign = -1 then smaller values has larger rank (e.g. (-5.4, 2.6, 7.2) => (3, 2, 1))
  rank_df = df %>%
    select(PEDIGREE_NAME, REPETITION, fitted_value_re) %>%
    arrange(PEDIGREE_NAME, REPETITION) %>%
    group_by(PEDIGREE_NAME) %>%
    summarise(mean_re = mean(fitted_value_re), .groups = "drop") %>%
    mutate(rank = rank(rank_sign*mean_re))
  return(rank_df)
}

rank_col = function(df, col, rank_sign = 1) {
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  rank_df = df %>%
    select(PEDIGREE_NAME, REPETITION, {{col}} ) %>%
    arrange(PEDIGREE_NAME, REPETITION) %>%
    group_by(PEDIGREE_NAME) %>%
    summarise(mean = mean({{col}}), .groups = "drop") %>%
    mutate(rank = rank(rank_sign*mean))
  return(rank_df)
}


expected_rank = function(df, rank_sign = 1) {
  # WIP
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  rank_df = df %>%
    select(PEDIGREE_NAME, REPETITION, fitted_value_re, condsd_re) %>%
    arrange(PEDIGREE_NAME, REPETITION) %>%
    group_by(PEDIGREE_NAME) %>%
    summarise(mean_re = mean(fitted_value_re), sd_re = sqrt(mean(condsd_re^2)), .groups = "drop") %>%
  mu = rank_df %>% pull(mean_re, PEDIGREE_NAME)
  sigma = rank_df %>% pull(sd_re, PEDIGREE_NAME)
  rank_df %>% left_join(
    get_expectedRank(mu, sigma, rank_sign),
    by = "PEDIGREE_NAME"
  ) %>% return()
}

get_expectedRank = function(names_vector, mu, sigma, rank_sign = 1) {
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  
  expected_rank = apply(p_matrix(mu, sigma), (-rank_sign + 3) %/%2, sum)
  
  rank_on_er = rank(expected_rank)
  sd_er = sqrt(v_fun(mu, sigma))
  
  percentile_er <- round(100 * (expected_rank - 0.5)/dim(ranef_df)[1])
  return(data.frame(PEDIGREE_NAME = names_vector,
                    expected_rank = expected_rank,
                    rank_on_er = rank_on_er,
                    percentile_er = percentile_er,
                    sd_er = sd_er))
}
