# Define Model 

fit_lme_poly = function(df, poly_max_degree, range_max_r) {
  max_repetition = df %>% summarise(max_repetition = max(as.integer(REPETITION))) %>% pull(max_repetition)
  poly_matrix = generate_poly(max_repetition, poly_max_degree)
  
  #formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "|| PEDIGREE_NAME) + ( 0 + poly0 ||
  #                       PLOT_BID)", collapse = "")
  formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "|| PEDIGREE_NAME)", collapse = "")
  
  formula_string_re = paste("~(0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "|| PEDIGREE_NAME)", collapse = "")
  
  pb <- progress_bar$new(
    format = ":spin [:bar] :percent eta: :eta elapse: :elapsedfull    ",
    total = length(range_max_r), clear = FALSE, width= 80)
  
  m_list = list()
  df_list = list()
  for(max_r in range_max_r){
    pb$tick()
    df_list[[max_r]] = df %>%
      mutate(REPETITION = as.integer(REPETITION)) %>%
      filter(REPETITION <= max_r) %>%
      left_join(poly_matrix %>%
                  data.frame() %>%
                  rownames_to_column("REPETITION") %>%
                  mutate(REPETITION = as.integer(REPETITION)),
                by = "REPETITION")
      m_list[[max_r]] = lmer(formula(formula_string),
               data = df_list[[max_r]])
      m_list[[max_r]]@frame = m_list[[max_r]]@frame %>% mutate(REPETITION = df_list[[max_r]]$REPETITION)
      df_list[[max_r]] = df_list[[max_r]] %>%
        mutate(full_harvest_repetition = as.integer(max_repetition),
               current_harvest_repetition = as.integer(max_r),
               skip = "SKIP_NONE",
               fitted_value = predict(m_list[[max_r]]),
               fitted_value_fe = predict(m_list[[max_r]], re.form = ~ 0),
               fitted_value_re = predict(m_list[[max_r]], re.form = formula(formula_string_re), random.only = T))
  }
  cat("\nDone!\n")
  df_output = Reduce(rbind, lapply(range_max_r, function(k) df_list[[k]]))
  return(list(df = df_output, model = m_list))
}

# Generate orthognal polynomial 
generate_poly = function(max_repetition, poly_max_degree) {
  if(poly_max_degree > 0) {
    df_poly = data.frame(cbind(1/sqrt(max_repetition),
                               poly(1:max_repetition, degree=poly_max_degree, simple = TRUE)))
    colnames(df_poly) = paste0("poly", 0:poly_max_degree)
  } else {
    df_poly = data.frame(1/sqrt(max_repetition))
    colnames(df_poly) = c("poly0")
  }
  rownames(df_poly) = 1:max_repetition
  df_poly %>% as.matrix() %>% return()
}

