# Define Model 

fit_lme_poly = function(df, poly_max_degree, range_max_r) {
  max_repetition = df %>% summarise(max_repetition = max(as.integer(REPETITION))) %>% pull(max_repetition)
  poly_matrix = generate_poly(max_repetition, poly_max_degree)
  
  formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME) + ( poly0 | PLOT_BID:PEDIGREE_NAME)", collapse = "")
  formula_string_re = paste("~(0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME)", collapse = "")
  
  pb <- progress_bar$new(
    format = ":spin [:bar] :percent eta: :eta elapse: :elapsedfull    ",
    total = length(range_max_r), clear = FALSE, width= 120)
  
  m_list = list()
  df_list = list()
  suppressWarnings({
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
  })
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




parallel_fit_lme_poly = function(df, poly_max_degree, range_max_r) {
  max_repetition = df %>% summarise(max_repetition = max(as.numeric(REPETITION))) %>% pull(max_repetition)
  poly_matrix = generate_poly(max_repetition, poly_max_degree)
  
  formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME) + ( poly0 | PLOT_BID:PEDIGREE_NAME)", collapse = "")
  formula_string_re = paste("~(0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME)", collapse = "")
  cat("range_max_r = ", range_max_r, "\n")
  print("hello")
  output = foreach(max_r=range_max_r, .packages = c("tidyverse", "stringr", "lme4")) %dopar% {
    cat("Fitting up to repetitions", max_r, "\n")
    df_max_r = df %>%
      mutate(REPETITION = as.integer(REPETITION)) %>%
      filter(REPETITION <= max_r) %>%
      left_join(poly_matrix %>%
                  data.frame() %>%
                  rownames_to_column("REPETITION") %>%
                  mutate(REPETITION = as.integer(REPETITION)),
                by = "REPETITION")
    m_max_r = lmer(formula(formula_string),
             data = df_max_r)
    m_max_r@frame = m_max_r@frame %>% mutate(REPETITION = df_max_r$REPETITION)
    df_max_r = df_max_r %>%
      mutate(current_max_r = max_r,
             skip = "SKIP_NONE",
             fitted_value = predict(m_max_r),
             fitted_value_fe = predict(m_max_r, re.form = ~ 0),
             fitted_value_re = predict(m_max_r, re.form = formula(formula_string_re), random.only = T))
    list(df = df_max_r, model = m_max_r)
  }
  return(output)
}


parallel_fit_lme_poly2 = function(df, poly_max_degree, range_max_r) {
  max_repetition = df %>% summarise(max_repetition = max(as.numeric(REPETITION))) %>% pull(max_repetition)
  poly_matrix = generate_poly(max_repetition, poly_max_degree)
  
  formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME) + ( poly0 | PLOT_BID:PEDIGREE_NAME)", collapse = "")
  formula_string_re = paste("~(0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME)", collapse = "")
  cat("range_max_r = ", range_max_r, "\n")
  start_time <- Sys.time()
  output = mclapply(range_max_r, function(max_r) {
    cat("Fitting up to repetitions", max_r, "\n")
    df_max_r = df %>%
      mutate(REPETITION = as.integer(REPETITION)) %>%
      filter(REPETITION <= max_r) %>%
      left_join(poly_matrix %>%
                  data.frame() %>%
                  rownames_to_column("REPETITION") %>%
                  mutate(REPETITION = as.integer(REPETITION)),
                by = "REPETITION")
    m_max_r = lmer(formula(formula_string),
                   data = df_max_r)
    m_max_r@frame = m_max_r@frame %>% mutate(REPETITION = df_max_r$REPETITION)
    df_max_r = df_max_r %>%
      mutate(current_max_r = max_r,
             skip = "SKIP_NONE",
             fitted_value = predict(m_max_r),
             fitted_value_fe = predict(m_max_r, re.form = ~ 0),
             fitted_value_re = predict(m_max_r, re.form = formula(formula_string_re), random.only = T))
    list(df = df_max_r, model = m_max_r)
  }, mc.cores = 8)
  end_time <- Sys.time()
  cat("time elapsed = ", end_time - start_time, "mins.\n")
  return(output)
}


fit_lme_single = function(df, max_r, poly_matrix, formula_string, formula_string_re) {
  cat("Fitting up to repetitions", max_r, "\n")
  df = df %>%
    mutate(REPETITION = as.integer(REPETITION)) %>%
    filter(REPETITION <= max_r) %>%
    left_join(poly_matrix %>%
                data.frame() %>%
                rownames_to_column("REPETITION") %>%
                mutate(REPETITION = as.integer(REPETITION)),
              by = "REPETITION")
  m = lmer(formula(formula_string),
           data = df)
  m@frame = m@frame %>% mutate(REPETITION = df$REPETITION)
  df = df %>%
    mutate(current_max_r = max_r,
           skip = "SKIP_NONE",
           fitted_value = predict(m),
           fitted_value_fe = predict(m, re.form = ~ 0),
           fitted_value_re = predict(m, re.form = formula(formula_string_re), random.only = T))
  return(list(df = df, model = m))
}


