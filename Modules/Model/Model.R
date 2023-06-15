# Define Model 

fit_lme_poly = function(df, poly_max_degree) {
  max_repetition = df %>% summarise(max_repetition = max(as.numeric(REPETITION))) %>% pull(max_repetition)
  poly_matrix = generate_poly(max_repetition, poly_max_degree)
  df = df %>%
    left_join(poly_matrix %>%
                data.frame() %>%
                rownames_to_column("REPETITION") %>%
                mutate(REPETITION = as.integer(REPETITION)),
              by = "REPETITION")
  
  formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME) + ( poly0 | PLOT_BID:PEDIGREE_NAME)", collapse = "")
  formula_string_re = paste("~(0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME)", collapse = "")
  
  m = lmer(formula(formula_string),
           data = df)
  m@frame = m@frame %>% mutate(REPETITION = df$REPETITION)
  
  df = df %>% mutate(fitted_value = predict(m),
                     fitted_value_fe = predict(m, re.form = ~ 0),
                     fitted_value_re = predict(m, re.form = formula(formula_string_re), random.only = T)
  )
  
  return(list(df = df, model = m))
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


