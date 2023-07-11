# Define Model 

fit_lme_poly = function(df,
                        poly_max_degree = 0,
                        range_max_r = NULL,
                        has_growseason = TRUE,
                        has_field_name = TRUE,
                        skip = "SKIP_NONE",
                        full_harvest_repetition = NULL,
                        verbose = TRUE) {
  if(is.null(full_harvest_repetition)){
    max_repetition = df %>% summarise(max_repetition = max(as.integer(REPETITION))) %>% pull(max_repetition)  
  } else {
    max_repetition = full_harvest_repetition
  }
  poly_matrix = generate_poly(max_repetition, poly_max_degree)
  if(is.null(range_max_r)) {range_max_r = 2:max_repetition}
  
  #formula_string = paste("TRAIT_VALUE ~ 0 + GROWSEASON + FIELD_NAME + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME) + ( 0 + poly0 ||
  #                       PLOT_BID)", collapse = "")

  growseason_string = ifelse(has_growseason, "+ GROWSEASON", "")
  field_name_string = ifelse(has_field_name, "+ FIELD_NAME", "")
  formula_string = paste("TRAIT_VALUE ~ 1 ",  growseason_string, field_name_string, " + (0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME)", collapse = "")  
  
  formula_string_re = paste("~(0 +", paste(paste0("poly", 0:poly_max_degree), collapse = " + "), "| PEDIGREE_NAME)", collapse = "")
  what_not_to_skip = switch(skip, SKIP_NONE = c(0,1), SKIP_EVEN = c(1), SKIP_ODD = c(2), c(0,1))

  pb <- progress_bar$new(
    format = ":spin [:bar] :percent eta: :eta elapse: :elapsedfull    ",
    total = length(range_max_r), clear = FALSE, width= 80)
  
  m_list = list()
  df_list = list()
  for(max_r in range_max_r){
    if(verbose) pb$tick()
    df_list[[max_r]] = df %>%
      mutate(REPETITION = as.integer(REPETITION)) %>%
      filter( (REPETITION <= max_r) & (REPETITION %% 2 %in% what_not_to_skip) ) %>%
      left_join(poly_matrix %>%
                  data.frame() %>%
                  rownames_to_column("REPETITION") %>%
                  mutate(REPETITION = as.integer(REPETITION)),
                by = "REPETITION")
      m_list[[max_r]] = lmer(formula(formula_string),
                             data = df_list[[max_r]], 
                             control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

      df_ranef = data.frame(ranef(m_list[[max_r]])) %>%
        filter(grpvar == "PEDIGREE_NAME", term == "poly0") %>%
        rename(PEDIGREE_NAME = grp, condval_u = condval, condsd_u = condsd) %>%
        select(PEDIGREE_NAME, condval_u, condsd_u)

      df_list[[max_r]] = df_list[[max_r]] %>%
        mutate(full_harvest_repetition = as.integer(max_repetition),
               current_harvest_repetition = as.integer(max_r),
               skip = skip,
               fitted_value = predict(m_list[[max_r]]),
               fitted_value_fe = predict(m_list[[max_r]], re.form = ~ 0),
               fitted_value_re = predict(m_list[[max_r]], re.form = formula(formula_string_re), random.only = T),
               sigma_u = (1/sqrt(max_repetition))*sqrt(VarCorr(m_list[[max_r]])$PEDIGREE_NAME[1,1]),
               sigma_e = sigma(m_list[[max_r]]),
               ) %>%
        left_join(df_ranef, by = "PEDIGREE_NAME")
  }
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
    df_poly = data.frame(rep(1/sqrt(max_repetition), max_repetition))
    colnames(df_poly) = c("poly0")
  }
  rownames(df_poly) = 1:max_repetition
  df_poly %>% as.matrix() %>% return()
}

h2 = function(max_repetition, model) {
  sigma_u = sqrt(max_repetition)*sqrt(VarCorr(model)$PEDIGREE_NAME[1,1])
  sigma_e = sigma(model)
  return(sigma_u^2/(sigma_u^2 + sigma_e^2))
}

by_trait_fit_lme_poly = function(df_input, trait_names, ...) {
  df_model = lapply(trait_names, function(trt) {
    fitted_model = df_input %>%
      filter(OBSRVTN_REF_CD == trt) %>%
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>%
      fit_lme_poly(...)
    fitted_model$df %>%
      select(OBSRVTN_REF_CD, TRAIT_VALUE, GROWSEASON, FIELD_NAME, PEDIGREE_NAME, PLOT_BID, REPETITION,
             full_harvest_repetition, current_harvest_repetition, skip,
             fitted_value, fitted_value_fe, fitted_value_re, sigma_u, sigma_e, condval_u, condsd_u) %>%
      return()
  })
  df_model = Reduce(rbind, df_model)
  return(df_model)
}
