library(glmmTMB)

# Function to fit a multi-trait model
fit_multitrait = function(df, trait_names, has_growseason = T, has_field_name = T){
  # Preprocessing the data
  df = df %>%
    select(OBSRVTN_REF_CD, PEDIGREE_NAME, TRAIT_VALUE, GROWSEASON, FIELD_NAME, PLOT_BID, REPETITION) %>%
    filter(OBSRVTN_REF_CD %in% trait_names) %>%
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>%
    mutate(OBSRVTN_REF_CD = as.factor(OBSRVTN_REF_CD))
  
  # Constructing the formula string based on the input
  growseason_string = ifelse(has_growseason, "+ GROWSEASON", "")
  field_name_string = ifelse(has_field_name, "+ FIELD_NAME", "")
  formula_string = paste("TRAIT_VALUE ~ 1 ",  growseason_string, field_name_string, " + (1 | PEDIGREE_NAME)", collapse = "")
  
  # Fitting the model for each trait and storing the results
  df_fitted = data.frame()
  for(trt in trait_names) {
    df_trt = df %>% filter(OBSRVTN_REF_CD == trt)
    m_trt = glmmTMB(formula(formula_string),
                    family = gaussian,
                    data = df_trt,
                    REML=TRUE)
    
    # Extracting fixed and random effects
    df_trt = df_trt %>% mutate(fe = predict(m_trt, re.form = ~ 0))
    df_re = ranef(m_trt)$cond$PEDIGREE_NAME %>%
      rownames_to_column("PEDIGREE_NAME") %>%
      rename(re = "(Intercept)")
    df_trt = df_trt %>% left_join(df_re, by = "PEDIGREE_NAME")
    df_trt = df_trt %>% mutate(fitted_value = fe + re)
    df_trt = df_trt %>% mutate(sigma_u = sqrt(VarCorr(m_trt)$cond$PEDIGREE_NAME[1,1]), sigma_e = sigma(m_trt))
    df_fitted = rbind(df_fitted, df_trt)
  }
  df_fitted
}

# Function to rank the principal components (PCs) of the fitted values
rank_pc = function(df_fitted, trait_names, fitted_col = "re", trait_rank_combination = NULL, which_quality_trait = 1) {
  # Grouping and summarizing the data
  df_fitted = df_fitted %>%
    group_by(OBSRVTN_REF_CD, PEDIGREE_NAME) %>%
    summarise(mean_value = mean(!! sym(fitted_col)), .groups="drop") %>%
    pivot_wider(names_from = OBSRVTN_REF_CD, values_from = mean_value) %>%
    ungroup()
  
  # Performing principal component analysis (PCA) and ranking
  if(is.null(trait_rank_combination)) {
    pc = prcomp(df_fitted %>% select(all_of(trait_names)) %>% data.matrix(), scale=TRUE)
    sign = (pc$rotation[which_quality_trait,1] > 0)*2 - 1
    df_fitted %>% mutate(rank = rank(sign*predict(pc)[,1]))
  } else {
    rk = rank((df_fitted %>% select(all_of(trait_names)) %>% data.matrix() %>% scale()) %*% trait_rank_combination)
    df_fitted %>% mutate(rank = rk)
  }
}

# Function to rank multi-trait data across different repetitions
rank_multitrait = function(df, trait_names, has_growseason, has_field_name, range_n_repetition, fitted_col = "re", trait_rank_combination = NULL, which_quality_trait = 1) {
  df_rank = data.frame()
  for(r in range_n_repetition){
    df_r = df %>% filter(REPETITION <= r)
    df_r_fitted = fit_multitrait(df_r, trait_names, has_growseason, has_field_name)
    df_rank = rbind(df_rank,
                    rank_pc(df_r_fitted, trait_names, fitted_col = fitted_col, trait_rank_combination = trait_rank_combination, which_quality_trait = which_quality_trait) %>%
                      mutate(current_harvest_repetition = as.integer(r))
    )
  }
  return(df_rank)
}

# Function to simulate multi-trait data based on the fitted values
simulate_multitrait = function(df_fitted){
  df_fitted %>%
    mutate(TRAIT_VALUE = fitted_value + rnorm(n(), sd = sigma_e)) %>%
    select(OBSRVTN_REF_CD, PEDIGREE_NAME, TRAIT_VALUE, GROWSEASON, FIELD_NAME, PLOT_BID, REPETITION)
}

# Function to perform bootstrapping by resampling on PLOT_BID
boot_on_group = function(df, g, prop = 0.9, replace=F) {
  df %>%
    nest(data = !all_of(g)) %>%
    slice_sample(prop=prop, replace=replace) %>%
    group_by(PLOT_BID) %>%
    mutate(id = 1:n()) %>%
    arrange(as.character(PLOT_BID)) %>%
    unite(PLOT_BID, PLOT_BID, id, sep="_") %>%
    unnest(data) %>%
    select(OBSRVTN_REF_CD, TRAIT_VALUE, GROWSEASON, FIELD_NAME, PEDIGREE_NAME, PLOT_BID, REPETITION)
}

# Function to simulate larger repetitions for multi-trait data
simulate_multitrait_more_repetition = function(df_fitted, more_repetition) {
  # Simulating more repetitions and combining with the original data
  rbind(
    df_fitted %>% select(-fe, -re, -fitted_value, -sigma_u, -sigma_e),
    df_fitted %>% mutate(max_repetition = max(REPETITION)) %>%
      distinct(OBSRVTN_REF_CD, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID, max_repetition) %>%
      rowwise() %>%
      mutate(REPETITION = list((max_repetition+1):(max_repetition+more_repetition))) %>%
      unnest(REPETITION) %>%
      left_join(df_fitted %>%
                  distinct(OBSRVTN_REF_CD, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, fitted_value, sigma_e),
                by = c("OBSRVTN_REF_CD", "PEDIGREE_NAME", "GROWSEASON", "FIELD_NAME")) %>%
      mutate(TRAIT_VALUE = fitted_value + rnorm(n(), sd = sigma_e)) %>% select(-fitted_value, -sigma_e, -max_repetition)
  ) %>% arrange(OBSRVTN_REF_CD, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, REPETITION)
}

simulate_more_repetition = function(df_fitted, more_repetition) {
  rbind(
    df_fitted %>% select(-fe, -re, -fitted_value, -sigma_u, -sigma_e),
    df_fitted %>% mutate(max_repetition = max(REPETITION)) %>%
      distinct(OBSRVTN_REF_CD, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID, max_repetition) %>%
      rowwise() %>%
      mutate(REPETITION = list((max_repetition+1):(max_repetition+more_repetition))) %>%
      unnest(REPETITION) %>%
      left_join(df_fitted %>%
                  distinct(OBSRVTN_REF_CD, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, fitted_value, sigma_e),
                by = c("OBSRVTN_REF_CD", "PEDIGREE_NAME", "GROWSEASON", "FIELD_NAME")) %>%
      mutate(TRAIT_VALUE = fitted_value + rnorm(n(), sd = sigma_e)) %>% select(-fitted_value, -sigma_e, -max_repetition)
  ) %>% arrange(OBSRVTN_REF_CD, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, REPETITION)
}

# Function to obtain uncertainty in metrics for multi-trait data using bootstrapping
smart_harvest_planning_multitrait = function(
  df, # Data frame containing the observations
  trait_names, # Vector of trait names to be analyzed
  range_n_repetition, # Range of repetitions to be considered
  metric_fn_list, # List of metric functions to be applied
  selection_intensity, # Selection intensity to be used in the metrics
  has_growseason = T, # Whether the growing season is considered
  has_field_name = T, # Whether the field name is considered
  fitted_col = "re", # Column name for the fitted values
  trait_rank_combination = NULL, # Combination of traits for ranking
  which_quality_trait = 1, # Quality trait to be considered
  nsim = 100, # Number of simulations for bootstrapping
  bootstrap_method = NULL, # Method for bootstrapping
  q_probs = c(0.1, 0.5, 0.9), # Quantile probabilities for summarizing metrics
  read_rds_file = NULL, # Optional path to read RDS file
  save_rds_file = NULL # Optional path to save RDS file
){
  
  # Initializing a data frame to store simulated ranks
  df_sim_rank = data.frame()
  
  # Creating a progress bar for bootstrapping
  pb <- progress_bar$new(
    format = "Bootstrapping [:bar] :percent eta: :eta elapsed :elapsed",
    total = nsim, clear = FALSE, width= 60)
  
  # Ranking the actual data
  df_rank = rank_multitrait(df, trait_names, has_growseason, has_field_name, range_n_repetition, fitted_col = fitted_col, trait_rank_combination = trait_rank_combination, which_quality_trait = which_quality_trait)
  
  # Extracting the truth rank
  rk_truth = df_rank %>%
    select(PEDIGREE_NAME, rank, current_harvest_repetition) %>%
    pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
    pull(get(as.character(max(range_n_repetition))))
  
  # Initializing a data frame to store metrics
  df_metric = data.frame()
  
  # Looping through metrics and selection intensities to compute metrics
  for(metric_name in names(metric_fn_list)) {
    for(s in selection_intensity) {
      df_metric = rbind(df_metric,
                        df_rank %>%
                          select(PEDIGREE_NAME, rank, current_harvest_repetition) %>%
                          pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
                          summarise(across(all_of(as.character(range_n_repetition)), ~ metric_fn_list[[metric_name]](s, .x, rk_truth)), .groups="drop") %>%
                          pivot_longer(all_of(as.character(range_n_repetition)), names_to = "REPETITION", values_to = "metric_value") %>% mutate(REPETITION = as.integer(REPETITION)) %>%
                          group_by(REPETITION) %>%
                          summarise(q = list(quantile(metric_value, probs = q_probs, na.rm=T)), .groups="drop") %>% 
                          unnest_wider(q) %>%
                          mutate(selection_intensity = s, metric = metric_name, note = "relative to full harvest")
      )
    }}
  
  # Fitting the multi-trait model
  df_fitted = fit_multitrait(df, trait_names, has_growseason, has_field_name)
  
  # Looping through bootstrapping simulations
  for(b in 1:nsim) {
    pb$tick()
    
    # Simulating data based on the bootstrap method
    if(is.null(bootstrap_method)) {
      df_sim = simulate_multitrait(df_fitted)
    } else {
      bootstrap_method_par = unlist(strsplit(bootstrap_method, split = "_"))
      if(bootstrap_method_par[1] == "resample"){
        df_sim = boot_on_group(df_fitted, "PLOT_BID",
                               prop = as.numeric(bootstrap_method_par[2]),
                               replace = as.logical(bootstrap_method_par[3]))
      }
    }
    
    # Ranking the simulated data and storing the results
    df_sim_rank = bind_rows(
      df_sim_rank,
      rank_multitrait(df_sim, trait_names, has_growseason, has_field_name, range_n_repetition, fitted_col = fitted_col, trait_rank_combination = trait_rank_combination, which_quality_trait = which_quality_trait) %>%
        mutate(b_sample = b)
    )
  }
  
  # Reading previous bootstrap results if provided
  if(!is.null(read_rds_file)) {
    df_sim_rank_read = readRDS(read_rds_file)
    last_b = df_sim_rank_read %>% pull(b_sample) %>% max()
  } else {
    df_sim_rank_read = data.frame()
    last_b = 0
  }
  # Combining previous and current bootstrap results
  df_sim_rank = bind_rows(df_sim_rank_read, df_sim_rank %>% mutate(b_sample = last_b + b_sample))
  
  # Creating a data frame that contains the truth rank for each bootstrap sample
  df_rank_truth_by_sim = df_sim_rank %>%
    select(PEDIGREE_NAME, current_harvest_repetition, b_sample) %>%
    left_join(df_rank %>% select(PEDIGREE_NAME, current_harvest_repetition, rank), by = c("PEDIGREE_NAME", "current_harvest_repetition")) %>%
    group_by(current_harvest_repetition, b_sample) %>%
    mutate(rank = rank(rank)) %>%
    select(PEDIGREE_NAME, rank, current_harvest_repetition, b_sample) %>%
    pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
    arrange(b_sample, PEDIGREE_NAME) %>%
    select(PEDIGREE_NAME, b_sample, !! sym(as.character(max(range_n_repetition)))) %>%
    rename(truth_rank=as.character(max(range_n_repetition)))
  
  # Creating a new progress bar for computing metrics
  pb <- progress_bar$new(
    format = "Computing Metric [:bar] :percent eta: :eta elapsed :elapsed",
    total = length(names(metric_fn_list))*length(selection_intensity), clear = FALSE, width= 60)
  
  # Initializing data frames to store metrics relative to truth and full harvest
  df_sim_metric_truth = data.frame()
  df_sim_metric_converge = data.frame()
  
  # Computing metrics relative to truth and full harvest for bootstrapped data
  # Looping through metric functions and selection intensities
  for(metric_name in names(metric_fn_list)) {
    for(s in selection_intensity) {
      pb$tick()
      
      # Computing metrics relative to truth for bootstrapped data
      df_sim_metric_truth = rbind(df_sim_metric_truth,
                                  df_sim_rank %>%
                                    select(PEDIGREE_NAME, rank, current_harvest_repetition, b_sample) %>%
                                    pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
                                    left_join(df_rank_truth_by_sim, by = c("PEDIGREE_NAME", "b_sample")) %>%
                                    group_by(b_sample) %>%
                                    summarise(across(all_of(as.character(range_n_repetition)), ~ metric_fn_list[[metric_name]](s, .x, truth_rank)), .groups="drop") %>%
                                    pivot_longer(all_of(as.character(range_n_repetition)), names_to = "REPETITION", values_to = "metric_value") %>% mutate(REPETITION = as.integer(REPETITION)) %>%
                                    group_by(REPETITION) %>%
                                    summarise(q = list(quantile(metric_value, probs = q_probs, na.rm=T)), .groups="drop") %>% 
                                    unnest_wider(q) %>%
                                    mutate(selection_intensity = s, metric = metric_name, note = "relative to truth")
      )
      
      # Computing metrics relative to full harvest for bootstrapped data
      df_sim_metric_converge = rbind(df_sim_metric_converge,
                                     df_sim_rank %>%
                                       select(PEDIGREE_NAME, rank, current_harvest_repetition, b_sample) %>%
                                       pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
                                       group_by(b_sample) %>%
                                       summarise(across(all_of(as.character(range_n_repetition)), ~ metric_fn_list[[metric_name]](s, .x, get(as.character(max(range_n_repetition))))), .groups="drop") %>%
                                       pivot_longer(all_of(as.character(range_n_repetition)), names_to = "REPETITION", values_to = "metric_value") %>% mutate(REPETITION = as.integer(REPETITION)) %>%
                                       group_by(REPETITION) %>%
                                       summarise(q = list(quantile(metric_value, probs = q_probs, na.rm=T)), .groups="drop") %>% 
                                       unnest_wider(q) %>%
                                       mutate(selection_intensity = s, metric = metric_name, note = "relative to full harvest")
      )
    }}
  
  # Saving the bootstrap results if a path is provided
  if(!is.null(save_rds_file)) {
    saveRDS(df_sim_rank, save_rds_file)
  }
  
  return(list(rank_actual=df_rank, rank_bootstrap=df_sim_rank,
              metric_actual=df_metric, metric_truth_bootstrap=df_sim_metric_truth, metric_converge_bootstrap=df_sim_metric_converge))
}
