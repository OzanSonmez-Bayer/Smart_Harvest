
library(tidyverse)
library(stringr)
library(lme4)

load("/mnt/Data/Pepper_data.RData")
df_pepper_netherlands <- df_pepper

df_pepper_netherlands = df_pepper_netherlands %>%
  mutate_at(c("FIELD_NAME", "PEDIGREE_NAME", "PLOT_BID"), as.factor)

df_pepper_mexico = read.csv("/mnt/Data/Pepper_Zeb_Digital_Pheno_Data.csv")
df_pepper_mexico = df_pepper_mexico %>%
  mutate_at(c("FIELD_NAME", "PEDIGREE_NAME", "PLOT_BID"), as.factor)


df_pepper_mexico %>%
  filter(OBSRVTN_REF_CD == "FQUAL") %>%
  select(TRAIT_VALUE, PEDIGREE_NAME, PLOT_BID, GROWSEASON, FIELD_NAME, REPETITION)




generate_poly_coef = function(u_PEDIGREE_NAME, poly_max_degree, poly_sigma) {
  matrix(rnorm(length(u_PEDIGREE_NAME)*(poly_max_degree+1), sd = poly_sigma),
                           nrow = length(u_PEDIGREE_NAME),
                           ncol=poly_max_degree+1,
                           dimnames = list(u_PEDIGREE_NAME, paste0("poly", 0:poly_max_degree)),
                           byrow = T)
}

generate_poly_fitted_value = function(max_repetition, poly_max_degree, v_PEDIGREE_NAME) {
  (generate_poly(max_repetition, poly_max_degree) %*% t(v_PEDIGREE_NAME)) %>%
    data.frame() %>%
    rownames_to_column("REPETITION") %>%
    mutate(REPETITION = as.integer(REPETITION)) %>%
    pivot_longer(cols = -c("REPETITION"), names_to = "PEDIGREE_NAME", values_to = "FITTED_VALUE") %>%
    return()
}





rank_random_effect = function(df) {
  rank_df = df %>%
    select(PEDIGREE_NAME, REPETITION, fitted_value_re) %>%
    arrange(PEDIGREE_NAME, REPETITION) %>%
    group_by(PEDIGREE_NAME) %>%
    summarise(mean_re = mean(fitted_value_re), .groups = "drop") %>%
    mutate(rank = rank(mean_re))
  return(rank_df)
}

rank_col = function(df, col) {
  rank_df = df %>%
    select(PEDIGREE_NAME, REPETITION, {{col}} ) %>%
    arrange(PEDIGREE_NAME, REPETITION) %>%
    group_by(PEDIGREE_NAME) %>%
    summarise(mean = mean({{col}}), .groups = "drop") %>%
    mutate(rank = rank(mean))
  return(rank_df)
}






rank_col(fitted_model$d, fitted_value_re)

fitted_model$df %>% select(PEDIGREE_NAME, REPETITION, "fitted_value_re") %>% arrange(PEDIGREE_NAME, REPETITION) %>%
  group_by(PEDIGREE_NAME) %>% summarise(mean = mean("fitted_value_re"), .groups = "drop")

  
set.seed(1234)

u_PEDIGREE_NAME = paste0("PE", str_pad(1:80, 2, pad = "0"))
u_GROWSEASON = paste0("2023:", str_pad(1:2, 2, pad = "0"))
u_REPETITION = 1:52
u_FIELD_NAME = paste0("FI", str_pad(1:3, 2, pad = "0"))


df_plot_bid = expand.grid(PEDIGREE_NAME = u_PEDIGREE_NAME,
                          GROWSEASON = u_GROWSEASON,
                          FIELD_NAME = u_FIELD_NAME) %>%
  arrange(FIELD_NAME, GROWSEASON, PEDIGREE_NAME) %>%
  mutate(PLOT_BID = paste0("PL", str_pad(1:n(), 4, pad = "0")))

df_sim_design = expand.grid(PEDIGREE_NAME = u_PEDIGREE_NAME,
                            GROWSEASON = u_GROWSEASON,
                            FIELD_NAME = u_FIELD_NAME,
                            REPETITION = u_REPETITION) %>%
  arrange(FIELD_NAME, GROWSEASON, PEDIGREE_NAME) %>%
  left_join(df_plot_bid, by = c("FIELD_NAME", "GROWSEASON", "PEDIGREE_NAME"))

v_GROWSEASON = setNames(rnorm(length(u_GROWSEASON), sd = 0.2),
                        u_GROWSEASON)

v_FIELD_NAME = setNames(rnorm(length(u_FIELD_NAME), mean = 4, sd = 0.2),
                        u_FIELD_NAME)
v_PLOT_BID = setNames(rnorm(n_distinct(df_plot_bid$PLOT_BID), sd = 0.1),
                      unique(df_plot_bid$PLOT_BID))

poly_max_degree = 3
poly_sigma = c(5.0, 0.5, 0.01, 0.01)
v_PEDIGREE_NAME = generate_poly_coef(u_PEDIGREE_NAME, poly_max_degree, poly_sigma)
f_PEDIGREE_NAME = generate_poly_fitted_value(max(u_REPETITION), poly_max_degree, v_PEDIGREE_NAME)

df_sim_reference = df_sim_design %>%
  left_join(f_PEDIGREE_NAME %>% rename(FITTED_VALUE_PEDIGREE_NAME = FITTED_VALUE), by = c("PEDIGREE_NAME", "REPETITION")) %>%
  left_join(v_GROWSEASON %>% enframe(name = "GROWSEASON", value = "FITTED_VALUE_GROWSEASON"), by = "GROWSEASON") %>%
  left_join(v_FIELD_NAME %>% enframe(name = "FIELD_NAME", value = "FITTED_VALUE_FIELD_NAME"), by = "FIELD_NAME") %>% 
  left_join(v_PLOT_BID %>% enframe(name = "PLOT_BID", value = "FITTED_VALUE_PLOT_BID"), by = "PLOT_BID") %>%
  mutate(FITTED_VALUE_TRAIT_VALUE = FITTED_VALUE_PEDIGREE_NAME + FITTED_VALUE_GROWSEASON + FITTED_VALUE_FIELD_NAME + FITTED_VALUE_PLOT_BID) %>%
  mutate(TRAIT_VALUE = pmin(pmax(round(FITTED_VALUE_TRAIT_VALUE + rnorm(n(), sd = 0.1), digits = 0), 0), 9) )

df_sim = df_sim_reference %>% select(TRAIT_VALUE, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID, REPETITION)


df_sim %>% left_join(generate_poly(max(u_REPETITION), 3) %>% data.frame() %>% rownames_to_column("REPETITION") %>% mutate(REPETITION = as.integer(REPETITION)), by = "REPETITION")



fitted_model = fit_lme_poly(df_sim, 3)

fitted_model@frame %>% group_by(GROWSEASON, FIELD_NAME) %>% summarise(fitted_value = mean(fitted_value_fe), .groups="drop")
fitted_model@frame %>% group_by(PEDIGREE_NAME, REPETITION) %>% summarise(fitted_value = mean(fitted_value_re), .groups="drop")

fitted_model@frame

fitted_model = fit_lme_poly(df_sim %>% slice_sample(prop = 0.95), 3)

df_new = df_sim %>%
  left_join(generate_poly(max(u_REPETITION), 3) %>%
              data.frame() %>%
              rownames_to_column("REPETITION") %>%
              mutate(REPETITION = as.integer(REPETITION)),
            by = "REPETITION")

df_new %>% mutate(p = predict(fitted_model$model, newdata = df_new))


rank_random_effect(fitted_model$df) %>% print(n=80)



metric_precision = function(selection_intensity, predict_rank, reference_rank){
  predict_rank %>% mutate(selected = rank <= selection_intensity * length(rank)) %>%
    left_join(reference_rank %>% mutate(selected = rank <= selection_intensity * length(rank)), by = "PEDIGREE_NAME") %>%
    summarise(precision = sum(selected.x * selected.y)/sum(selected.x)) %>%
    pull(precision)
}




df_sim_reference_rank = df_sim_reference %>%
  group_by(PEDIGREE_NAME, REPETITION) %>%
  summarise(mean_re = mean(FITTED_VALUE_PEDIGREE_NAME), .groups = "drop") %>%
  group_by(PEDIGREE_NAME) %>%
  summarise(mean_re = mean(mean_re), .groups = "drop") %>%
  mutate(rank = rank(mean_re))

fitted_model = fit_lme_poly(df_sim %>% filter(REPETITION <= 10) %>% slice_sample(prop = 0.95), 3)



plot(seq(0.1, 1,by=0.1), sapply(seq(0.1, 1,by=0.1), function(p) metric_precision(p, rank_random_effect(fitted_model$df), df_sim_reference_rank)), type = "l")
plot(seq(0.1, 1,by=0.1), sapply(seq(0.1, 1,by=0.1), function(p) metric_ndcg(p, rank_random_effect(fitted_model$df), df_sim_reference_rank)), type = "l")


rank_random_effect(fitted_model$df) %>% left_join(df_sim_reference_rank, by = "PEDIGREE_NAME") %>%
  ggplot() + aes(x = rank.x, y = rank.y) + geom_point()



metric_ndcg = function(selection_intensity, predict_rank, reference_rank){
  predict_rank %>% filter(rank <= selection_intensity * length(rank)) %>%
    arrange(rank) %>%
    left_join(reference_rank %>% mutate(score = rank <= selection_intensity * length(rank)) %>% select(PEDIGREE_NAME, score), by = "PEDIGREE_NAME") %>%
    summarise(nDCG = sum(score/log2(rank + 1))/sum(1/log2(rank + 1)) ) %>%
    pull(nDCG)
}


# Function to calculate DCG
calc_DCG <- function(scores) {
  DCG = sum(scores / log2(1:length(scores) + 1))
  return(DCG)
}

metric_ndcg(0.2, rank_random_effect(fitted_model$df), df_sim_reference_rank)

metric_spearman(0.2, rank_random_effect(fitted_model$df), df_sim_reference_rank)

plot(sapply(seq(0.1, 1, by = 0.05), function(p) metric_spearman(p, rank_random_effect(fitted_model$df), df_sim_reference_rank)))


rank_random_effect(fitted_model$df) %>% left_join(df_sim_reference_rank, by = "PEDIGREE_NAME")













library(tidyverse)
library(stringr)
library(lme4)

library(foreach)

setwd("/repos/Smart_Harvest/Modules/Test")
source("../Model/Model.R")
source("../Model/rank.R")
source("../Metrics/precision.R")
source("../Metrics/ndcg.R")
source("../Metrics/jaccard_similarity.R")
source("../Metrics/spearman.R")

library(parallel)
detectCores()

test_model = fit_lme_poly(df_sim, 3, 5:52)


fitted_model = fit_lme_poly(df_sim, 3, 5:52)

tibble(fitted_model$df)

x <- foreach(i=1:3) %do% sqrt(i)
x

library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
output_parallel = parallel_fit_lme_poly2(df_sim, 3, 5:52)



output_parallel[[1]]



get_precision = function(df, sel_seq) {
  df_rank = df %>%
    group_by(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip) %>%
    summarise(mean_re = mean(fitted_value_re), .groups = "drop") %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    mutate(rank = rank(mean_re)) %>% ungroup() %>%
    arrange(current_harvest_repetition) %>%
    select(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip, rank)
  
  Reduce(rbind, lapply(sel_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      summarise(selection_intensity = s, across("5":"52", ~ metric_precision_rank(s, ., `52`)))
  })) %>%
    group_by(full_harvest_repetition, skip, selection_intensity) %>%
    pivot_longer(cols = "5":"52", names_to = "current_harvest_repetition", values_to = "precision") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition))
}

get_precision(fitted_model$df, seq(0.1, 0.5, by = 0.1)) %>% print(n = 100)



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

rank_to_precision = function(df_rank, s_seq) {
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                precision = map(data,
                                function(d) summarise(d,
                                                      across(-PEDIGREE_NAME, ~ metric_precision_rank(s, ., d[[as.character(full_harvest_repetition)]]))
                                )
                ), .groups="drop") %>%
      unnest(cols = precision)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, selection_intensity), names_to = "current_harvest_repetition", values_to = "precision") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}

rank_to_ndcg = function(df_rank, s_seq) {
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                precision = map(data,
                                function(d) summarise(d,
                                                      across(-PEDIGREE_NAME, ~ metric_ndcg(s, ., d[[as.character(full_harvest_repetition)]]))
                                )
                ), .groups="drop") %>%
      unnest(cols = precision)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, selection_intensity), names_to = "current_harvest_repetition", values_to = "precision") %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}


fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_precision(seq(0.1, 0.5, by = 0.1)) %>% print(n = 100)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_ndcg(seq(0.1, 0.5, by = 0.1)) %>% print(n = 500)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_jaccard_similarity(seq(0.1, 0.5, by = 0.1)) %>% print(n = 500)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>% print(n = 100, width=Inf)

fitted_model$df %>%
  group_by(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip) %>%
  summarise(mean_re = mean(fitted_value_re), .groups = "drop") %>%
  left_join(fitted_model$df %>% fitted_to_rank(fitted_value_re),
            by = c("PEDIGREE_NAME", "full_harvest_repetition", "current_harvest_repetition", "skip")) %>%
  group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
  mutate(selected = rank <= 0.1*length(rank), mean_re_all = mean(mean_re)) %>%
  filter(selected == TRUE) %>%
  summarise(mean_re_selected = mean(mean_re), mean_re_all = mean(mean_re_all)) %>%
  print(n = 200)
  


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

fitted_to_genetic_gain(fitted_model$df, seq(0.1, 0.5, by=0.1)) %>% print(n = 300)



%>%
  ggplot() + aes(x = current_harvest_repetition, y = normalized_genetic_gain) + geom_line()











%>% arrange(current_harvest_repetition) %>% print(n = 100)



%>%
  pivot_wider(names_from = current_harvest_repetition, values_from = mean_re)

%>%
  pivot_wider(names_from = current_harvest_repetition, values_from = rank)
  

metric_ndcg(0.5, c(1,2,3,4,5), c(1,3,4,5,2))

c(1,3,4,5,2) <= 0.5*5
c(1,2,3,4,5) <= 0.5*5


hello_list = list()
hello_list[[10]] = c(1,2,3)
hello_list[[11]] = "Hello World"


names(hello_list) = NULL


str(hello_list)

helnames(list(a = 1, b = 2))

for(m in fitted_model$model[1:10]){
  print(summary(m))
}

sapply(fitted_model$model, is.null)

fitted_model$model[[52]]@theta
sapply(fitted_model$model[5:52], function(m) m@theta)

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

df_rank = fitted_model$df %>%
  group_by(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip) %>%
  summarise(mean_value = mean(fitted_value_re), .groups = "drop") %>%
  group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
  mutate(rank = rank(mean_value)) %>% ungroup() %>%
  arrange(current_harvest_repetition) %>%
  select(PEDIGREE_NAME, full_harvest_repetition, current_harvest_repetition, skip, mean_value, rank)

Reduce(rbind, lapply(seq(0.1, 0.5, by = 0.1), function(selection_intensity) {
  df_rank %>% left_join(df_rank %>%
                          filter(current_harvest_repetition == full_harvest_repetition) %>%
                          rename(full_harvest_rank = rank, full_harvest_mean_value = mean_value) %>%
                          select(PEDIGREE_NAME, full_harvest_repetition, skip, full_harvest_rank, full_harvest_mean_value),
                        by = c("PEDIGREE_NAME", "full_harvest_repetition", "skip")) %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    filter( (rank <= selection_intensity*length(rank)) | (full_harvest_rank <= selection_intensity*length(full_harvest_rank)) ) %>%
    summarise(selection_intensity = selection_intensity, cor = cor(mean_value, full_harvest_mean_value), .groups = "drop")
}))



  


rank_to_metric = function(df_rank, s_seq, metric_name, metric_fun, benchmark_harvest_repetition = NULL) {
  if(is.null(benchmark_harvest_repetition)) {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = full_harvest_repetition)
  } else {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = benchmark_harvest_repetition)
  }
  
  Reduce(rbind, lapply(s_seq, function(s) {
    df_rank %>%
      group_by(full_harvest_repetition, skip, benchmark_harvest_repetition) %>%
      pivot_wider(names_from = current_harvest_repetition, values_from = rank) %>%
      nest() %>%
      summarise(selection_intensity = s,
                metric_name = map(data,
                           function(d) summarise(d,
                                                 across(-PEDIGREE_NAME, ~ metric_fun(s, ., d[[as.character(benchmark_harvest_repetition)]]))
                           )
                ), .groups="drop") %>%
      unnest(cols = metric_name)
  })) %>%
    pivot_longer(cols = -c(full_harvest_repetition, skip, benchmark_harvest_repetition, selection_intensity), names_to = "current_harvest_repetition", values_to = metric_name) %>%
    mutate(current_harvest_repetition = as.integer(current_harvest_repetition)) %>%
    return()
}


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






df_sim %>% head()



df_sim %>% slice_sample(prop = 0.97) %>%
  group_by(PEDIGREE_NAME, GROWSEASON, FIELD_NAME) %>%
  summarise(n = n(), .groups="drop") %>% 
  filter(n > 1) %>% print(n = 500)



df_sim %>% slice_sample(prop = 0.97) %>%
  group_by(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
  summarise(n = n(), .groups="drop") %>% print(n = 500)


df_sim %>% slice_sample(prop = 0.97) %>%
  distinct(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
  arrange(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
  mutate(count = 52) %>%
  uncount(count) %>%
  group_by(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
  mutate(REPETITION = 1:52) %>% print(n = 500)



df_sim %>% slice_sample(prop = 0.97) %>%
  distinct(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
  arrange(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
  mutate(REPETITION = list(1:52)) %>% unnest(cols = REPETITION) %>%
  print(n = 500)




full_design = function(df, max_repitition) {
  df %>%
    distinct(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
    arrange(PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID) %>%
    mutate(REPETITION = list(1:max_repitition)) %>%
    unnest(cols = REPETITION) %>%
    return()
}

full_design(df_sim %>% slice_sample(prop = 0.97), 3)


predict(fitted_model$model[[52]], 


        