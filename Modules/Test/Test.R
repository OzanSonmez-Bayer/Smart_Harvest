# test the workflow for a given data set

library(tidyverse)
library(stringr)
library(lme4)
library(progress)

setwd("/repos/Smart_Harvest/Modules/Test")
source("../Model/Model.R")
source("../Model/rank.R")
source("../Metrics/metric.R")
source("../Metrics/precision.R")
source("../Metrics/ndcg.R")
source("../Metrics/jaccard_similarity.R")
source("../Metrics/spearman.R")
source("../Metrics/kendall.R")
source("../Metrics/marginal_rank_change.R")
source("../Metrics/genetic_gain.R")
source("../Metrics/genetic_correlation.R")
source("../Metrics/heritability.R")

# Functions for generating random polynomials.
#

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

simulate_data = function(poly_sigma, 
         u_PEDIGREE_NAME,
         u_GROWSEASON,
         u_FIELD_NAME,
         u_REPETITION,
         sigma_GROWSEASON,
         simga_FIELD_NAME,
         sigma_PLOT_BID,
         mean_FIELD,
         sigma_noise,
         seed = 1234) {
  set.seed(seed)
  poly_max_degree = length(poly_sigma) - 1
  df_plot_bid = expand.grid(PEDIGREE_NAME = u_PEDIGREE_NAME,
                            GROWSEASON = u_GROWSEASON,
                            FIELD_NAME = u_FIELD_NAME) %>%
    arrange(FIELD_NAME, GROWSEASON, PEDIGREE_NAME) %>%
    mutate(count = sample.int(3, size = n(), prob = c(0.8, 0.15, 0.05), replace=T)) %>%
    uncount(count) %>%
    mutate(PLOT_BID = paste0("PL", str_pad(1:n(), 4, pad = "0")))

  df_sim_design = expand.grid(PEDIGREE_NAME = u_PEDIGREE_NAME,
                              GROWSEASON = u_GROWSEASON,
                              FIELD_NAME = u_FIELD_NAME,
                              REPETITION = u_REPETITION) %>%
    arrange(FIELD_NAME, GROWSEASON, PEDIGREE_NAME) %>%
    left_join(df_plot_bid, by = c("FIELD_NAME", "GROWSEASON", "PEDIGREE_NAME"))
  
  v_GROWSEASON = setNames(rnorm(length(u_GROWSEASON), sd = sigma_GROWSEASON),
                          u_GROWSEASON)
  
  v_FIELD_NAME = setNames(rnorm(length(u_FIELD_NAME), mean = mean_FIELD, sd = simga_FIELD_NAME),
                          u_FIELD_NAME)
  v_PLOT_BID = setNames(rnorm(n_distinct(df_plot_bid$PLOT_BID), sd = sigma_PLOT_BID),
                        unique(df_plot_bid$PLOT_BID))
  
  v_PEDIGREE_NAME = generate_poly_coef(u_PEDIGREE_NAME, poly_max_degree, poly_sigma)
  f_PEDIGREE_NAME = generate_poly_fitted_value(max(u_REPETITION), poly_max_degree, v_PEDIGREE_NAME)
  
  df_sim_reference = df_sim_design %>%
    left_join(f_PEDIGREE_NAME %>% rename(FITTED_VALUE_PEDIGREE_NAME = FITTED_VALUE), by = c("PEDIGREE_NAME", "REPETITION")) %>%
    left_join(v_GROWSEASON %>% enframe(name = "GROWSEASON", value = "FITTED_VALUE_GROWSEASON"), by = "GROWSEASON") %>%
    left_join(v_FIELD_NAME %>% enframe(name = "FIELD_NAME", value = "FITTED_VALUE_FIELD_NAME"), by = "FIELD_NAME") %>% 
    left_join(v_PLOT_BID %>% enframe(name = "PLOT_BID", value = "FITTED_VALUE_PLOT_BID"), by = "PLOT_BID") %>%
    mutate(FITTED_VALUE_TRAIT_VALUE = FITTED_VALUE_PEDIGREE_NAME + FITTED_VALUE_GROWSEASON + FITTED_VALUE_FIELD_NAME + FITTED_VALUE_PLOT_BID) %>%
    mutate(TRAIT_VALUE = pmin(pmax(round(FITTED_VALUE_TRAIT_VALUE + rnorm(n(), sd = sigma_noise), digits = 0), 1), 9))
  
  df_sim_reference_rank = df_sim_reference %>%
    group_by(PEDIGREE_NAME, REPETITION) %>%
    summarise(mean_re = mean(FITTED_VALUE_PEDIGREE_NAME), .groups = "drop") %>%
    group_by(PEDIGREE_NAME) %>%
    summarise(mean_re = mean(mean_re), .groups = "drop") %>%
    mutate(rank = rank(mean_re))
  
  df_sim = df_sim_reference %>%
    select(TRAIT_VALUE, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID, REPETITION)
  
  return(list(df_sim = df_sim, df_sim_reference = df_sim_reference, df_sim_reference_rank = df_sim_reference_rank))
}


# Preparing simulated data set
#

sigma_GROWSEASON = 0.2
simga_FIELD_NAME = 0.2
sigma_PLOT_BID = 0.1
mean_FIELD = 4


poly_sigma = c(5.0, 0.5, 0.02, 0.02)
u_PEDIGREE_NAME = paste0("PE", str_pad(1:80, 2, pad = "0"))
u_GROWSEASON = paste0("2023:", str_pad(1:2, 2, pad = "0"))
u_REPETITION = 1:52
u_FIELD_NAME = paste0("FI", str_pad(1:3, 2, pad = "0"))

sigma_GROWSEASON = 1.0
simga_FIELD_NAME = 1.0
sigma_PLOT_BID = 1.0
mean_FIELD = 4
sigma_noise = 1.0

df_sim_list = simulate_data(
                poly_sigma,
                u_PEDIGREE_NAME,
                u_GROWSEASON,
                u_FIELD_NAME,
                u_REPETITION,
                sigma_GROWSEASON,
                simga_FIELD_NAME,
                sigma_PLOT_BID,
                mean_FIELD,
                sigma_noise
              )

df_sim = df_sim_list$df_sim
df_sim_reference = df_sim_list$df_sim_reference
df_sim_reference_rank = df_sim_list$df_sim_reference_rank

df_sim$TRAIT_VALUE %>% hist()

# Fit model
#

fitted_model = fit_lme_poly(df_sim %>% slice_sample(prop = 0.97), 3, 5:52)
fitted_model0 = fit_lme_poly(df_sim %>% slice_sample(prop = 0.97), 0, 5:52)
fitted_model1 = fit_lme_poly(df_sim %>% slice_sample(prop = 0.97), 1, 5:52)
fitted_model2 = fit_lme_poly(df_sim %>% slice_sample(prop = 0.97), 2, 5:52)

fitted_model0$df %>% filter(current_harvest_repetition >= 5) %>%
  select(full_harvest_repetition, current_harvest_repetition, skip, sigma_u, sigma_e) %>%
  group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
  summarise(h2=mean(sigma_u^2/(sigma_u^2+sigma_e^2)), .groups="drop")

fitted_model0$df %>% fitted_to_heritability() %>%
  ggplot() + aes(x = current_harvest_repetition, y = h2) + geom_line()

fitted_model0$df %>% fitted_to_heritability()





get_notnull_index = function(x) which(!sapply(x, is.null))


sapply(5:52, function(i) fitted_model$model[[i]]@theta)

str(VarCorr(fitted_model$model[[52]]))
as.data.frame(VarCorr(fitted_model$model[[52]]))

fitted_model$model[[52]]@theta

cbind(as.data.frame(VarCorr(fitted_model0$model[[52]]))[,1:3],
      sapply(5:52, function(i) as.data.frame(VarCorr(fitted_model0$model[[i]]))[,"vcov"]))
cbind(as.data.frame(VarCorr(fitted_model1$model[[52]]))[,1:3],
      sapply(5:52, function(i) as.data.frame(VarCorr(fitted_model1$model[[i]]))[,"vcov"]))
cbind(as.data.frame(VarCorr(fitted_model2$model[[52]]))[,1:3],
      sapply(5:52, function(i) as.data.frame(VarCorr(fitted_model2$model[[i]]))[,"vcov"]))
cbind(as.data.frame(VarCorr(fitted_model$model[[52]]))[,1:3],
      sapply(5:52, function(i) as.data.frame(VarCorr(fitted_model$model[[i]]))[,"vcov"]))

(VarCorr(fitted_model$model[[52]])$PEDIGREE_NAME)[1,1]

data.frame(apply(generate_poly(52, 3), 2, cumsum)) %>% mutate(REPETITION = 1:n()) %>% pivot_longer("poly0":"poly3", names_to = "poly", values_to = "cumsum") %>%
  ggplot() + aes(x = REPETITION, y = cumsum, color = poly) + geom_line()

Var(sqrt(52)*X) = 52*Var(X)
sd(sqrt(52)*X) = sqrt(52)*sd(X)


(max_r/max_repetition)*sqrt(max_repetition)

Var(mar_r/sqrt(52)*X) = mar_r^2/52*Var(X)
sd(mar_r/sqrt(52)*X) = mar_r/sqrt(52)*sd(X)



Var(mar_r/sqrt(52)*X) = mar_r^2/52*Var(X)


as.data.frame(ranef(fitted_model$model[[52]]))


(VarCorr(fitted_model2$model[[52]])$PEDIGREE_NAME)[1,1]

(52*(VarCorr(fitted_model$model[[52]])$PEDIGREE_NAME)[1,1])/(52*(VarCorr(fitted_model$model[[52]])$PEDIGREE_NAME)[1,1] + sigma(fitted_model$model[[52]])^2)



h2(52, fitted_model$model[[52]])

matrix(fitted_model$model[[5]]@theta, nrow =3, ncol=3)

fitted_model$model[[5]]
fitted_model$model[[5]]@theta



fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_precision(seq(0.1, 0.5, by = 0.1)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = precision) %>% print(width=Inf)



fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_precision(c(1.0)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = precision) %>% print(width=Inf)


fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_ndcg(seq(0.1, 0.5, by = 0.1)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = ndcg) %>% print(width=Inf)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_jaccard_similarity(seq(0.1, 0.5, by = 0.1)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = jaccard_similarity) %>% print(width=Inf)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_spearman(seq(0.1, 0.5, by = 0.1)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = spearman) %>% print(width=Inf)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_kendall(seq(0.1, 0.5, by = 0.1)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = kendall) %>% print(width=Inf)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_marginal_rank_change(seq(0.1, 0.5, by = 0.1)) %>%
  pivot_wider(names_from = current_harvest_repetition, values_from = marginal_rank_change) %>% print(width=Inf)

fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_marginal_rank_change(seq(0.1, 0.5, by = 0.1))

fitted_to_genetic_gain(fitted_model$df, seq(0.1, 0.5, by=0.1))



fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_genetic_correlation(seq(0.1, 0.5, by=0.1), benchmark_harvest_repetition = 25) %>%
  mutate(selection_intensity = as.character(selection_intensity)) %>%
  ggplot() + aes(x = current_harvest_repetition, y = genetic_correlation, color=selection_intensity) + geom_line()
  
fitted_model$df %>%
  fitted_to_rank(fitted_value_re) %>%
  rank_to_genetic_gain(seq(0.1, 0.5, by=0.1), benchmark_harvest_repetition = NULL) %>%
  mutate(selection_intensity = as.character(selection_intensity)) %>%
  ggplot() + aes(x = current_harvest_repetition, y = genetic_gain, color=selection_intensity) + geom_line()



