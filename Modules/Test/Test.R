# test the workflow for a given data set

library(tidyverse)
library(stringr)
library(lme4)

setwd("/repos/Smart_Harvest/Modules/Test")
source("../Model/Model.R")
source("../Model/rank.R")
source("../Metrics/precision.R")
source("../Metrics/ndcg.R")
source("../Metrics/jaccard_similarity.R")
source("../Metrics/spearman.R")

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


# Preparing simulated data set
#
set.seed(1234)

poly_max_degree = 3
poly_sigma = c(5.0, 0.5, 0.01, 0.01)

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

v_PEDIGREE_NAME = generate_poly_coef(u_PEDIGREE_NAME, poly_max_degree, poly_sigma)
f_PEDIGREE_NAME = generate_poly_fitted_value(max(u_REPETITION), poly_max_degree, v_PEDIGREE_NAME)

df_sim_reference = df_sim_design %>%
  left_join(f_PEDIGREE_NAME %>% rename(FITTED_VALUE_PEDIGREE_NAME = FITTED_VALUE), by = c("PEDIGREE_NAME", "REPETITION")) %>%
  left_join(v_GROWSEASON %>% enframe(name = "GROWSEASON", value = "FITTED_VALUE_GROWSEASON"), by = "GROWSEASON") %>%
  left_join(v_FIELD_NAME %>% enframe(name = "FIELD_NAME", value = "FITTED_VALUE_FIELD_NAME"), by = "FIELD_NAME") %>% 
  left_join(v_PLOT_BID %>% enframe(name = "PLOT_BID", value = "FITTED_VALUE_PLOT_BID"), by = "PLOT_BID") %>%
  mutate(FITTED_VALUE_TRAIT_VALUE = FITTED_VALUE_PEDIGREE_NAME + FITTED_VALUE_GROWSEASON + FITTED_VALUE_FIELD_NAME + FITTED_VALUE_PLOT_BID) %>%
  mutate(TRAIT_VALUE = pmin(pmax(round(FITTED_VALUE_TRAIT_VALUE + rnorm(n(), sd = 0.1), digits = 0), 0), 9) )

df_sim_reference_rank = df_sim_reference %>%
  group_by(PEDIGREE_NAME, REPETITION) %>%
  summarise(mean_re = mean(FITTED_VALUE_PEDIGREE_NAME), .groups = "drop") %>%
  group_by(PEDIGREE_NAME) %>%
  summarise(mean_re = mean(mean_re), .groups = "drop") %>%
  mutate(rank = rank(mean_re))

df_sim = df_sim_reference %>%
  select(TRAIT_VALUE, PEDIGREE_NAME, GROWSEASON, FIELD_NAME, PLOT_BID, REPETITION)


# Fit model
#

fitted_model = fit_lme_poly(df_sim %>% filter(REPETITION <= 20) %>% slice_sample(prop = 0.95),
                            poly_max_degree)

metric_precision(0.2, rank_random_effect(fitted_model$df), df_sim_reference_rank)
metric_ndcg(0.2, rank_random_effect(fitted_model$df), df_sim_reference_rank)
metric_jaccard_similarity(0.2, rank_random_effect(fitted_model$df), df_sim_reference_rank)
metric_spearman(0.2, rank_random_effect(fitted_model$df), df_sim_reference_rank)

