#' Rank Fitted Values
#'
#' This function ranks the fitted values based on the specified column and sign.
#'
#' @param df A data frame containing the fitted values.
#' @param fitted_col The column name containing the fitted values to rank.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @return A data frame containing the ranked values.
#' 
fitted_to_rank = function(df, fitted_col, rank_sign = 1) {
  df_output = df %>%
    group_by(PEDIGREE_NAME, full_harvest_repetition, skip, current_harvest_repetition) %>%
    summarise(mean_value = mean(!! sym(fitted_col)), .groups = "drop") %>%
    group_by(full_harvest_repetition, current_harvest_repetition, skip) %>%
    mutate(rank = rank(rank_sign*mean_value)) %>%
    mutate(percentile_rank = round(100 * (rank - 0.5)/length(rank))) %>%
    ungroup() %>%
    arrange(current_harvest_repetition) %>%
    select(PEDIGREE_NAME, full_harvest_repetition, skip, current_harvest_repetition, mean_value, rank, percentile_rank)
  return(df_output)
}

#' Rank Fitted Values by Trait
#'
#' This function ranks the fitted values by trait based on the specified column and sign.
#'
#' @param d A data frame containing the fitted values.
#' @param fitted_col The column name containing the fitted values to rank.
#' @param rank_signs A vector of numeric values indicating the direction of ranking for each trait.
#' @return A data frame containing the ranked values by trait.
#' 
by_trait_fitted_to_rank = function(d, fitted_col, rank_signs) {
  d %>%
    group_by(OBSRVTN_REF_CD) %>%
    nest() %>%
    mutate(r = map(data, function(d) {
      d %>% fitted_to_rank(fitted_col, rank_signs[OBSRVTN_REF_CD])
    })) %>% select(-data) %>% unnest(r)
}

#' Compute Expected Rank from Fitted Values
#'
#' This function computes the expected rank for the given fitted values using the empirical Bayes ranking methods as described by Laird NM and Louis TA.
#'
#' @param df A data frame containing the fitted values.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @param mc.cores The number of cores to use for parallel processing (default is NULL).
#' @return A data frame containing the expected rank, rank on expected rank, percentile of expected rank, standard deviation of expected rank, original rank, and percentile of original rank.
#' @references Laird NM and Louis TA. Empirical Bayes Ranking Methods.
#' 
fitted_to_expected_rank = function(df, rank_sign = 1, mc.cores = NULL) {
  df_output = df %>%
    select(PEDIGREE_NAME, full_harvest_repetition, skip, current_harvest_repetition,
           condval_u,
           condsd_u) %>%
    group_by(PEDIGREE_NAME, full_harvest_repetition, skip, current_harvest_repetition) %>%
    summarise(condval_u = mean(rank_sign*condval_u), condsd_u = mean(condsd_u), .groups = "drop") %>%
    ungroup() %>%
    group_by(full_harvest_repetition, skip, current_harvest_repetition) %>%
    nest() %>%
    mutate(rank_df = map(data, function(d) get_expectedRank(d$condval_u, d$condsd_u, mc.cores = mc.cores))) %>%
    unnest(cols = c(data, rank_df)) %>%
    ungroup() %>%
    select(PEDIGREE_NAME, full_harvest_repetition, skip, current_harvest_repetition, 
           condval_u,
           condsd_u,
           rank,
           expected_rank,
           rank_on_er,
           percentile_er,
           sd_er)
  return(df_output)
}

#' Compute Expected Rank for Each Trait
#'
#' This function computes the expected rank for each trait in the given data frame using the empirical Bayes ranking methods.
#'
#' @param d A data frame containing the fitted values for different traits.
#' @param rank_signs A named vector of numeric values indicating the direction of ranking for each trait (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @param ... Additional arguments to be passed to the `fitted_to_expected_rank` function.
#' @return A data frame containing the expected rank, rank on expected rank, percentile of expected rank, standard deviation of expected rank, original rank, and percentile of original rank for each trait.
#' @references Laird NM and Louis TA. Empirical Bayes Ranking Methods.
#' 
by_trait_fitted_to_expected_rank = function(d, rank_signs, ...) {
  d %>%
    group_by(OBSRVTN_REF_CD) %>%
    nest() %>%
    mutate(r = map(data, function(d) {
      d %>% fitted_to_expected_rank(rank_signs[OBSRVTN_REF_CD], ...)
    })) %>% select(-data) %>% unnest(r)
}

#' Trace Rank
#'
#' This function traces the rank of the selected pedigree names based on the specified selection intensity and benchmark harvest repetition.
#'
#' @param df_rank A data frame containing the rank information.
#' @param selection_intensity A numeric value indicating the selection intensity (default is 1.0).
#' @param benchmark_harvest_repetition The benchmark harvest repetition (default is NULL).
#' @param which_rank The column name containing the rank values to trace (default is NULL).
#' @return A data frame containing the traced rank.
#' 
rank_trace = function(df_rank, selection_intensity = c(1.0), benchmark_harvest_repetition=NULL, which_rank = NULL) {
  if(is.null(benchmark_harvest_repetition)) {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = full_harvest_repetition)
  } else {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = benchmark_harvest_repetition)
  }
  
  if(is.null(which_rank)) {which_rank = "rank"}
  
  selected_PEDIGREE_NAME = df_rank %>%
    filter(current_harvest_repetition == benchmark_harvest_repetition) %>%
    filter(get(which_rank) <= selection_intensity*length(get(which_rank)) ) %>%
    distinct(PEDIGREE_NAME) %>%
    pull(PEDIGREE_NAME)
  
  df_rank %>%
    filter(PEDIGREE_NAME %in% selected_PEDIGREE_NAME) %>%
    return()
}

#' Rank Random Effect
#'
#' This function ranks the random effect based on the specified sign.
#'
#' @param df A data frame containing the random effect values.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @return A data frame containing the ranked random effect.
#' 
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

#' Rank Specified Column
#'
#' This function ranks the values in the specified column based on the specified sign.
#'
#' @param df A data frame containing the values to rank.
#' @param col The column name containing the values to rank.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @return A data frame containing the ranked values.
#' 
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


# Additional functions for computing probabilities and expected ranks
# These functions are used internally to support the ranking process.

p_fun = function(mu0, sigma0, mu1, sigma1) ifelse(sigma0^2+sigma1^2 > 0, pnorm((mu1-mu0)/sqrt(sigma0^2+sigma1^2)), 0)
p_matrix = function(mu, sigma) {
  n = length(mu)
  p = matrix(NA, nrow=n, ncol=n)
  for(k in 1:n){
    for(j in 1:n) {
      if(k != j) {
        p[j,k] = p_fun(mu[k], sigma[k], mu[j], sigma[j])
      } else {
        p[j,k] = 1
      }
    }
  }
  return(p)
}

c_fun = function(mu0, sigma0, mu1, sigma1, mu2, sigma2) {
  ifelse(sigma0^2+sigma1^2+sigma2^2 > 0,
  integrate(function(a) {
    pnorm(a, mean=mu1, sd=sigma1, lower.tail=F)*pnorm(a, mean=mu2, sd=sigma2, lower.tail=F)*dnorm(a, mean=mu0, sd=sigma0)
  },
  min(mu0, mu1, mu2)-4*sqrt(sigma0^2+sigma1^2+sigma2^2),
  max(mu0, mu1, mu2)+4*sqrt(sigma0^2+sigma1^2+sigma2^2)
  )$value, 0)
}

c_help_fun = function(k, n, mu, sigma) {
  cp = array(NA, dim = c(n,n))
  for(j in 1:n) {
    for(l in 1:n) {
      if(!is.na(cp[l,j])){
        cp[j,l] = cp[l,j]
      } else {
        if( j == k && l == k) {
          cp[j,l] = 1
        } else {
          if (j == l) {
            cp[j,l] = p_fun(mu[k], sigma[k], mu[j], sigma[j])
          } else { if(j == k) {
            cp[j,l] = p_fun(mu[k], sigma[k], mu[l], sigma[l])
          } else { if(l == k) {
            cp[j,l] = p_fun(mu[k], sigma[k], mu[j], sigma[j])
          } else {
            cp[j,l] = c_fun(mu[k], sigma[k], mu[j], sigma[j], mu[l], sigma[l])
          }
          }
          }
        }
      }
    }
  }
  return(cp)
}

c_matrix = function(mu, sigma, mc.cores = NULL) {
  n = length(mu)
  if(is.null(mc.cores)) {
    out = lapply(1:n, function(k) c_help_fun(k, n, mu, sigma))
  } else {
    out = mclapply(1:n, function(k) c_help_fun(k, n, mu, sigma), mc.cores = mc.cores)
  }
  return(array(unlist(out), dim=c(n, n, n)))
}

v_fun = function(mu, sigma, mc.cores = NULL) {
  n = length(mu)
  p = p_matrix(mu, sigma)
  cp = c_matrix(mu, sigma, mc.cores = mc.cores)
  out = sapply(1:n, function(k) sum(cp[,,k] -(p[,k] %o% p[,k])))
  return(out)
}

cv_fun = function(mu, sigma) {
  n = length(mu)
  p = p_matrix(mu, sigma)
  cp_star = c_matrix(-mu, sigma)
  Reduce("+", lapply(1:n, function(j) cp_star[,,j] - (p[j,] %o% p[j,])))
}

#' Get Expected Rank
#'
#' This function computes the expected rank using the empirical Bayes ranking methods as described by Laird NM and Louis TA.
#' It ranks the values based on the specified sign.
#'
#' @param mu A vector of mean values.
#' @param sigma A vector of standard deviations.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @param mc.cores The number of cores to use for parallel processing (default is NULL).
#' @return A data frame containing the expected rank, rank on expected rank, percentile of expected rank, standard deviation of expected rank, original rank, and percentile of original rank.
#' @references Laird NM and Louis TA. Empirical Bayes Ranking Methods.
#' 
get_expectedRank = function(mu, sigma, rank_sign = 1, mc.cores = NULL) {
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  
  expected_rank = apply(p_matrix(mu, sigma), (-rank_sign + 3) %/%2, sum)
  rank_on_er = rank(expected_rank)
  sd_er = sqrt(v_fun(mu, sigma, mc.cores = mc.cores))
  percentile_er <- round(100 * (expected_rank - 0.5)/length(mu))
  rank_org = rank(rank_sign*mu)
  percentile_rank_org = round(100 * (rank_org - 0.5)/length(rank_org))
  return(data.frame(expected_rank = expected_rank,
                    rank_on_er = rank_on_er,
                    percentile_er = percentile_er,
                    sd_er = sd_er,
                    rank = rank_org,
                    percentile_rank = percentile_rank_org
                    ))
}

#' Get Expected Rank with Covariance
#'
#' This function computes the expected rank along with the variance and covariance using the empirical Bayes ranking methods.
#'
#' @param mu A vector of mean values.
#' @param sigma A vector of standard deviations.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @return A list containing the expected rank (er), variance of expected rank (var_er), and covariance of expected rank (cov_er).
#' @references Laird NM and Louis TA. Empirical Bayes Ranking Methods.
#' 
get_expectedRank_cov = function(mu, sigma, rank_sign = 1) {
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  p = apply(p_matrix(mu, sigma), (-rank_sign + 3) %/%2, sum)
  v = v_fun(mu, sigma)
  cv = cv_fun(mu, sigma)
  for(i in 1:length(v)) {
    cv[i, i] = v[i]
  }
  return(list(er = p, var_er = v, cov_er = cv))
}

#' Get Expected Rank (Simple Version)
#'
#' This function computes the expected rank without computing the standard error and covariance to reduce the running time.
#'
#' @param mu A vector of mean values.
#' @param sigma A vector of standard deviations.
#' @param rank_sign A numeric value indicating the direction of ranking (1 for larger values having larger ranks, -1 for smaller values having larger ranks).
#' @return A vector containing the expected rank.
#' @references Laird NM and Louis TA. Empirical Bayes Ranking Methods.
get_expectedRank_simple = function(mu, sigma, rank_sign = 1) {
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  
  expected_rank = apply(p_matrix(mu, sigma), (-rank_sign + 3) %/%2, sum)
  return(expected_rank)
}


