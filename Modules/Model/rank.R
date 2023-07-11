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

by_trait_fitted_to_rank = function(d, fitted_col, rank_signs) {
  d %>%
    group_by(OBSRVTN_REF_CD) %>%
    nest() %>%
    mutate(r = map(data, function(d) {
      d %>% fitted_to_rank(fitted_col, rank_signs[OBSRVTN_REF_CD])
    })) %>% select(-data) %>% unnest(r)
}

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

by_trait_fitted_to_expected_rank = function(d, rank_signs, ...) {
  d %>%
    group_by(OBSRVTN_REF_CD) %>%
    nest() %>%
    mutate(r = map(data, function(d) {
      d %>% fitted_to_expected_rank(rank_signs[OBSRVTN_REF_CD], ...)
    })) %>% select(-data) %>% unnest(r)
}


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

get_expectedRank_simple = function(mu, sigma, rank_sign = 1) {
  # If rank_sign = 1 then larger values has larger rank
  # If rank_sign = -1 then smaller values has larger rank
  
  expected_rank = apply(p_matrix(mu, sigma), (-rank_sign + 3) %/%2, sum)
  return(expected_rank)
}


