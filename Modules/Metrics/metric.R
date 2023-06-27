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
      pivot_wider(id_cols = c(PEDIGREE_NAME, full_harvest_repetition, skip, benchmark_harvest_repetition),
                  names_from = current_harvest_repetition,
                  values_from = rank) %>%
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

by_trait_metric = function(d, f, s_seq, ...) {
  d %>% group_by(OBSRVTN_REF_CD) %>% nest() %>%
    mutate(r = map(data, function(d){
      d %>% f(s_seq, ...)
    })) %>% select(-data) %>% unnest(r)
}

by_trait_metric_nos = function(d, f, ...) {
  d %>% group_by(OBSRVTN_REF_CD) %>% nest() %>%
    mutate(r = map(data, function(d){
      d %>% f(...)
    })) %>% select(-data) %>% unnest(r)
}

