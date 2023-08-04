#' A general function to compute a metric. It implements the common pattern of rank and selection metric,
#' given rank and selection intensity. The function takes a data frame of ranks, a sequence of selection
#' intensities, a metric name, a metric function, and an optional benchmark harvest repetition.
#' It returns a data frame with the computed metric.
#'
#' @param df_rank A data frame containing the ranks.
#' @param s_seq A sequence of selection intensities.
#' @param metric_name A string representing the name of the metric.
#' @param metric_fun A function to compute the metric.
#' @param benchmark_harvest_repetition An optional integer representing the benchmark harvest repetition.
#' @return A data frame with the computed metric.
#' 
rank_to_metric = function(df_rank, s_seq, metric_name, metric_fun, benchmark_harvest_repetition = NULL) {
  # If benchmark_harvest_repetition is NULL, use full_harvest_repetition as the benchmark
  if(is.null(benchmark_harvest_repetition)) {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = full_harvest_repetition)
  } else {
    df_rank = df_rank %>%
      mutate(benchmark_harvest_repetition = benchmark_harvest_repetition)
  }
  
  # Compute the metric for each selection intensity and combine the results
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

#' A general function to compute a metric for each trait in the input data frame and selection intensity.
#' It groups the data by observation reference code and applies the given function to compute the metric.
#'
#' @param d A data frame containing the input data.
#' @param f A function to compute the metric.
#' @param s_seq A sequence of selection intensities.
#' @param ... Additional arguments to pass to the function f.
#' @return A data frame with the computed metric for each trait.
#'
by_trait_metric = function(d, f, s_seq, ...) {
  # Group by observation reference code, nest the data, and apply the function f
  d %>% group_by(OBSRVTN_REF_CD) %>% nest() %>%
    mutate(r = map(data, function(d){
      d %>% f(s_seq, ...)
    })) %>% select(-data) %>% unnest(r)
}

#' A general function to compute a metric for each trait in the input data frame. No selection intensity version for those metrics that do not use selection intensity.
#' It groups the data by observation reference code and applies the given function to compute the metric.
#'
#' @param d A data frame containing the input data.
#' @param f A function to compute the metric.
#' @param ... Additional arguments to pass to the function f.
#' @return A data frame with the computed metric for each trait.
#'
by_trait_metric_nos = function(d, f, ...) {
  # Group by observation reference code, nest the data, and apply the function f
  d %>% group_by(OBSRVTN_REF_CD) %>% nest() %>%
    mutate(r = map(data, function(d){
      d %>% f(...)
    })) %>% select(-data) %>% unnest(r)
}
