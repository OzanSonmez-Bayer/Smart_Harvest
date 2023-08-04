# Precision Metric
#
#' This function calculates the Precision metric, a measure of the number of true positive results divided by the number of true positives and false positives.
#' It is a useful measure in classification tasks.
#' The function takes a data frame of ranks and a vector of selection intensities and returns the Precision.
#'
#' @param df_rank A data frame containing the ranks.
#' @param s_seq A vector of selection intensities, defaulting to c(1.0). Represents the percentage of total pedigrees one intends to select. Should be between 0 and 1.
#' @param ... Additional arguments to pass to the rank_to_metric function.
#'
#' @return A data frame containing the Precision metric.
#' 
rank_to_precision = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "precision", metric_precision, ...)
}

metric_precision = function(selection_intensity, predict_rank, reference_rank){
  predict_selected = predict_rank <= selection_intensity*length(predict_rank)
  reference_selected = reference_rank <= selection_intensity*length(reference_rank)
  return(sum(predict_selected * reference_selected, na.rm = T)/sum(predict_selected, na.rm = T))
}
