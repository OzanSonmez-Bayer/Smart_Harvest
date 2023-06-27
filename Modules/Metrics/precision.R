# Precision Metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of precision.
#' 
#' 
#' 
rank_to_precision = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "precision", metric_precision, ...)
}

metric_precision = function(selection_intensity, predict_rank, reference_rank){
  predict_selected = predict_rank <= selection_intensity*length(predict_rank)
  reference_selected = reference_rank <= selection_intensity*length(reference_rank)
  return(sum(predict_selected * reference_selected, na.rm = T)/sum(predict_selected, na.rm = T))
}
