# Spearman Rank Correlation metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Spearman Rank Correlation.
rank_to_spearman = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "spearman", metric_spearman, ...)
}

metric_spearman = function(selection_intensity, predict_rank, reference_rank){
  if(sd(predict_rank) < 1e-12 || sd(reference_rank) < 1e-12) return(NA)
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    filter(predict_rank <= selection_intensity*length(predict_rank) | reference_rank <= selection_intensity*length(reference_rank)) %>%
    summarise(spearman = cor(predict_rank, reference_rank, method = "spearman") ) %>%
    pull(spearman)
}
