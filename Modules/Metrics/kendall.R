# Kendall Rank Correlation metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Kendall Rank Correlation.
rank_to_kendall = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "kendall", metric_kendall, ...)
}

metric_kendall = function(selection_intensity, predict_rank, reference_rank){
  if(sd(predict_rank) < 1e-12 || sd(reference_rank) < 1e-12) return(NA)
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    filter(predict_rank <= selection_intensity*length(predict_rank) | reference_rank <= selection_intensity*length(reference_rank)) %>%
    summarise(kendall = cor(predict_rank, reference_rank, method = "kendall")) %>%
    pull(kendall)
}
