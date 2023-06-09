# nDCG metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of nDCG.
#' 
#'
#'
rank_to_ndcg = function(df_rank, s_seq = c(1.0), ...) {
  rank_to_metric(df_rank, s_seq, "ndcg", metric_ndcg, ...)
}

metric_ndcg = function(selection_intensity, predict_rank, reference_rank){
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    mutate(scores = reference_rank <= selection_intensity * length(reference_rank)) %>%
    filter(predict_rank <= selection_intensity * length(predict_rank)) %>%
    arrange(predict_rank) %>%
    summarise(nDCG = sum(scores/log2(predict_rank + 1))/sum(1/log2(predict_rank + 1)) ) %>%
    pull(nDCG)
}

# Function to calculate DCG
calc_DCG <- function(scores) {
  DCG = sum(scores / log2(1:length(scores) + 1))
  return(DCG)
}
