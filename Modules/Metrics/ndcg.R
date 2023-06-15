# nDCG metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of nDCG.
#' 
#' 
metric_ndcg = function(selection_intensity, predict_rank, reference_rank){
  predict_rank %>% filter(rank <= selection_intensity * length(rank)) %>%
    arrange(rank) %>%
    left_join(reference_rank %>% mutate(scores = rank <= selection_intensity * length(rank)) %>% select(PEDIGREE_NAME, scores), by = "PEDIGREE_NAME") %>%
    summarise(nDCG = sum(scores/log2(rank + 1))/sum(1/log2(rank + 1)) ) %>%
    pull(nDCG)
}

# Function to calculate DCG
calc_DCG <- function(scores) {
  DCG = sum(scores / log2(1:length(scores) + 1))
  return(DCG)
}
