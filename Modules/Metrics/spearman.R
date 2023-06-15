# Spearman Rank Correlation metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank A dataframe with the two columns: PEDIGREE_NAME, rank. The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Spearman Rank Correlation.
#' 
#' 
metric_spearman = function(selection_intensity, predict_rank, reference_rank){
  predict_rank %>%
    left_join(reference_rank, by = "PEDIGREE_NAME") %>%
    summarise(spearman = cor(rank.x, rank.y, method = "spearman") ) %>%
    pull(spearman)
}
