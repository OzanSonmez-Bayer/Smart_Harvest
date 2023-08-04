# Marginal Rank Change metric
#

#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1.
#' @param predict_rank The rank estimation for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' @param reference_rank The rank reference for each PEDIGREE_NAME. The PEDIGREE_NAME with smaller rank (i.e. 1, 2, etc) are selected.
#' 
#' @return the numeric of Marginal Rank Change metric.

rank_to_marginal_rank_change = function(df_rank, s_seq, ...) {
  rank_to_metric(df_rank, s_seq, "marginal_rank_change", metric_marginal_rank_change, ...)
}

# Marginal Rank Change Metric Function
#
#' This function computes the Marginal Rank Change between predicted and reference ranks.
#' It quantifies the difference in ranks and normalizes it by the total number of possible changes.
#'
#' @param selection_intensity A numeric value representing the selection intensity.
#' @param predict_rank A numeric vector representing the predicted ranks.
#' @param reference_rank A numeric vector representing the reference ranks.
#' @return A numeric value representing the Marginal Rank Change.
#' 
metric_marginal_rank_change = function(selection_intensity, predict_rank, reference_rank) {
  data.frame(predict_rank = predict_rank, reference_rank = reference_rank) %>%
    mutate(fake_reference_rank = pmin(reference_rank, as.integer(selection_intensity*length(reference_rank) + 1))) %>%
    mutate(abs_rc = abs(predict_rank - fake_reference_rank)) %>%
    filter(predict_rank <= selection_intensity*length(predict_rank)) %>%
    summarise(rc = 1 - sum(abs_rc)/(length(abs_rc) * (length(abs_rc) - 1) / 2)) %>%
    pull(rc)
}
