
# Multi harvest optimization metric calculation function. Given a class of lme4 class model this function will calculate
# Following metrics:
#     - Proportion of top x% pedigrees selected common with upto harvest i compared to full harvest
#     - Proportion of top x% pedigrees ranking changes with upto harvest i compared to full harvest
#     - Jaccard similarity beween two consequitive harvests
#     - Jaccard similarity beween upto harvest i to full harvest
#     - Genetic gain with each additional harvest
#     - Ranking correlations upto harvest i to full harvest


#' @param df Phenotype Data with single trait
#' @param selection_intensity Selection intensity for the downstream analysis. What % of total pedigrees one inteds to select. Should be bewteen 0 and 1
#' @param model lme4 class model 
#' 
#' @return Metrics: A data frame with the above metrics calculated
#' @return Top_Pedigrees : Top x% selected Pedigree
#' 
#' 
metrics_fun <- function(df, selection_intensity, model){
  trt = unique(df$OBSRVTN_REF_CD)
  # get the fitted values
  df_i <- df %>% mutate(Fitted_Values = fitted(model))
  # Find the maximum repetition
  max_repetition <- max(df_i$REPETITION)
  
  # Create a list to store top x% PEDIGREE_NAMEs for each repetition --------------------------------------------------
  top_x_percent_pedigrees <- vector("list", max_repetition)
  top_x_percent <- floor(n_distinct(df$PEDIGREE_NAME) * selection_intensity)
  
  # Loop through each repetition
  for (r in 1:max_repetition) {
    # Calculate the average fitted values for each PEDIGREE_NAME up to the current repetition
    temp_df <- df_i[df_i$REPETITION <= r, ]
    avg_fitted_values <- aggregate(Fitted_Values ~ PEDIGREE_NAME, data = temp_df, FUN = mean)
    
    # Rank the PEDIGREE_NAMEs based on the average fitted values and extract the top x% PEDIGREE_NAMEs
    top_pedigrees <- avg_fitted_values[order(avg_fitted_values$Fitted_Values, decreasing = TRUE), ]$PEDIGREE_NAME[1:top_x_percent]
    
    # Store the top x% PEDIGREE_NAMEs in the list
    top_x_percent_pedigrees[[r]] <- top_pedigrees
  }
  
  # Compare the intersection of the top x% PEDIGREE_NAMEs for each repetition with the full/max repetition
  intersection_counts <- sapply(top_x_percent_pedigrees, function(pedigrees) {
    sum(pedigrees %in% top_x_percent_pedigrees[[max_repetition]])
  })
  
  # Compute the proportion of top x% PEDIGREE_NAMEs that are in the full/max repetition
  proportion_in_final_top_x <- intersection_counts / floor(n_distinct(df_i$PEDIGREE_NAME) * selection_intensity)
  
  # Incorporate the ranking changes within the top x% -------------------------------------------------------------------
  ranking_changes <- numeric(max_repetition)
  
  for (r in 1:max_repetition) {
    current_top_pedigrees <- top_x_percent_pedigrees[[r]]
    final_top_pedigrees <- top_x_percent_pedigrees[[max_repetition]]
    
    # Calculate the ranking changes within the top x%
    # Calculate the ranking changes within the top x%
    ranking_change <- sum(abs(ifelse(is.na(match(current_top_pedigrees, final_top_pedigrees)),
                                     length(final_top_pedigrees) + 1,
                                     match(current_top_pedigrees, final_top_pedigrees)) - seq_along(current_top_pedigrees)))
    
    
    # Weight the ranking changes and store them in the vector
    ranking_changes[r] <- ranking_change / (length(current_top_pedigrees) * (length(current_top_pedigrees) - 1) / 2)
  }
  
  jaccard_similarity <- function(set1, set2) {
    intersection_size <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    
    if (union_size == 0) {
      return(1)
    }
    
    return(intersection_size / union_size)
  }
  
  jaccard_indices_consequtive <- jaccard_indices_overall <- c()
  for (i in 1:(max_repetition - 1)){
    jaccard_indices_consequtive[i] <- jaccard_similarity(top_x_percent_pedigrees[[i]], top_x_percent_pedigrees[[i+1]])
    jaccard_indices_overall[i] <- jaccard_similarity(top_x_percent_pedigrees[[i]], top_x_percent_pedigrees[[max_repetition]])
    
  }
  
  # Function to get top x% pedigrees based on fitted values
  get_top_pedigrees <- function(df, REPETITION) {
    df_filtered <- df[df$REPETITION <= REPETITION, ]
    unique_pedigrees <- unique(df_filtered$PEDIGREE_NAME)
    average_fitted_values <- sapply(unique_pedigrees, function(ped) {
      mean(df_filtered[df_filtered$PEDIGREE_NAME == ped, "Fitted_Values"])
    })
    top_n <- top_x_percent
    top_pedigrees_indices <- order(average_fitted_values, decreasing = TRUE)[1:top_n]
    top_pedigrees <- data.frame(Pedigree = unique_pedigrees[top_pedigrees_indices],
                                Fitted_Values = average_fitted_values[top_pedigrees_indices])
    return(top_pedigrees)
  }
  
  
  # Define the function for calculating genetic gain
  calculate_genetic_gain <- function(df, repetition) {
    top_pedigrees <- get_top_pedigrees(df, repetition)
    mean_top_pedigrees <- mean(top_pedigrees$Fitted_Values)
    mean_population <- mean(df$Fitted_Values)
    genetic_gain <- mean_top_pedigrees - mean_population
    return(genetic_gain)
  }
  
  # Calculate the genetic gain for each harvest
  genetic_gains <- sapply(1:max(df_i$REPETITION), function(r) {
    calculate_genetic_gain(df_i, repetition = r)
  })
  
  # Prepare an empty vector to store rank correlation coefficients
  rank_correlations <- vector("numeric", length = max_repetition - 1)
  
  # Calculate the rank correlation for each cumulative harvest
  for (harvest in 1:(max_repetition - 1)) {
    
    # Get the top x% unique pedigrees for each dataset
    top_current_pedigrees <- get_top_pedigrees(df, harvest)
    top_next_pedigrees <- get_top_pedigrees(df, harvest+1)
    
    # Calculate the rank correlation between the two sets of top x% unique pedigrees
    rank_correlations[harvest] <- cor(rank(top_current_pedigrees), rank(top_next_pedigrees), method = "spearman")
  }
  
  res <- data.frame(
    Harvest = c(1:max_repetition),
    Trait = rep(trt, max_repetition), 
    selection_intensity = rep(selection_intensity, max_repetition),
    Prop_final_top = proportion_in_final_top_x,
    Ranking_change = ranking_changes,
    Jaccard_Consequtive = c(NA, jaccard_indices_consequtive),
    Jaccard_Overall = c(NA, jaccard_indices_overall),
    Genetic_Gain = genetic_gains,
    rank_correlations = c(NA,rank_correlations)
  )
   list(Metrics = res, 
       Top_Pedigrees = top_x_percent_pedigrees)
}