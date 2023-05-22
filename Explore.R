library(tidyverse)
library(plotly)
library(lme4)
source("/repos/Smart_Harvest/Model_Metrics.R")
select <- dplyr::select

# Input Data
load("/mnt/Data/Pepper_data.RData")
data0 <- df_pepper


data0 <- data0 %>% mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE))
model_formula <- formula("TRAIT_VALUE ~  (1 | FIELD_NAME) + (1 | PEDIGREE_NAME) + (1 | REPETITION)")
trt_list <- unique(data0$OBSRVTN_REF_CD)

metric_all <- data.frame()
sel_list <- list()
for (trt in trt_list){
  df <- data0 %>% filter(OBSRVTN_REF_CD == trt)
  model <- lmer(model_formula, data = df)
  metric_df_trt <- data.frame()
  sel_list[[trt]] = list()
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    metric_df_i <- metrics_fun(df, sel, model)
    metric_df_trt <- rbind(metric_df_trt, metric_df_i$Metrics)
    sel_list[[trt]][[as.character(sel)]] = metric_df_i$Top_Pedigrees
    print(paste0("Trait: ", trt, "---- Selection Intensity: ",sel))
    
  }
  metric_all <- rbind(metric_all, metric_df_trt)
}

ggplotly(
 metric_all %>%
  filter(selection_intensity == 0.1) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)

ggplotly(
 metric_all %>%
  filter(selection_intensity == 0.5) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)


# Assume that we skip every other harvest:

data1 <- data0 %>%
  filter(REPETITION %% 2 == 1) %>%
  mutate(REPETITION = (REPETITION + 1)/2)

metric_all2 <- data.frame()
sel_list2 <- list()
for (trt in trt_list){
  df <- data1 %>% filter(OBSRVTN_REF_CD == trt)
  model <- lmer(model_formula, data = df)
  metric_df_trt <- data.frame()
  sel_list2[[trt]] = list()
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    metric_df_i <- metrics_fun(df, sel, model)
    metric_df_trt <- rbind(metric_df_trt, metric_df_i$Metrics)
    sel_list2[[trt]][[as.character(sel)]] = metric_df_i$Top_Pedigrees
    print(paste0("Trait: ", trt, "---- Selection Intensity: ",sel))
    
  }
  metric_all2 <- rbind(metric_all2, metric_df_trt)
}


ggplotly(
 metric_all2 %>%
  filter(selection_intensity == 0.1) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)

jaccard_similarity <- function(set1, set2) {
    intersection_size <- length(intersect(set1, set2))
    union_size <- length(union(set1, set2))
    
    if (union_size == 0) {
      return(1)
    }
    
    return(intersection_size / union_size)
  }


res_comp <- data.frame()
for (trt in trt_list){
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    full_harvest <- sel_list[[trt]][[as.character(sel)]]
    skip_harvest = sel_list2[[trt]][[as.character(sel)]]
    comp_df <- data.frame()
    for (i in 1:length(skip_harvest)){
      set1 <- skip_harvest[[i]]
      set2 <- full_harvest[[(2*i)-1]]
      proportion_common <- n_distinct(intersect(set1, set2))/n_distinct(set1)
      
      ranking_change <- sum(abs(ifelse(is.na(match(set1, set2)),
                                     length(set2) + 1,
                                     match(set1, set2)) - seq_along(set1)))
      ranking_changes <- ranking_change / (length(set1) * (length(set2) - 1) / 2)
      jaccard_sim <- jaccard_similarity(set1, set2)
      res <- data.frame(
        Trait = trt,
        selection_intensity = sel,
        Harvest = i,
        proportion_common = proportion_common,
        ranking_changes = ranking_changes,
        jaccard_sim = jaccard_sim
      )
      comp_df <- rbind(comp_df, res)
    }
    res_comp <- rbind(res_comp, comp_df)
  }
}

ggplotly(
 res_comp %>%
  filter(selection_intensity == 0.1) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)

## lets skip every other starting from 1, so we will skip 1, 3, 5, ...
# Assume that we skip every other harvest:

data2 <- data0 %>%
  filter(REPETITION %% 2 == 0) %>%
  mutate(REPETITION = (REPETITION )/2)

metric_all3 <- data.frame()
sel_list3 <- list()
for (trt in trt_list){
  df <- data2 %>% filter(OBSRVTN_REF_CD == trt)
  model <- lmer(model_formula, data = df)
  metric_df_trt <- data.frame()
  sel_list3[[trt]] = list()
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    metric_df_i <- metrics_fun(df, sel, model)
    metric_df_trt <- rbind(metric_df_trt, metric_df_i$Metrics)
    sel_list3[[trt]][[as.character(sel)]] = metric_df_i$Top_Pedigrees
    print(paste0("Trait: ", trt, "---- Selection Intensity: ",sel))
    
  }
  metric_all3 <- rbind(metric_all3, metric_df_trt)
}


ggplotly(
 metric_all3 %>%
  filter(selection_intensity == 0.1) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)

res_comp2 <- data.frame()
for (trt in trt_list){
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    full_harvest <- sel_list[[trt]][[as.character(sel)]]
    skip_harvest = sel_list3[[trt]][[as.character(sel)]]
    comp_df <- data.frame()
    for (i in 1:length(skip_harvest)){
      set1 <- skip_harvest[[i]]
      set2 <- full_harvest[[(2*i)]]
      proportion_common <- n_distinct(intersect(set1, set2))/n_distinct(set1)
      
      ranking_change <- sum(abs(ifelse(is.na(match(set1, set2)),
                                     length(set2) + 1,
                                     match(set1, set2)) - seq_along(set1)))
      ranking_changes <- ranking_change / (length(set1) * (length(set2) - 1) / 2)
      jaccard_sim <- jaccard_similarity(set1, set2)
      res <- data.frame(
        Trait = trt,
        selection_intensity = sel,
        Harvest = i,
        proportion_common = proportion_common,
        ranking_changes = ranking_changes,
        jaccard_sim = jaccard_sim
      )
      comp_df <- rbind(comp_df, res)
    }
    res_comp2 <- rbind(res_comp, comp_df)
  }
}

ggplotly(
 res_comp2 %>%
  filter(selection_intensity == 0.1) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)


## Lets say we take every 5 th harvest
# Assume that we skip every other harvest:

data5 <- data0 %>%
  filter(REPETITION %% 5 == 0) %>%
  mutate(REPETITION = (REPETITION )/5)

metric_all5 <- data.frame()
sel_list5 <- list()
for (trt in trt_list){
  df <- data5 %>% filter(OBSRVTN_REF_CD == trt)
  model <- lmer(model_formula, data = df)
  metric_df_trt <- data.frame()
  sel_list5[[trt]] = list()
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    metric_df_i <- metrics_fun(df, sel, model)
    metric_df_trt <- rbind(metric_df_trt, metric_df_i$Metrics)
    sel_list5[[trt]][[as.character(sel)]] = metric_df_i$Top_Pedigrees
    print(paste0("Trait: ", trt, "---- Selection Intensity: ",sel))
    
  }
  metric_all5 <- rbind(metric_all5, metric_df_trt)
}


ggplotly(
 metric_all5 %>%
  filter(selection_intensity == 0.1) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric, scales = "free", ncol = 3) +
  theme_bw()
)

res_comp5 <- data.frame()
for (trt in trt_list){
  for (sel in c(0.1, 0.2, 0.3, 0.5)){
    full_harvest <- sel_list[[trt]][[as.character(sel)]]
    skip_harvest = sel_list5[[trt]][[as.character(sel)]]
    comp_df <- data.frame()
    for (i in 1:length(skip_harvest)){
      set1 <- skip_harvest[[i]]
      set2 <- full_harvest[[(5*i)]]
      proportion_common <- n_distinct(intersect(set1, set2))/n_distinct(set1)
      
      ranking_change <- sum(abs(ifelse(is.na(match(set1, set2)),
                                     length(set2) + 1,
                                     match(set1, set2)) - seq_along(set1)))
      ranking_changes <- ranking_change / (length(set1) * (length(set2) - 1) / 2)
      jaccard_sim <- jaccard_similarity(set1, set2)
      res <- data.frame(
        Trait = trt,
        selection_intensity = sel,
        Harvest = i,
        proportion_common = proportion_common,
        ranking_changes = ranking_changes,
        jaccard_sim = jaccard_sim
      )
      comp_df <- rbind(comp_df, res)
    }
    res_comp5 <- rbind(res_comp, comp_df)
  }
}

ggplotly(
 res_comp5 %>%
  filter(selection_intensity == 0.3) %>% 
  gather("Metric", "Value", -c(Harvest, Trait, selection_intensity)) %>%
  ggplot(., aes(x = Harvest, y = Value, col = Trait)) +
  geom_line() + 
  facet_wrap( ~ Metric,  ncol = 3) +
  theme_bw()
)


