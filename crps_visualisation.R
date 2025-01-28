## clear environment
rm(list = ls())

## set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))

## load packages
library(rstudioapi)
library(arrow)
library(dplyr)
library(scoringRules)
library(tidyverse)
library(ggplot2)

# read-in predictive samples
# renv::snapshot()
# renv::update() 
# renv::restore()

## -----
# read-in predictive samples from benchmark models and submitted forecasts
## -----
benchmark_names <- c("boot_240", "conflictology", "last", "zero")
submissions_names <- c("bodentien_rueter_negbin", "bodentien_rueter_neuralnet", "conflictforecast_v2", "Neg_Bin_GAM", "Neg_Bin_GLMM",
                       "P_GAM", "P_GLMM", "quantile_forecast", "ShapeFinder", "submission_final_gpcmm", "submission_final_hpmm",
                       "submission_final_omm", "submission_muchlinski_thornhill", "tft", "TW_GAM", "TW_GLMM", "unito_transformer")

model_names <- c(benchmark_names, submissions_names)
n_models <- length(model_names)

appendix_names <- paste0("_cm_Y20", c(18, 19, 20, 21, 22, 23, 24))
predictive_samples <- list()

for (m in 1:length(model_names)) {
  model_files <- paste0("../Data/predictions/", model_names[m], "/", model_names[m], appendix_names, ".parquet")
  predictive_samples[[m]] <- lapply(model_files, arrow::read_parquet)
  predictive_samples[[m]] <- bind_rows(predictive_samples[[m]])
}
names(predictive_samples) <- model_names


## -----
# read-in actual observations
## -----
observations_files <- paste0("../Data/observations/cm_actuals_20",  c(18, 19, 20, 21, 22, 23, 24), ".parquet")
observations <- lapply(observations_files, arrow::read_parquet)
# names(observations) <- paste0("Y20", c(18, 19, 20, 21, 22, 23, 24))
# create one dataframe from the list
observations <- bind_rows(observations)


## -----
# compute crps on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----
country_ids <- unique(observations$country_id)
actuals_ids <- 457:528 # month_ids for 01-2018, 02-2018, ..., 12-2023
cm_pairs <- cbind(rep(country_ids, each = length(actuals_ids)),
                  rep(actuals_ids, length(country_ids)))

models_crps <- list()

# for (m in 1:n_models) {
#   crps_m <- apply(cm_pairs, 1, function(cm_pair) {
#     print(paste0("Benchmark/model (", m, "/", n_models, "): ", model_names[m], ", country: ", cm_pair[1], ", month: ", cm_pair[2]))
#     true_observation <- observations %>%
#       filter(country_id == cm_pair[1] & month_id == cm_pair[2]) %>%
#       select(outcome)
#     pred_sample <- predictive_samples[[m]] %>%
#       filter(country_id == cm_pair[1] & month_id == cm_pair[2]) %>%
#       select(outcome)
#     crps_sample(y = unlist(true_observation),
#                 dat = unlist(pred_sample))
#   })
#   models_crps[[m]] <- data.frame("country_id" = cm_pairs[,1], 
#                                  "month_id" = cm_pairs[,2],
#                                  "crps" = crps_m)
# }
# names(models_crps) <- model_names
# save(models_crps, file = "output/models_crps.RData")
load("output/models_crps.RData")


## -----
# compute Brier score on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----

# compute empirical probabilities
models_predictive_probabilities <- lapply(predictive_samples, function(pred_sample) {
  pred_sample %>%
    mutate("predicted_conflict" = outcome > 0) %>%
    group_by(country_id, month_id) %>%
    summarise(predictive_probability = mean(predicted_conflict))
})

# merge with observations and compute summands of brier score for each model, month and country
models_brier <- lapply(models_predictive_probabilities, function(pred_probs) {
  observations %>% inner_join(pred_probs,
                              by=c("country_id"="country_id", "month_id"="month_id")) %>%
    mutate("actual_conflict" = outcome > 0) %>%
    mutate("brier" = (actual_conflict - predictive_probability)^2) %>%
    select("country_id", "month_id", "outcome", "actual_conflict", "predictive_probability", "brier")
})

# save(models_brier, file = "output/models_brier.RData")
# load("output/models_brier.RData")


## -----
# label observations as either peace ("peace"), conflict onset ("onset"), ongoing conflict ("conflict") or end of conflict ("deescalation")
# VARIANT 1: reference period is previous month
# VARIANT 2: reference period is previous year (previous 12 months)
## -----
# additionally read in all (previous) observations from 01-2017 to 12-2017 for defining the conflict situation of respective countries
# since month 457 is 01-2018, we keep months 445:456 (01-2017 to 12-2017)
observations_previous <- arrow::read_parquet("../Data/observations/cm_features.parquet") %>%
  select(month_id, country_id, ged_sb) %>%
  filter(month_id %in% 445:456) %>%
  rename(outcome = ged_sb)

observations_all <- rbind(observations_previous[,colnames(observations)], observations) %>%
  arrange(country_id, month_id)

actual_conflict <- observations_all %>%
  filter(country_id %in% country_ids) %>%
  filter(month_id %in% actuals_ids) %>%
  mutate(conflict = outcome>0)

prev_conflict <- list_cbind(lapply(1:12, function(prev_month) {
  observations_all %>%
    filter(country_id %in% country_ids) %>%
    filter(month_id %in% (actuals_ids - prev_month)) %>%
    transmute(conflict = outcome>0)
}))

names(prev_conflict) <- paste0("lag_", 1:12)

conflict_prev_month <- prev_conflict$lag_1
conflict_prev_year <- rowSums(prev_conflict)>0 # label previous period as conflict period, if at least one month had >0 fatalities

situation_month <- ifelse(!actual_conflict$conflict & !conflict_prev_month, "peace", # no conflict, no conflict
                          ifelse(actual_conflict$conflict & !conflict_prev_month, "onset", # no conflict, conflict
                                 ifelse(actual_conflict$conflict & conflict_prev_month, "conflict", # conflict, conflict
                                        "deescalation"))) # conflict, no conflict

situation_year <- ifelse(!actual_conflict$conflict & !conflict_prev_year, "peace",
                         ifelse(actual_conflict$conflict & !conflict_prev_year, "onset",
                                ifelse(actual_conflict$conflict & conflict_prev_year, "conflict",
                                       "deescalation")))

sum(situation_month != situation_year)

conflict_situations <- data.frame(actual_conflict[,c(2,3)], "situation_month" = as.vector(situation_month), "situation_year" = as.vector(situation_year))


## -----
## Plot CRPS values by conflict situation
## -----
models_crps_conflict <- models_crps %>% reduce(left_join, c("country_id", "month_id"))
names(models_crps_conflict) <- c("country_id", "month_id", model_names)
models_crps_conflict <- list(models_crps_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

models_crps_conflict_month <- models_crps_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum)
models_crps_conflict_year <- models_crps_conflict %>% select(!c("country_id", "month_id", "situation_month")) %>% group_by(situation_year) %>% summarise_all(sum)

models_crps_conflict_month[,2:ncol(models_crps_conflict_month)] <- models_crps_conflict_month[,2:ncol(models_crps_conflict_month)] / nrow(models_crps_conflict) # compute contributions to average CRPS
models_crps_conflict_year[,2:ncol(models_crps_conflict_year)] <- models_crps_conflict_year[,2:ncol(models_crps_conflict_year)] / nrow(models_crps_conflict) # compute contributions to average CRPS

colSums(models_crps_conflict_month[,2:ncol(models_crps_conflict_month)])
# create ggplot data frames
crps_month <- data.frame("CRPS" = unlist(c(models_crps_conflict_month[,2:ncol(models_crps_conflict_month)])),
                         "Situation" = rep(models_crps_conflict_month$situation_month, ncol(models_crps_conflict_month)-1),
                         "Model" = rep(names(models_crps_conflict_month)[2:ncol(models_crps_conflict_month)], each = 4))
crps_year <- data.frame("CRPS" = unlist(c(models_crps_conflict_year[,2:ncol(models_crps_conflict_year)])),
                        "Situation" = rep(models_crps_conflict_year$situation_year, ncol(models_crps_conflict_year)-1),
                        "Model" = rep(names(models_crps_conflict_year)[2:ncol(models_crps_conflict_year)], each = 4))

# exclude 3 largest models with insane CRPS values "Neg_Bin_GAM", "P_GAM", "TW_GAM"
ggplot(crps_month[-which(crps_month$Model %in% c("Neg_Bin_GAM", "P_GAM", "TW_GAM")), ], aes(fill = Situation, y = Model, x = CRPS)) + 
  geom_bar(position = "stack", stat = "identity") +
  # xlim(0.0, 1) + 
  scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
  ggtitle("Contribution to average CRPS (01-2018 to 12-2023, all countries) per conflict situation, reference: previous month") +
  xlab("Contribution to average CRPS per conflict situation")

ggplot(crps_year[-which(crps_year$Model %in% c("Neg_Bin_GAM", "P_GAM", "TW_GAM")), ], aes(fill = Situation, y = Model, x = CRPS)) + 
  geom_bar(position = "stack", stat = "identity") +
  # xlim(0.0, 1) + 
  scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
  ggtitle("Contribution to average CRPS (01-2018 to 12-2023, all countries) per conflict situation, reference: previous year") +
  xlab("Contribution to average CRPS per conflict situation")


## -----
## Plot Brier values by conflict situation
## -----
models_brier_conflict <- lapply(models_brier, function(m) m %>% select("country_id", "month_id", "brier")) %>%
  reduce(left_join, c("country_id", "month_id"))
names(models_brier_conflict) <- c("country_id", "month_id", model_names)
models_brier_conflict <- list(models_brier_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

## -----
# a) create ggplots of contributions to average Brier scores for all conflict situations
## -----
models_brier_conflict_month <- models_brier_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum)
models_brier_conflict_year <- models_brier_conflict %>% select(!c("country_id", "month_id", "situation_month")) %>% group_by(situation_year) %>% summarise_all(sum)

models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] <- models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] / nrow(models_brier_conflict) # compute contributions to average brier
models_brier_conflict_year[,2:ncol(models_brier_conflict_year)] <- models_brier_conflict_year[,2:ncol(models_brier_conflict_year)] / nrow(models_brier_conflict) # compute contributions to average brier

brier_month <- data.frame("Brier" = unlist(c(models_brier_conflict_month[,2:ncol(models_brier_conflict_month)])),
                          "Situation" = rep(models_brier_conflict_month$situation_month, ncol(models_brier_conflict_month)-1),
                          "Model" = rep(names(models_brier_conflict_month)[2:ncol(models_brier_conflict_month)], each = 4))
brier_year <- data.frame("Brier" = unlist(c(models_brier_conflict_year[,2:ncol(models_brier_conflict_year)])),
                         "Situation" = rep(models_brier_conflict_year$situation_year, ncol(models_brier_conflict_year)-1),
                         "Model" = rep(names(models_brier_conflict_year)[2:ncol(models_brier_conflict_year)], each = 4))

# exclude 3 largest models with insane Brier values "Neg_Bin_GAM", "P_GAM", "TW_GAM"
ggplot(brier_month[-which(brier_month$Model %in% c("Neg_Bin_GAM", "P_GAM", "TW_GAM")), ], aes(fill = Situation, y = Model, x = Brier)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlim(0.0, 1) + 
  scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
  ggtitle("Contribution to average Brier score (01-2018 to 12-2023, all countries) per conflict situation, reference: previous month") +
  xlab("Contribution to average Brier score per conflict situation")

ggplot(brier_year[-which(brier_year$Model %in% c("Neg_Bin_GAM", "P_GAM", "TW_GAM")), ], aes(fill = Situation, y = Model, x = Brier)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlim(0.0, 1) + 
  ggtitle("Contribution to average Brier score (01-2018 to 12-2023, all countries) per conflict situation, reference: previous year") +
  scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
  xlab("Contribution to average Brier score per conflict situation")


## -----
# b) create ggplots of contributions to average Brier scores for previously no conflict, i.e. "peace" and "onset"
## -----
models_brier_prev_peace_month <- models_brier_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum) %>% filter(situation_month %in% c("peace", "onset"))
models_brier_prev_peace_year <- models_brier_conflict %>% select(!c("country_id", "month_id", "situation_month")) %>% group_by(situation_year) %>% summarise_all(sum) %>% filter(situation_year %in% c("peace", "onset"))

models_brier_prev_peace_month[,2:ncol(models_brier_prev_peace_month)] <- models_brier_prev_peace_month[,2:ncol(models_brier_prev_peace_month)] / nrow(models_brier_conflict %>% filter(situation_month %in% c("peace", "onset"))) # compute contributions to average brier
models_brier_prev_peace_year[,2:ncol(models_brier_prev_peace_year)] <- models_brier_prev_peace_year[,2:ncol(models_brier_conflict_year)] / nrow(models_brier_conflict %>% filter(situation_year %in% c("peace", "onset"))) # compute contributions to average brier

brier_prev_peace_month <- data.frame("Brier" = unlist(c(models_brier_prev_peace_month[,2:ncol(models_brier_prev_peace_month)])),
                                     "Situation" = rep(models_brier_prev_peace_month$situation_month, ncol(models_brier_prev_peace_month)-1),
                                     "Model" = rep(names(models_brier_prev_peace_month)[2:ncol(models_brier_prev_peace_month)], each = 2))
brier_prev_peace_year <- data.frame("Brier" = unlist(c(models_brier_prev_peace_year[,2:ncol(models_brier_prev_peace_year)])),
                                    "Situation" = rep(models_brier_prev_peace_year$situation_year, ncol(models_brier_prev_peace_year)-1),
                                    "Model" = rep(names(models_brier_prev_peace_year)[2:ncol(models_brier_prev_peace_year)], each = 2))

# exclude 3 largest models with insane Brier values "Neg_Bin_GAM", "P_GAM", "TW_GAM"
ggplot(brier_prev_peace_month[-which(brier_prev_peace_month$Model %in% c("Neg_Bin_GAM", "P_GAM", "TW_GAM")), ], aes(fill = Situation, y = Model, x = Brier)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlim(0.0, 1) + 
  scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
  ggtitle("Contribution to average Brier score (01-2018 to 12-2023, all countries) in case of previous peace, reference: previous month") +
  xlab("Contribution to average Brier score per conflict situation")

ggplot(brier_prev_peace_year[-which(brier_prev_peace_year$Model %in% c("Neg_Bin_GAM", "P_GAM", "TW_GAM")), ], aes(fill = Situation, y = Model, x = Brier)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlim(0.0, 1) + 
  ggtitle("Contribution to average Brier score (01-2018 to 12-2023, all countries) in case of previous peace, reference: previous year") +
  scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
  xlab("Contribution to average Brier score per conflict situation")


## -----
# MANUAL CHECKS
## -----
# Brier score values for previous peace
sort(colSums(models_brier_prev_peace_month[,2:ncol(models_brier_prev_peace_month)]))
sort(colSums(models_brier_prev_peace_year[,2:ncol(models_brier_prev_peace_year)]))

# Poisson baseline: Where do differing Brier scores compared to zero benchmark come from?
last_comb <- observations %>% mutate("pois_param" = outcome) %>% mutate("pois_month" = month_id+1) %>%
  select(c("country_id", "pois_month", "pois_param")) %>%
  inner_join(predictive_samples$last, by=c("country_id"="country_id", "pois_month"="month_id")) %>%
  filter(pois_param == 0)

sum(last_comb$outcome>0, decreasing = TRUE) # some values are much greater than 0 although 0 should be the parameter...
sum(last_comb$outcome==0, decreasing = TRUE) # some values are much greater than 0 although 0 should be the parameter...