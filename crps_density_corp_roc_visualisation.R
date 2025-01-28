## clear environment
rm(list = ls())

## set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))

#### -------------------------------------------------------------------------------
# visualisation of the predictions from the VIEWS competition of the years 2018-2023 
#### -------------------------------------------------------------------------------


## -------------
# load packages
## -------------
library(arrow)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(pbapply)
library(scoringRules)


## ---------
# load observational data (actuals)
## ---------

# path to directory on "share-alle"
# data_path <- "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/"
data_path <- "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/" # MacOS version

# actual data from 2018-2023 in the directory
files_actuals_from18 <- list.files(data_path, pattern = "cm_actuals_\\d{4}\\.parquet", full.names = TRUE)
# read all files and bind them into a single data frame
observations_18_24 <- do.call(rbind, lapply(files_actuals_from18, arrow::read_parquet)) # observations from 2018 - 2024

## -----
# load predictive samples from benchmark models and submitted forecasts
## -----

benchmark_names <- c("boot_240", "conflictology", "last", "zero")
submissions_names <- c("bodentien_rueter_negbin", "bodentien_rueter_neuralnet", "conflictforecast_v2", "Neg_Bin_GAM", "Neg_Bin_GLMM",
                       "P_GAM", "P_GLMM", "quantile_forecast", "ShapeFinder", "submission_final_gpcmm", "submission_final_hpmm",
                       "submission_final_omm", "submission_muchlinski_thornhill", "tft", "TW_GAM", "TW_GLMM", "unito_transformer")

model_names <- c(benchmark_names, submissions_names)
n_models <- length(model_names)

appendix_names <- paste0("_cm_Y20", c(18, 19, 20, 21, 22, 23, 24))
subfolder_names <- c(18, 19, 20, 21, 22, 23, 24)
predictive_samples <- list()

for (m in 1:length(model_names)) {
  # model_files <- paste0("../Data/predictions/", model_names[m], "/", model_names[m], appendix_names, ".parquet")
  model_files <- paste0(data_path, "all_available_cm_predictions/", model_names[m], "/cm/window=Y20", subfolder_names, "/", model_names[m], appendix_names, ".parquet")
  predictive_samples[[m]] <- lapply(model_files, arrow::read_parquet)
  predictive_samples[[m]] <- bind_rows(predictive_samples[[m]])
}
names(predictive_samples) <- model_names


### CRPS/Brier-score: visualisation (Lotta)-------------------------------------------------------------------------

## -----
# compute crps on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----

country_ids <- unique(observations_18_24$country_id)
actuals_ids <- 457:528 # month_ids for 01-2018, 02-2018, ..., 12-2023
cm_pairs <- cbind(rep(country_ids, each = length(actuals_ids)),
                  rep(actuals_ids, length(country_ids)))
 
models_crps <- list()

# for (m in 1:n_models) {
#   crps_m <- apply(cm_pairs, 1, function(cm_pair) {
#     print(paste0("Benchmark/model (", m, "/", n_models, "): ", model_names[m], ", country: ", cm_pair[1], ", month: ", cm_pair[2]))
#     true_observation <- observations_18_24 %>%
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
  observations_18_24 %>% inner_join(pred_probs,
                              by=c("country_id"="country_id", "month_id"="month_id")) %>%
    mutate("actual_conflict" = outcome > 0) %>%
    mutate("brier" = (actual_conflict - predictive_probability)^2) %>%
    select("country_id", "month_id", "outcome", "actual_conflict", "predictive_probability", "brier")
})

# save(models_brier, file = "output/models_brier.RData")
load("output/models_brier.RData")


## -----
# label observations as either peace ("peace"), conflict onset ("onset"), ongoing conflict ("conflict") or end of conflict ("deescalation")
# VARIANT 1: reference period is previous month
# VARIANT 2: reference period is previous year (previous 12 months)
## -----
# additionally read in all (previous) observations from 01-2017 to 12-2017 for defining the conflict situation of respective countries
# since month 457 is 01-2018, we keep months 445:456 (01-2017 to 12-2017)
observations_17 <- arrow::read_parquet(paste0(data_path,"cm_features.parquet")) %>%
  select(month_id, country_id, ged_sb) %>%
  filter(month_id %in% 445:456) %>%
  rename(outcome = ged_sb)

observations_17_24 <- rbind(observations_17[,colnames(observations_18_24)], observations_18_24) %>%
  arrange(country_id, month_id)

actual_conflict <- observations_17_24 %>%
  filter(country_id %in% country_ids) %>%
  filter(month_id %in% actuals_ids) %>%
  mutate(conflict = outcome>0)

prev_conflict <- list_cbind(lapply(1:12, function(prev_month) {
  observations_17_24 %>%
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
## Generic plot function
## -----
create_score_decomposition_plot <- function(data, excluded_models, scoringRule, reference, conflictORpeace, xlimits = NULL) {
  scoringRuleString <- scoringRule
  if(scoringRule == "Brier"){
    scoringRuleString <- "Brier score"
  }
  p <- ggplot(data[-which(data$Model %in% excluded_models), ], aes(fill = Situation, y = Model, x = .data[[scoringRule]])) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual("legend", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa")) +
    ggtitle(paste0("Contribution to the average ", scoringRuleString ," ", conflictORpeace, " (01-2018 to 12-2023, all countries)"), 
            subtitle = paste0( "Reference: ", reference)) +
    xlab(paste0("Contribution to average ", scoringRuleString, " ", conflictORpeace))
  
  if (!is.null(xlimits)) {
    p <- p + xlim(xlimits)
  }
  return(p)
}


## -----
## Excluded models
## -----

# exclude 3 largest models with insane CRPS values "Neg_Bin_GAM", "P_GAM", "TW_GAM" for all plots
excluded_models <- c("Neg_Bin_GAM", "P_GAM", "TW_GAM")

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


# create ggplot data frames
crps_month <- data.frame("CRPS" = unlist(c(models_crps_conflict_month[,2:ncol(models_crps_conflict_month)])),
                         "Situation" = rep(models_crps_conflict_month$situation_month, ncol(models_crps_conflict_month)-1),
                         "Model" = rep(names(models_crps_conflict_month)[2:ncol(models_crps_conflict_month)], each = 4))
crps_year <- data.frame("CRPS" = unlist(c(models_crps_conflict_year[,2:ncol(models_crps_conflict_year)])),
                        "Situation" = rep(models_crps_conflict_year$situation_year, ncol(models_crps_conflict_year)-1),
                        "Model" = rep(names(models_crps_conflict_year)[2:ncol(models_crps_conflict_year)], each = 4))

print(colSums(models_crps_conflict_month[4,2:ncol(models_crps_conflict_month)] ))


# CRPS decomposition monthly
create_score_decomposition_plot(crps_month, excluded_models, "CRPS", "previous month", "per conflict situation")

# CRPS decomposition yearly
create_score_decomposition_plot(crps_year, excluded_models, "CRPS", "previous year", "per conflict situation")


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


# Brier score decomposition monthly
create_score_decomposition_plot(brier_month, excluded_models, "Brier", "previous month", "per conflict situation", c(0.0, 1))

# Brier score decomposition yearly
create_score_decomposition_plot(brier_year, excluded_models, "Brier", "previous year", "per conflict situation", c(0.0, 1))


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

# Brier score decomposition monthly
create_score_decomposition_plot(brier_prev_peace_month, excluded_models, "Brier", "previous month", "in case of previous peace", c(0.0, 1))

# Brier score decomposition yearly
create_score_decomposition_plot(brier_prev_peace_year, excluded_models, "Brier", "previous year", "in case of previous peace", c(0.0, 1))

#colSums(models_brier_prev_peace_year[,2:ncol(models_crps_conflict_month)] )



### distribution of onset probabilities : visualisation (Tobi)-------------------------------------------------------------------------

# dataframe to store onset probabilities in long format
prob_models_onset_long <- data.frame(
  model = character(),
  month_id = integer(),
  country_id = integer(),
  predictive_probability = numeric(),
  situation_month = character(),
  situation_year = character()
)

# iterate over all models except the excluded_model in models_predictive_probabilities
for (model_name in setdiff(names(models_predictive_probabilities), excluded_models)) {
  # only keep peace situations
  joined_data <- models_predictive_probabilities[[model_name]] %>%
    inner_join(conflict_situations, by = c("month_id", "country_id")) %>%
    mutate(model = model_name) # add model name
  
  # add to dataframe
  prob_models_onset_long <- bind_rows(prob_models_onset_long, joined_data)
}
tail(prob_models_onset_long)  # ÜBERPRÜFEN WAS DER LETZTE MONAT IST (NICHT DASS 2024 MIT REIN ZÄHLT)!!!!

## -----
## Generic plot functions
## -----

# density for all models -------------
create_density_all_plot <- function(data, individual, reference, onsetORpeace) {
  
  individual_string <- "Individual Models"
  
  if(individual == FALSE){
    individual_string <- "All Models"
  }
  
  p <- ggplot(data, aes(x = predictive_probability)) +
    geom_density(fill = "#90EE90", alpha = 0.6, size = 0.8) + 
    labs(
      title = paste0("Distribution of Fatality Probabilities in case of ", onsetORpeace, " (01-2018 to 12-2023, all countries)"),
      subtitle = paste0( "Reference: ", reference, ", ", individual_string),
      x = "onset probability",
      y = "density"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  if(individual == TRUE){
    p <- p + facet_wrap(~ model, scales = "free_y")  # facet for each model
  }
  
  return(p)
}


# density for each model ----------------




# boxplot for each model -----------------


## ongoing peace --------------------------------------------------------

# filter for "peace" month
ongoing_peace_prob_month_long <- prob_models_onset_long %>%
  filter(situation_month == "peace")


ongoing_peace_prob_year_long <- prob_models_onset_long %>%
  filter(situation_year == "peace")

onset_prob_year_long <- prob_models_onset_long %>%
  filter(situation_year == "onset")


create_density_all_plot(ongoing_peace_prob_month_long,individual = TRUE, "Previous month", "peace")
create_density_all_plot(ongoing_peace_prob_year_long,individual = TRUE, "Previous year", "peace")


ggplot(ongoing_peace_prob_month_long, aes(x = predictive_probability)) +
  geom_density(fill = "#90EE90", alpha = 0.6, size = 0.8) + 
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: Individual Models",
    subtitle = "Ongoing Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# density plot for every observation over all models
ggplot(ongoing_peace_prob_month_long, aes(x = predictive_probability)) +
  geom_density(color="black",fill="#90EE90", size = 1, alpha = 0.6) +
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: All Models",
    subtitle = "Ongoing Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_bw()


# density for each model
ggplot(prob_long_ongoing_peace, aes(x = probability_gr_0)) +
  geom_density(fill = "#90EE90", alpha = 0.6, size = 0.8) + 
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: Individual Models",
    subtitle = "Ongoing Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )
