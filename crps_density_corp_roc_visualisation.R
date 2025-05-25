## clear environment
rm(list = ls())

## set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))

# renv::snapshot()
# renv::update() 
# renv::restore()

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
library(reliabilitydiag)
library(gridExtra)
library(grid)
library(precrec)
library(triptych)
library(ggpattern)
library(cowplot)

#devtools::install_github("aijordan/reliabilitydiag")

## ---------
# load observational data (actuals)
## ---------

# path to directory on "share-alle"
# data_path <- "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/"
data_path <- "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/" # MacOS version


# Erkenne das Betriebssystem
os <- Sys.info()["sysname"]

# Setze den Pfad abhängig vom Betriebssystem
data_path <- ifelse(os == "Windows",
                    "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/",  # Windows
                    "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/")  # macOS/Linux

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


## -----
## Models
## -----
model_names <- c(benchmark_names, submissions_names)
n_models <- length(model_names)


## -----
## load data
## -----

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


## -----
## excluded models
## -----

# exclude some models for all plots and analysis
excluded_models <- c("P_GLMM","TW_GLMM","submission_final_gpcmm","submission_final_hpmm","Neg_Bin_GAM", "P_GAM", "TW_GAM")

# later on in the CORP, Murphy diagram section there will be even more models removed.
# but this is only for the detailed analyis. 
# For the general analysis only "excluded_models" are removed




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

# compute empirical probabilities for the binary onset event
models_predictive_probabilities <- lapply(predictive_samples, function(pred_sample) {
  pred_sample %>%
    mutate("predicted_conflict" = outcome > 0) %>%
    group_by(country_id, month_id) %>%
    summarise(predictive_probability = mean(predicted_conflict))
})

# merge models_predictive_probabilities and models_crps into new list "models_scoring_rules"
models_scoring_rules <- list()
for (m in 1:n_models) {
  models_scoring_rules[[m]] <- models_crps[[m]] %>%
    left_join(models_predictive_probabilities[[m]], by = c("country_id", "month_id")) %>%
    rename(onset_prob_pred = predictive_probability)
}
names(models_scoring_rules) <- model_names

# add the acutal observations to list
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    left_join(observations_18_24 %>%
                select(country_id, month_id, actual = outcome),
              by = c("country_id", "month_id"))
})

# compute summands of brier score for each model, month and country
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    mutate(
      actual_conflict = actual > 0,
      brier_onset = (actual_conflict - onset_prob_pred)^2
    )
})

# remove lists that are not longer needed
rm(models_predictive_probabilities, models_crps)




## -----
# compute log-score on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----
# modified log-score to cure log(0)
safe_log_score <- function(actual, p, eps = 1e-12) {
  p_clipped <- pmin(pmax(p, eps), 1 - eps)
  return(- (actual * log(p_clipped) +
              (1 - actual) * log(1 - p_clipped)))
}

# compute log-scores
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    mutate(
      log_score_onset = - (actual_conflict * log(onset_prob_pred) +
                             (1 - actual_conflict) * log(1 - onset_prob_pred)),
      log_score_eps_onset = safe_log_score(actual_conflict, onset_prob_pred)
    )
})






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
## Plot CRPS values by conflict situation
## -----
models_crps_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "crps")) %>%
  reduce(left_join, c("country_id", "month_id"))


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
models_brier_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "brier_onset")) %>%
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







## -----
## Plot log-score values by conflict situation for the onset problem (y \in {0,1})
## -----
models_logscore_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "log_score_eps_onset")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_logscore_conflict) <- c("country_id", "month_id", model_names)
models_logscore_conflict <- list(models_logscore_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

models_logscore_conflict_month <- models_logscore_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum)
models_logscore_conflict_year <- models_logscore_conflict %>% select(!c("country_id", "month_id", "situation_month")) %>% group_by(situation_year) %>% summarise_all(sum)

models_logscore_conflict_month[,2:ncol(models_logscore_conflict_month)] <- models_logscore_conflict_month[,2:ncol(models_logscore_conflict_month)] / nrow(models_logscore_conflict) # compute contributions to average CRPS
models_logscore_conflict_year[,2:ncol(models_logscore_conflict_year)] <- models_logscore_conflict_year[,2:ncol(models_logscore_conflict_year)] / nrow(models_logscore_conflict) # compute contributions to average CRPS


# create ggplot data frames
logscore_month <- data.frame("logscore" = unlist(c(models_logscore_conflict_month[,2:ncol(models_logscore_conflict_month)])),
                         "Situation" = rep(models_logscore_conflict_month$situation_month, ncol(models_logscore_conflict_month)-1),
                         "Model" = rep(names(models_logscore_conflict_month)[2:ncol(models_logscore_conflict_month)], each = 4))
logscore_year <- data.frame("logscore" = unlist(c(models_logscore_conflict_year[,2:ncol(models_logscore_conflict_year)])),
                        "Situation" = rep(models_logscore_conflict_year$situation_year, ncol(models_logscore_conflict_year)-1),
                        "Model" = rep(names(models_logscore_conflict_year)[2:ncol(models_logscore_conflict_year)], each = 4))

print(colSums(models_logscore_conflict_month[4,2:ncol(models_logscore_conflict_month)] ))


# CRPS decomposition monthly
create_score_decomposition_plot(logscore_month, excluded_models, "logscore", "previous month", "per conflict situation")

# CRPS decomposition yearly
create_score_decomposition_plot(logscore_year, excluded_models, "logscore", "previous year", "per conflict situation")









### distribution of onset probabilities : visualisation (Tobi)-------------------------------------------------------------------------

# dataframe to store onset probabilities in long format
prob_models_all_situation_long <- data.frame(
  model = character(),
  month_id = integer(),
  country_id = integer(),
  onset_prob_pred = numeric(),
  situation_month = character(),
  situation_year = character()
)

# iterate over all models except the excluded_model in models_predictive_probabilities
for (model_name in setdiff(names(models_scoring_rules), excluded_models)) {
  
  joined_data <- models_scoring_rules[[model_name]] %>%
    inner_join(conflict_situations, by = c("month_id", "country_id")) %>%
    mutate(model = model_name) # add model name
  
  # add to dataframe
  prob_models_all_situation_long <- bind_rows(prob_models_all_situation_long, joined_data)
}
# tail(prob_models_all_situation_long)


#length(which(prob_models_all_situation_long$situation_year == "onset"))/(length(model_names)-3)


## -----
## Generic plot functions
## -----

## density -----
##
## Achtung: dichte wird geglättet (siehe zero modell) -> adjust einstellen
create_density_plot <- function(data, individual, reference, onsetORpeace) {
  
  individual_string <- "individual models"
  fill_color <- ""
  
  if(individual == FALSE){
    individual_string <- "all models"
  }
  
  if(onsetORpeace == "peace"){
    fill_color <- "#4664aa"
  } else if(onsetORpeace == "onset"){
    fill_color <- "#df9b1b"
  } else if(onsetORpeace == "previous peace") {
    fill_color <- "#556B2F" #  #5F2F4F
  }
  
  p <- ggplot(data, aes(x = onset_prob_pred)) +
    geom_density(fill = fill_color, alpha = 0.6, linewidth = 0.8, adjust = 0.2) + 
    labs(
      title = paste0("Distribution of fatality probabilities in case of ", onsetORpeace, " (01-2018 to 12-2023, all countries)"),
      subtitle = paste0( "Reference: ", reference, ", ", individual_string),
      x = "predicted probability",
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

## boxpolot -----
create_box_plot <- function(data, individual, reference, onsetORpeace) {
  
  individual_string <- "individual models"
  fill_color <- ""
  
  if(individual == FALSE){
    individual_string <- "all models"
  }
  
  if(onsetORpeace == "peace"){
    fill_color <- "#4664aa"
  } else if(onsetORpeace == "onset"){
    fill_color <- "#df9b1b"
  } else if(onsetORpeace == "previous peace") {
    fill_color <- "#556B2F" #  #5F2F4F
  }
  
  p <- ggplot(data, aes(x = onset_prob_pred)) +
    geom_boxplot(fill = fill_color, alpha = 0.6, color = "black", linewidth = 0.6, outlier.color = "red", outlier.size = 1, width = 0.4) +
    labs(
      title = paste0("Boxplot of fatality probabilities in case of ", onsetORpeace, " (01-2018 to 12-2023, all countries)"),
      subtitle = paste0( "Reference: ", reference, ", ", individual_string),
      x = "predicted probability",
      y = NULL
    ) +
    scale_y_continuous(
      limits = c(-0.4, 0.4)
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



## -----
## ongoing peace
## -----

# filter for "peace" month
ongoing_peace_prob_month_long <- prob_models_all_situation_long %>%
  filter(situation_month == "peace")

ongoing_peace_prob_year_long <- prob_models_all_situation_long %>%
  filter(situation_year == "peace")

## density-plots -----
# previous month
#create_density_plot(ongoing_peace_prob_month_long,individual = TRUE, "previous month", "peace")
# previous year
#create_density_plot(ongoing_peace_prob_year_long,individual = TRUE, "previous year", "peace")

# ## boxplots -----
# # previous month
# create_box_plot(ongoing_peace_prob_month_long,individual = TRUE, "previous month", "peace")
# # previous year
# create_box_plot(ongoing_peace_prob_year_long,individual = TRUE, "previous year", "peace")



## -----
## onset
## -----

# filter for "onset" month
onset_prob_month_long <- prob_models_all_situation_long %>%
  filter(situation_month == "onset")

onset_prob_year_long <- prob_models_all_situation_long %>%
  filter(situation_year == "onset")

## density-plots -----
# previous month
#create_density_plot(onset_prob_month_long,individual = TRUE, "previous month", "onset")
# previous year
#create_density_plot(onset_prob_year_long,individual = TRUE, "previous year", "onset")

# ## boxplots -----
# # previous month
# create_box_plot(onset_prob_month_long,individual = TRUE, "previous month", "onset")
# # previous year
# create_box_plot(onset_prob_year_long,individual = TRUE, "previous year", "onset")



## -----
## previous peace 
## -----

# filter for "onset" or "peace" month
prev_peace_prob_month_long <- prob_models_all_situation_long %>%
  filter(situation_month == "peace" | situation_month == "onset")

prev_peace_prob_year_long <- prob_models_all_situation_long %>%
  filter(situation_year == "peace" | situation_year == "onset")

## density-plots -----
# previous month
density_previous_peace_month_curves <- create_density_plot(prev_peace_prob_month_long,individual = TRUE, "previous month", "previous peace")
density_previous_peace_month_curves
# previous year
#create_density_plot(prev_peace_prob_year_long,individual = TRUE, "previous year", "previous peace")

# ggsave("plots_tobi/density_previous_peace_month.png", 
#        plot = density_previous_peace_month_curves, width = 20, height = 12, dpi = 300, bg = "white")

# ## boxplots -----
# # previous month
# create_box_plot(prev_peace_prob_month_long,individual = TRUE, "previous month", "previous peace")
# # previous year
# create_box_plot(prev_peace_prob_year_long,individual = TRUE, "previous year", "previous peace")













## -----
## previous peace density-decomposition-month for predicted prob. > 0!
## -----
# filter dataset
decomp_prev_peace_prob_month_long <- prev_peace_prob_month_long %>%
  filter(onset_prob_pred > 0)

# ideal forecast would have everything green for p close to 0 and 
# everything yellow for high probabilities (p close to 1).
density_decomposition_previous_peace_month_curves <- ggplot(data = decomp_prev_peace_prob_month_long, aes(x = onset_prob_pred)) +
  # "peace" and "onset" density (i.e. "previous peace")
  geom_density(aes(y = after_stat(density), fill = "previous peace"), 
               alpha = 0.9, 
               size = 0.8,
               adjust = 0.2) +
  
  # "onset" density weighted relative to overall dist of prev_peace_prob_month_long
  geom_density(data = filter(decomp_prev_peace_prob_month_long, situation_month == "onset"),
               aes(y = after_stat(density) * nrow(filter(decomp_prev_peace_prob_month_long, situation_month == "onset")) / nrow(decomp_prev_peace_prob_month_long),
                   fill = "onset"),
               alpha = 0.8,
               size = 0.8,
               adjust = 0.2) +
  
  facet_wrap(~ model, scales = "free_y") +
  
  scale_y_sqrt() +
  
  scale_fill_manual(name="conflict situation",
                    values=c("previous peace"="#556B2F","onset"="#df9b1b")) + 
  
  labs(
    title = "Distribution of fatality probabilities greater zero in case of preavious peace (01-2018 to 12-2023, all countries)",
    subtitle = "Reference: previous month, individual models",
    x = "predicted probability",
    y = "density"
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

density_decomposition_previous_peace_month_curves

# ggsave("plots_tobi/density_decomposition_previous_peace_month.png", 
#        plot = density_decomposition_previous_peace_month_curves, width = 20, height = 12, dpi = 300, bg = "white")





## -----
## previous peace density-decomposition-month for predicted prob. >= 0 
## -----

# ideal forecast would have everything green for p close to 0 and 
# everything yellow for high probabilities (p close to 1).
density_decomposition_previous_peace_month_curves_with_zero <- ggplot(data = prev_peace_prob_month_long, aes(x = onset_prob_pred)) +
  # "peace" and "onset" density (i.e. "previous peace")
  geom_density(aes(y = after_stat(density), fill = "previous peace"), 
               alpha = 0.9, 
               size = 0.8,
               adjust = 0.2) +
  
  # "onset" density weighted relative to overall dist of prev_peace_prob_month_long
  geom_density(data = filter(decomp_prev_peace_prob_month_long, situation_month == "onset"),
               aes(y = after_stat(density) * nrow(filter(decomp_prev_peace_prob_month_long, situation_month == "onset")) / nrow(decomp_prev_peace_prob_month_long),
                   fill = "onset"),
               alpha = 0.8,
               size = 0.8,
               adjust = 0.2) +
  
  facet_wrap(~ model, scales = "free_y") +
  
  scale_y_sqrt() +
  
  scale_fill_manual(name="conflict situation",
                    values=c("previous peace"="#556B2F","onset"="#df9b1b")) + 
  
  labs(
    title = "Distribution of fatality probabilities in case of preavious peace (01-2018 to 12-2023, all countries)",
    subtitle = "Reference: previous month, individual models",
    x = "predicted probability",
    y = "density"
  ) +
  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

density_decomposition_previous_peace_month_curves_with_zero



















## Proportion of fatality probability exceeding x% (from previous peace for >0%, >1%, >2%, >5%) -----
thresholds <- c(0, 0.01, 0.02, 0.05)

prob_proportion_exceeding_thresh_df <- prev_peace_prob_month_long %>%
  group_by(model) %>%
  summarise(across(onset_prob_pred, list(
    `t1` = ~ mean(. > thresholds[1], na.rm = TRUE), # . stands for column predictive_probability
    `t2` = ~ mean(. > thresholds[2], na.rm = TRUE), # mean is sufficient due to True = 1, False = 0
    `t3` = ~ mean(. > thresholds[3], na.rm = TRUE),
    `t4` = ~ mean(. > thresholds[4], na.rm = TRUE)
  ))) %>%
  pivot_longer(cols = -model, names_to = "threshold", values_to = "percentage") %>%
  pivot_wider(names_from = model, values_from = percentage) %>%
  mutate(threshold = thresholds)

for (i in seq_along(thresholds)) {
  thresh <- thresholds[i]
  
  cat("\nProb. exceeding", thresh, ":\n")
  print(setNames(colSums(prob_proportion_exceeding_thresh_df[i, 2:ncol(prob_proportion_exceeding_thresh_df)]),
                 colnames(prob_proportion_exceeding_thresh_df)[2:ncol(prob_proportion_exceeding_thresh_df)]))
}








## -----
## plot analysis inspired by the triptych (Dimitriades et. al.)
## -----
## at first plots for all 13 models are made, after that the final
# 8 models for the in depth analysis are chosen

## CORP (reliability diagram) -----
prev_peace_prob_month_long_binary_actual <- prev_peace_prob_month_long %>%
  mutate(actual = ifelse(actual >= 1, 1, 0))

# list to store the CORP
corp_plots_list <- list()

for (model_name in setdiff(names(models_scoring_rules), excluded_models)) {
  
  reliability_data <- prev_peace_prob_month_long_binary_actual %>%
    filter(model == model_name)
  
  r <- reliabilitydiag(x = reliability_data$onset_prob_pred, y = reliability_data$actual)
  plot <- autoplot(r) + ggplot2::labs(
    title = model_name,
    x = "predicted probability",
    y = "CEP") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13)
    )
  
  # hier potentiell noch mehr abspeichern
  corp_plots_list[[model_name]] <- list(plot = plot)
  
}

# extract all plots from list -> returns list
plots_corp <- lapply(corp_plots_list, function(x) x$plot)

# plot with grid arrange
tg_corp <- textGrob('Reliability diagram of predicted fatality prob. in case of previous peace (01-2018 to 12-2023, all countries)', gp = gpar(fontsize = 18, fontface = 'bold'))
sg_corp <- textGrob('Reference: previous month, individual models', gp = gpar(fontsize = 15))
margin <- unit(0.5, "line")
grided_corp <- gridExtra::grid.arrange(grobs = plots_corp, ncol = 5)
corp_previous_peace_month_curves <- gridExtra::grid.arrange(tg_corp, sg_corp, grided_corp,
                                                            heights = unit.c(grobHeight(tg_corp) + 1.2*margin, 
                                                                             grobHeight(sg_corp) + margin, 
                                                                             unit(1,"null")))
                                       
#ggsave("plots_tobi/corp_previous_peace_month.png", plot = corp_previous_peace_month_curves, width = 20, height = 12, dpi = 300)

## ROC -----
# list to store the roc diagrams
roc_plots_list <- list()

for (model_name in setdiff(names(models_scoring_rules), excluded_models)) {
  
  roc_data <- prev_peace_prob_month_long_binary_actual %>%
    filter(model == model_name)
  
  mm <- mmdata(roc_data$onset_prob_pred, roc_data$actual)
  
  plot <- autoplot(evalmod(mm), curvetype = "ROC") + 
    ggplot2::labs(
      title = model_name,
      x = "FPR",
      y = "TPR") +
    ggplot2::geom_line(size = 1) + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13)
    ) + 
    ggplot2::scale_color_manual(values=c("#556B2F"))
  
  roc_plots_list[[model_name]] <- list(plot = plot)
  
}

# extract all plots from list -> returns list
plots_roc <- lapply(roc_plots_list, function(x) x$plot)

# plot with grid arrange
tg_roc <- textGrob('ROC diagram of predicted fatality prob. in case of previous peace (01-2018 to 12-2023, all countries)', gp = gpar(fontsize = 18, fontface = 'bold'))
sg_roc <- textGrob('Reference: previous month, individual models', gp = gpar(fontsize = 15))
grided_roc <- gridExtra::grid.arrange(grobs = plots_roc, ncol = 5)
roc_previous_peace_month_curves <- gridExtra::grid.arrange(tg_roc, sg_roc, grided_roc,
                                                           heights = unit.c(grobHeight(tg_roc) + 1.2*margin, 
                                                                            grobHeight(sg_roc) + margin, 
                                                                            unit(1,"null")))
                                      
#ggsave("plots_tobi/roc_previous_peace_month.png", plot = roc_previous_peace_month_curves, width = 20, height = 12, dpi = 300)                       



## -----
## murphy diagramm
## -----
murphy_data <- prev_peace_prob_month_long_binary_actual %>%
  select(1:6,8) %>%
  pivot_wider(names_from = model, values_from = onset_prob_pred) %>%
  select(5:ncol(.)) 

mr_conflict <- murphy(murphy_data, 
                      y_var = "actual")


df_est_conflict <- estimates(mr_conflict)
murphy_plot <- ggplot(df_est_conflict) + 
  geom_path(aes(x = knot, y = mean_score, col = forecast), linewidth = 0.6, alpha = 1) + 
  labs(title = "Murphy diagram of predicted fatality prob. in case of previous peace (01-2018 to 12-2023, all countries)",
       subtitle = "Reference: previous month, individual models")
murphy_plot

#ggsave("plots_tobi/murphy_previous_peace_month__selected.png", plot = selected_murphy_plot, width = 10, height = 9, dpi = 300)                       





## -----------------------------------------------------
## plots for selection of 8 models
##------------------------------------------------------
selected_models <- c("last","conflictology",
                     "conflictforecast_v2", "submission_muchlinski_thornhill", "submission_final_omm", "unito_transformer", "Neg_Bin_GLMM")
#"zero",

model_colors <- c(
  "zero"                        = "#009682",  # green  
  "last"                        = "#4664aa",  # blue
  "conflictology"              = "#23a1e0",  # maygreen
  "conflictforecast_v2"        = "black", # grey
  "submission_muchlinski_thornhill" = "#df9b1b",  # orange
  "submission_final_omm"       = "#8cb63c",  # yellow
  "unito_transformer"          = "#a22223",  # red
  "Neg_Bin_GLMM"               = "#a3107c"   # purple
)

model_labels <- c(
  "zero"                          = "VIEWS Zero",
  "last"                          = "VIEWS Last Month",
  "conflictology"                = "VIEWS Conflictology",
  "conflictforecast_v2"          = "Conflict Forecast",
  "submission_muchlinski_thornhill" = "MT Hurdle GAM",
  "submission_final_omm"         = "Observed MM",
  "unito_transformer"            = "UNITO Transformer",
  "Neg_Bin_GLMM"                 = "NegBin GLMM",
  #
  "tft" = "CCEW tft",
  "ShapeFinder" = "Shape Finder",
  "quantile_forecast" = "Quantile Forecast",
  "boot_240" = "VIEWS Bootstrap",
  "bodentien_rueter_neuralnet" = "BR Neural Net",
  "bodentien_rueter_negbin" = "BR NegBin"
)

theme_setup <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
  axis.title = ggplot2::element_text(size = 13),
  legend.title = ggplot2::element_text(size = 13),
  axis.text = ggplot2::element_text(size = 12),
  legend.text = ggplot2::element_text(size = 11)
)


## -----
## CRPS decomposition
## -----

## Attention: not only for selected models but for all models without excluded ones!!

crps_month_selected_models <- crps_month

# Rename models
crps_month_selected_models$Model <- recode(
  crps_month_selected_models$Model,
  !!!model_labels
)

crps_conflict_situation_plot <- ggplot(crps_month_selected_models[-which(crps_month_selected_models$Model %in% excluded_models), ], 
                                       aes(fill = Situation, y = Model, x = .data[["CRPS"]])) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean CRPS per Conflict Situation",
       x = "Mean CRPS") +
  scale_fill_manual("Situation", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa"),
                    labels = c("conflict" = "Conflict", "deescalation" = "Deescalation", "onset" = "Onset", "peace" = "Peace")) +
  theme_setup +
  theme(
    axis.ticks.y = element_blank(),
  )

crps_conflict_situation_plot







## -----
## murphy diagram
## -----
mr_conflict_subset <- murphy(subset(murphy_data, select = c(selected_models,"actual")), 
                      y_var = "actual")

df_est_conflict_subset <- estimates(mr_conflict_subset)
selected_murphy_plot <- ggplot(df_est_conflict_subset) + 
  geom_path(aes(x = knot, y = mean_score, col = forecast), linewidth = 1, alpha = 0.8) + 
  scale_color_manual(
    values = model_colors,
    breaks = names(model_labels),
    labels = model_labels
  ) +
  theme_bw() +
  labs(
    title = "Murphy Diagram",
    x = expression(paste("Parameter ", theta)),
    y = "Mean elementary score",
    color = ""
  ) + 
  theme_setup + 
  theme(legend.position = "bottom")

selected_murphy_plot




### RECHTECKIG ABSPEICHERN





## -----
## ROC curves
## -----

## -
## https://cran.r-project.org/web/packages/precrec/vignettes/introduction.html
## -

roc_data <- prev_peace_prob_month_long_binary_actual

# filter relevant models
roc_data_filtered <- roc_data %>% 
  filter(model %in% selected_models)

# create list of scores
score_list <- lapply(selected_models, function(m) {
  roc_data_filtered %>%
    filter(model == m) %>%
    pull(onset_prob_pred)
})

# create list of labels
label_list <- lapply(selected_models, function(m) {
  roc_data_filtered %>%
    filter(model == m) %>%
    pull(actual)
})

# join lists
scores_joined <- join_scores(score_list)
labels_joined <- join_labels(label_list)

# roc-curve data
mm_selected_models <- mmdata(
  scores   = scores_joined,
  labels   = labels_joined,
  modnames = selected_models
)


roc_curve_selected_models <- autoplot(evalmod(mm_selected_models), curvetype = "ROC") + 
  ggplot2::geom_line(size = 1, alpha = 0.8) + 
  ggplot2::scale_color_manual(
    values = model_colors,
    breaks = names(model_labels),
    labels = model_labels
  ) +
  ggplot2::labs(
    title = "ROC Curve",
    x = "FAR",
    y = "HR",
    color = ""
  ) +
  theme_setup +
  ggplot2::theme(legend.position = "bottom")

roc_curve_selected_models





### RECHTECKIG ABSPEICHERN






## -----
## Reliability diagram
## -----

## function to create custom reliability diagram
create_reliability_diag <- function(data, forecast_model) {
  
  reliability_data_selected <- data %>%
    filter(model == forecast_model)
  
  r_selected <- reliabilitydiag(x = reliability_data_selected$onset_prob_pred, y = reliability_data_selected$actual)
  
  data_estim <- estimates(reliability(x = reliability_data_selected$onset_prob_pred, y = reliability_data_selected$actual))
  
  # delete duplicate rows
  df_clean <- data_estim %>% distinct()
  
  # make sure that CEP is in right order
  df_clean <- df_clean[order(df_clean$CEP), ]
  
  ## create segment data frame
  df_segments <- data.frame(
    x = numeric(),
    x_end = numeric(),
    CEP = numeric(),
    CEP_end = numeric()
  )
  
  path_color = model_colors[forecast_model]
  #path_color = "green"
  
  
  # extract all segments (slope 0)
  for (i in 1:(nrow(df_clean) - 1)) {
    if (df_clean$CEP[i] == df_clean$CEP[i + 1]) {
      x_min <- min(df_clean$x[i], df_clean$x[i + 1])
      x_max <- max(df_clean$x[i], df_clean$x[i + 1])
      cep <- df_clean$CEP[i]
      cep_end <- df_clean$CEP[i+1]
      
      df_segments <- rbind(df_segments, data.frame(
        x = x_min,
        x_end = x_max,
        CEP = cep,
        CEP_end = cep_end
      ))
    }
  }
  
  p <- autoplot(r_selected) + ggplot2::labs(
    title = model_labels[model_name],
    x = "Forecast value",
    y = "CEP") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$x, 
        y = .data$CEP, 
        xend = .data$x_end, 
        yend = .data$CEP_end),
      data = df_segments,
      linewidth = 1.8,
      colour = path_color
    ) +
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$CEP
      ),
      data = data_estim,
      linewidth = 1,
      colour = path_color
    ) 
  #+ 
   # ggplot2::scale_color_manual(values = c("value" = path_color))
  
  return(p)
  
}


# list to store the CORP
corp_plots_list <- list()

for (model_name in selected_models) {
  
  # reliability_data_selected <- prev_peace_prob_month_long_binary_actual %>%
  #   filter(model == model_name)
  # 
  # r_selected <- reliabilitydiag(x = reliability_data_selected$onset_prob_pred, y = reliability_data_selected$actual)
  # plot <- autoplot(r_selected) + ggplot2::labs(
  #   title = model_labels[model_name],
  #   x = "Forecast value",
  #   y = "CEP") +
  #   ggplot2::theme(
  #     plot.title = ggplot2::element_text(size = 14, hjust = 0.5)
  #   ) +
  #   ggplot2::geom_path(
  #     mapping = ggplot2::aes(
  #       x = .data$x,
  #       y = .data$CEP,
  #       col = .data$forecast
  #     ),
  #     data = data
  #   )
  
  
  plot <- create_reliability_diag(prev_peace_prob_month_long_binary_actual, model_name)
  
  # hier potentiell noch mehr abspeichern
  corp_plots_list[[model_name]] <- list(plot = plot)
  
}

# extract all plots from list -> returns list
plots_corp <- lapply(corp_plots_list, function(x) x$plot)

plots_corp["submission_muchlinski_thornhill"]












reliability_data_selected <- prev_peace_prob_month_long_binary_actual %>%
  filter(model == "submission_muchlinski_thornhill")

r_selected <- reliabilitydiag(x = reliability_data_selected$onset_prob_pred, y = reliability_data_selected$actual)

data_estim <- estimates(reliability(x = reliability_data_selected$onset_prob_pred, y = reliability_data_selected$actual))



df_clean <- data_estim %>% distinct()

# make sure that CEP is in right order
df_clean <- df_clean[order(df_clean$CEP), ]

df_segments <- data.frame(
  x = numeric(),
  x_end = numeric(),
  CEP = numeric(),
  CEP_end = numeric()
)

# extract all segments (slope 0)
for (i in 1:(nrow(df_clean) - 1)) {
  if (df_clean$CEP[i] == df_clean$CEP[i + 1]) {
    x_min <- min(df_clean$x[i], df_clean$x[i + 1])
    x_max <- max(df_clean$x[i], df_clean$x[i + 1])
    cep <- df_clean$CEP[i]
    cep_end <- df_clean$CEP[i+1]
    
    df_segments <- rbind(df_segments, data.frame(
      x = x_min,
      x_end = x_max,
      CEP = cep,
      CEP_end = cep_end
    ))
  }
}


autoplot(r_selected) + ggplot2::labs(
  title = model_labels[model_name],
  x = "Forecast value",
  y = "CEP") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
    legend.position = "none"
  ) +
  ggplot2::geom_segment(
    mapping = ggplot2::aes(
      x = .data$x, 
      y = .data$CEP, 
      xend = .data$x_end, 
      yend = .data$CEP_end),
    data = df_segments,
    linewidth = 1.8,
    colour = "green"
  ) +
  ggplot2::geom_path(
    mapping = ggplot2::aes(
      x = .data$x,
      y = .data$CEP,
      col = .data$forecast
    ),
    data = data_estim,
    linewidth = 1
  ) + 
  ggplot2::scale_color_manual(values = c("value" = "green"))
  
  



















# Make three plots.
# We set left and right margins to 0 to remove unnecessary spacing in the
# final plot arrangement.
p1 <- plots_corp["conflictforecast_v2"]$conflictforecast_v2 +
  theme(plot.margin = unit(c(6,0,6,0), "pt"))
p2 <- plots_corp["zero"]$zero +
  theme(plot.margin = unit(c(6,0,6,0), "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  ylab("")
p3 <- plots_corp["last"]$last +
  theme(plot.margin = unit(c(6,0,6,0), "pt")) + ylab("")


# arrange the three plots in a single row
prow <- plot_grid( p1 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   p3 + theme(legend.position="none"),
                   align = 'vh',
                   hjust = -1,
                   nrow = 1
)

# extract the legend from one of the plots
# (clearly the whole thing only makes sense if all plots
# have the same legend, so we can arbitrarily pick one.)
legend_b <- get_legend(p1 + theme(legend.position="bottom"))

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
p




















## -----
## MSC-DSC-plots
## -----

brier_decomposition_results <- map2(
  .x = models_scoring_rules,
  .y = names(models_scoring_rules),
  .f = function(df, model_name) {
    # caclulate MSC-DSC-UNC 
    #@param score 
    # A string specifying the score function.
    #'   One of: `"Brier_score"` (default), `"log_score"`, `"MR_score"`.
    result <- mcbdsc(df %>% select(onset_prob_pred), y = df$actual_conflict, score="Brier_score")
    
    # Extrahiere estimates + Modellname
    estimates(result) %>%
      mutate(model = model_name)
  }
)

# filter for selected models
brier_decomposition_results <- brier_decomposition_results[selected_models]

## important remark for BRIER score:
# UNC = variance of X~Ber(p=E(y)=mean(y)) -> p(1-p)
#p <- models_scoring_rules$boot_240$actual_conflict
#mean(p)*(1-mean(p))

# combine to one dataframe
brier_score_decomposition <- bind_rows(brier_decomposition_results) %>%
  select(model, mean_score, MCB, DSC, UNC) %>%
  rename(MeanScore = mean_score)

brier_score_decomposition_barplot <- brier_score_decomposition %>%
  mutate(score_invisible = MeanScore, gap = 0)

# data frame with format for the barchart
brier_score_decomposition_barplot <- brier_score_decomposition_barplot %>%
  arrange(MeanScore) %>%  # sort by MeanScore
  pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap),
               names_to = "component",
               values_to = "value") %>%
  mutate(class = case_when(
    component %in% c("MeanScore") ~ "SCORE",
    component %in% c("MCB", "UNC") ~ "MCB_UNC",
    component %in% c("DSC", "score_invisible") ~ "DSC_score", 
    component %in% c("gap") ~ "GAP"
  )) %>%
  select(model, class, component, value)

brier_score_decomposition_barplot <- brier_score_decomposition_barplot %>%
  mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))


# dataset for model orderbased on MeanScore
brier_levs <- brier_score_decomposition_barplot %>%
  filter(component == "MeanScore") %>%
  arrange(desc(value)) %>%
  pull(model)

# data for the plot
brier_score_decomposition_barplot <- brier_score_decomposition_barplot %>%
  mutate(
    class_group = case_when(
      class == "MCB_UNC"   ~ 1,
      class == "DSC_score"       ~ 2,
      class == "SCORE" ~ 3,
      class == "GAP" ~ 4
    ),  
    model = factor(model, levels = brier_levs, ordered = TRUE)
  )

brier_score_decomposition_barplot <- brier_score_decomposition_barplot %>%
  mutate(
    component = factor(component, levels = c("gap","MCB", "DSC","score_invisible","UNC", "MeanScore"))
  )

brier_score_decomposition_barplot <- brier_score_decomposition_barplot %>%
  mutate(
    wdth = case_when(
      class_group == 3 ~ 4,
      class_group %in% c(1, 2) ~ 1,
      class_group == 4 ~ 10
    )
  )

brier_score_decomposition_barplot$model <- recode(
  brier_score_decomposition_barplot$model,
  !!!model_labels
)



brier_score_decomposition <- ggplot(brier_score_decomposition_barplot, aes(x = class_group, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack", width = brier_score_decomposition_barplot$wdth) +
  facet_grid(model ~ ., switch = "y") +
  coord_flip() +
  scale_fill_manual(
    values = c("gap"="white","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),
    ###################################
    breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
    ##################################
    name = ""
  ) +
  labs(title = "Mean Brier Score Decomposition Plot") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 11),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "points"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
    plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
    legend.title = ggplot2::element_text(size = 13),
    legend.text = ggplot2::element_text(size = 11),
    axis.line.x = element_blank(),
    axis.line.y = element_blank()
  )

brier_score_decomposition


