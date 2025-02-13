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
prob_models_all_situation_long <- data.frame(
  model = character(),
  month_id = integer(),
  country_id = integer(),
  predictive_probability = numeric(),
  situation_month = character(),
  situation_year = character()
)

# iterate over all models except the excluded_model in models_predictive_probabilities
for (model_name in setdiff(names(models_predictive_probabilities), excluded_models)) {
  
  joined_data <- models_predictive_probabilities[[model_name]] %>%
    inner_join(conflict_situations, by = c("month_id", "country_id")) %>%
    mutate(model = model_name) # add model name
  
  # add to dataframe
  prob_models_all_situation_long <- bind_rows(prob_models_all_situation_long, joined_data)
}
# tail(prob_models_all_situation_long)

# add outcome values
prob_models_all_situation_long <- prob_models_all_situation_long %>%
  left_join(select(actual_conflict, month_id, country_id, outcome), by = c("month_id","country_id"))


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
  
  p <- ggplot(data, aes(x = predictive_probability)) +
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
  
  p <- ggplot(data, aes(x = predictive_probability)) +
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
create_density_plot(ongoing_peace_prob_month_long,individual = TRUE, "previous month", "peace")
# previous year
create_density_plot(ongoing_peace_prob_year_long,individual = TRUE, "previous year", "peace")

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
create_density_plot(onset_prob_month_long,individual = TRUE, "previous month", "onset")
# previous year
create_density_plot(onset_prob_year_long,individual = TRUE, "previous year", "onset")

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
# previous year
create_density_plot(prev_peace_prob_year_long,individual = TRUE, "previous year", "previous peace")

# ggsave("plots_tobi/density_previous_peace_month.png", 
#        plot = density_previous_peace_month_curves, width = 20, height = 12, dpi = 300, bg = "white")

# ## boxplots -----
# # previous month
# create_box_plot(prev_peace_prob_month_long,individual = TRUE, "previous month", "previous peace")
# # previous year
# create_box_plot(prev_peace_prob_year_long,individual = TRUE, "previous year", "previous peace")



## previous peace density-decomposition-month for predicted prob. > 0! -----
# filter dataset
decomp_prev_peace_prob_month_long <- prev_peace_prob_month_long %>%
  filter(predictive_probability > 0)

# ideal forecast would have everything green for p close to 0 and 
# everything yellow for high probabilities (p close to 1).
density_decomposition_previous_peace_month_curves <- ggplot(data = decomp_prev_peace_prob_month_long, aes(x = predictive_probability)) +
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


# ggsave("plots_tobi/density_decomposition_previous_peace_month.png", 
#        plot = density_decomposition_previous_peace_month_curves, width = 20, height = 12, dpi = 300, bg = "white")

## Proportion of fatality probability exceeding x% (from previous peace for >0%, >1%, >2%, >5%) -----
thresholds <- c(0, 0.01, 0.02, 0.05)

prob_proportion_exceeding_thresh_df <- prev_peace_prob_month_long %>%
  group_by(model) %>%
  summarise(across(predictive_probability, list(
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


## CORP (reliability diagram) -----
prev_peace_prob_month_long_binary_out <- prev_peace_prob_month_long %>%
  mutate(outcome = ifelse(outcome >= 1, 1, 0))

# list to store the CORP
corp_plots_list <- list()

for (model_name in setdiff(names(models_predictive_probabilities), excluded_models)) {
  
  reliability_data <- prev_peace_prob_month_long_binary_out %>%
    filter(model == model_name)
  
  r <- reliabilitydiag(x = reliability_data$predictive_probability, y = reliability_data$outcome)
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

for (model_name in setdiff(names(models_predictive_probabilities), excluded_models)) {
  
  roc_data <- prev_peace_prob_month_long_binary_out %>%
    filter(model == model_name)
  
  mm <- mmdata(roc_data$predictive_probability, roc_data$outcome)
  
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
murphy_data <- prev_peace_prob_month_long_binary_out %>%
  pivot_wider(names_from = model, values_from = predictive_probability) %>%
  select(5:ncol(.))


mr_conflict <- murphy(subset(murphy_data, select = 
                               c(outcome,
                                 bodentien_rueter_negbin,
                                 bodentien_rueter_neuralnet,
                                 conflictforecast_v2,
                                 submission_final_omm,
                                 conflictology,
                                 last,
                                 quantile_forecast,
                                 ShapeFinder,
                                 zero)), 
                      y_var = "outcome")


df_est_conflict <- estimates(mr_conflict)
selected_murphy_plot <- ggplot(df_est_conflict) + 
  geom_path(aes(x = knot, y = mean_score, col = forecast), linewidth = 0.6, alpha = 1) + 
  labs(title = "Murphy diagram of predicted fatality prob. in case of previous peace (01-2018 to 12-2023, all countries)",
       subtitle = "Reference: previous month, individual models")
  

#ggsave("plots_tobi/murphy_previous_peace_month__selected.png", plot = selected_murphy_plot, width = 10, height = 9, dpi = 300)                       



