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
observations_18_23 <- do.call(rbind, lapply(files_actuals_from18, arrow::read_parquet)) # observations from 2018 - 2023

country_id_list <- read_csv(
  "country_list.csv",
  col_types = cols(
    country_id = col_integer(),
    id         = col_integer(),
    name       = col_character()
  )
)

## -----
# load predictive samples from benchmark models and submitted forecasts
## -----

benchmark_names <- c("boot_240", "conflictology", "last", "zero")
submissions_names <- c("bodentien_rueter_negbin", "conflictforecast_v2", "Neg_Bin_GAM", "Neg_Bin_GLMM",
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



# set country ids and observations etc
country_ids <- unique(observations_18_23$country_id)
actuals_ids <- 457:528 # month_ids for 01-2018, 02-2018, ..., 12-2023
cm_pairs <- cbind(rep(country_ids, each = length(actuals_ids)),
                  rep(actuals_ids, length(country_ids)))

date_seq <- seq(
  from = as.Date("2018-01-01"),
  to   = as.Date("2023-12-01"),
  by   = "month"
)

# format as "MM-YYYY"
month_labels <- format(date_seq, "%m-%Y")

# Named character vector: names are month_ids
month_lookup_vec <- setNames(month_labels, actuals_ids)


### CRPS -BRIER - LOG Score Calculation -------------------------------------------------------------------------

## -----
# compute crps on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----

models_crps <- list()


# for (m in 1:n_models) {
#   crps_m <- apply(cm_pairs, 1, function(cm_pair) {
#     print(paste0("Benchmark/model (", m, "/", n_models, "): ", model_names[m], ", country: ", cm_pair[1], ", month: ", cm_pair[2]))
#     true_observation <- observations_18_23 %>%
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
    left_join(observations_18_23 %>%
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
## REMOVE excluded models
## -----
# exclude some models for all plots and analysis
excluded_models <- c("Neg_Bin_GAM", "P_GAM", "TW_GAM")

# delete models
models_scoring_rules <- models_scoring_rules[
  ! names(models_scoring_rules) %in% excluded_models
]

model_names <- names(models_scoring_rules)
n_models <- length(model_names)





### CRPS -BRIER - LOG Score Plot Data Preparation -------------------------------------------------------------------------

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

observations_17_23 <- rbind(observations_17[,colnames(observations_18_23)], observations_18_23) %>%
  arrange(country_id, month_id)

actual_conflict <- observations_17_23 %>%
  filter(country_id %in% country_ids) %>%
  filter(month_id %in% actuals_ids) %>%
  mutate(conflict = outcome>0)

prev_conflict <- list_cbind(lapply(1:12, function(prev_month) {
  observations_17_23 %>%
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
## Plot Data: CRPS values by conflict situation
## -----
models_crps_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "crps")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_crps_conflict) <- c("country_id", "month_id", model_names)
models_crps_conflict <- list(models_crps_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

models_crps_conflict_month <- models_crps_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum)

models_crps_conflict_month[,2:ncol(models_crps_conflict_month)] <- models_crps_conflict_month[,2:ncol(models_crps_conflict_month)] / nrow(models_crps_conflict) # compute contributions to average CRPS


# create ggplot data frames
crps_month <- data.frame("CRPS" = unlist(c(models_crps_conflict_month[,2:ncol(models_crps_conflict_month)])),
                         "Situation" = rep(models_crps_conflict_month$situation_month, ncol(models_crps_conflict_month)-1),
                         "Model" = rep(names(models_crps_conflict_month)[2:ncol(models_crps_conflict_month)], each = 4))

# print CRPS contribution of peace months
print(colSums(models_crps_conflict_month[4,2:ncol(models_crps_conflict_month)] ))

# print overall CRPS per model
print(colSums(models_crps_conflict_month[1:nrow(models_crps_conflict_month),2:ncol(models_crps_conflict_month)] ))


## -----
## Plot Data: Brier values by conflict situation
## -----
models_brier_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "brier_onset")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_brier_conflict) <- c("country_id", "month_id", model_names)
models_brier_conflict <- list(models_brier_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

## -----
# a) create ggplots of contributions to average Brier scores for all conflict situations
## -----
models_brier_conflict_month <- models_brier_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum)

models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] <- models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] / nrow(models_brier_conflict) # compute contributions to average brier

brier_month <- data.frame("Brier" = unlist(c(models_brier_conflict_month[,2:ncol(models_brier_conflict_month)])),
                          "Situation" = rep(models_brier_conflict_month$situation_month, ncol(models_brier_conflict_month)-1),
                          "Model" = rep(names(models_brier_conflict_month)[2:ncol(models_brier_conflict_month)], each = 4))



## ---
## Plot-Data: log-score values by conflict situation for the onset problem (y \in {0,1})
## ---
models_logscore_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "log_score_eps_onset")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_logscore_conflict) <- c("country_id", "month_id", model_names)
models_logscore_conflict <- list(models_logscore_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

models_logscore_conflict_month <- models_logscore_conflict %>% select(!c("country_id", "month_id", "situation_year")) %>% group_by(situation_month) %>% summarise_all(sum)

models_logscore_conflict_month[,2:ncol(models_logscore_conflict_month)] <- models_logscore_conflict_month[,2:ncol(models_logscore_conflict_month)] / nrow(models_logscore_conflict) # compute contributions to average CRPS


## ---
## Data for Onset Prediction
## ---

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
for (model_name in names(models_scoring_rules)) {
  
  joined_data <- models_scoring_rules[[model_name]] %>%
    inner_join(conflict_situations, by = c("month_id", "country_id")) %>%
    mutate(model = model_name) # add model name
  
  # add to dataframe
  prob_models_all_situation_long <- bind_rows(prob_models_all_situation_long, joined_data)
}

# filter for "onset" or "peace" month
prev_peace_prob_month_long <- prob_models_all_situation_long %>%
  filter(situation_month == "peace" | situation_month == "onset")

# onset actuals
prev_peace_prob_month_long_binary_actual <- prev_peace_prob_month_long %>%
  mutate(actual = ifelse(actual >= 1, 1, 0))




################################################################################
## PLOTS
################################################################################
model_labels <- c(
  "zero"                         = "VIEWS Zero",
  "last"                         = "VIEWS Last",
  "conflictology"                = "VIEWS Conflictology",
  "conflictforecast_v2"          = "CFLT RF",
  "submission_muchlinski_thornhill" = "MT ZeroInfl GAM",
  "submission_final_omm"         = "RV O MM",
  "submission_final_gpcmm"       = "RV GPC MM",
  "submission_final_hpmm"        = "RV HP MM",
  "unito_transformer"            = "UNITO NB Transformer",
  "Neg_Bin_GLMM"                 = "BDT NB GLMM",
  "tft"                          = "CCEW TFT",
  "ShapeFinder"                  = "PACE ShapeFinder",
  "quantile_forecast"            = "DB Quantile",
  "boot_240"                     = "VIEWS Bootstrap",
  "bodentien_rueter_negbin"      = "BR NB",
  "P_GLMM"                       = "BDT P GLMM",
  "TW_GLMM"                      = "BDT TW GLMM" 
  
)


## -----
## Selected models
## -----

## --
## change this if needed. everything else is dynamic!
selected_models <- c("zero","last",
                     "conflictology","conflictforecast_v2",
                     "submission_muchlinski_thornhill", "submission_final_omm",
                     "unito_transformer", "P_GLMM")
##
## --



selected_colors <- c(
  "#009682",  # green  
  "#4664aa",  # blue
  "#23a1e0",  # maygreen
  "black", # grey
  "#df9b1b",  # orange
  "#8cb63c",  # yellow
  "#a22223",  # red
  "#a3107c"   # purple
)

selected_model_colors <- setNames(selected_colors, selected_models)

selected_model_labels <- model_labels[selected_models]

## ---
## Remaining models
## ---

remaining_models <- setdiff(model_names, selected_models)

remaining_colors <- c(
  "#004B41",  # von "#009682"
    "#233357",  # von "#4664aa"
  "#19719D",  # von "#23a1e0"
  "#4D4D4D",  # von "black"
  "#6F4D0D",  # von "#df9b1b"
  "#46631E",  # von "#8cb63c"
  "#511111",  # von "#a22223"
  "#51103C",   # von "#a3107c"
  "#2B1A05"    # additional
)

remaining_model_colors <- setNames(remaining_colors, remaining_models)

remaining_model_labels <- model_labels[remaining_models]

model_colors <- c(
  selected_model_colors,
  remaining_model_colors
)

theme_setup <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
  axis.title = ggplot2::element_text(size = 13),
  legend.title = ggplot2::element_text(size = 13),
  axis.text = ggplot2::element_text(size = 12),
  legend.text = ggplot2::element_text(size = 11)
)


## -----------------------------------------------------------------------------
## Number of total fatalities worldwide per month, stacked bar plot
## -----------------------------------------------------------------------------

# 3. Top-Länder bestimmen 
top_country_fatalities <- observations_18_23 %>%
  group_by(country_id) %>%
  summarise(
    total_fatalities = sum(outcome, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  slice_max(total_fatalities, n = 5) 

top_countries <- top_country_fatalities %>%
  pull(country_id)


fatalities_worldwide_plot_data <- observations_18_23 %>%
  mutate(
    # if_else liefert immer einen Vektor gleicher Länge
    country_category = if_else(
      country_id %in% top_countries,
      # TRUE-Zweig: ID als Zeichenkette, damit "others" passt
      as.character(country_id),
      # FALSE-Zweig
      "others"
    )
  ) %>%
  group_by(month_id, country_category) %>%
  summarise(
    n_fatalities = sum(outcome, na.rm = TRUE),
    .groups      = "drop"
  )

# change country_id's to names
fatalities_worldwide_plot_data <- fatalities_worldwide_plot_data %>%
  mutate(
    country_category = if_else(
      country_category == "others",
      "others",
      # match liefert für jede ID die Position in country_id_list
      country_id_list$name[
        match(
          as.character(country_category),
          as.character(country_id_list$country_id)
        )
      ]
    )
  )


fatalities_months <- sort(unique(fatalities_worldwide_plot_data$month_id))

# to determin maximum y range
fatalities_worldwide_y_range <- fatalities_worldwide_plot_data %>%
  group_by(month_id) %>%
  summarise(
    total_fatalities = sum(n_fatalities)
  )
#max(fatalities_worldwide_y_range$total_fatalities)


fatalities_worldwide_plot <- ggplot(fatalities_worldwide_plot_data, aes(fill=country_category, y=n_fatalities, x=month_id)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(
    values = c("Ukraine" = "#E69F00", "Yemen" = "#0072B2",
               "Afghanistan"="#D55E00","Syria" = "#51103C",
               "Ethiopia" = "grey30", "others" = "#009E73")
  ) +
  labs(title = "Aggregated Fatalities for Top 5 Countries and All Others",
       x = "Month") + 
  scale_y_continuous(
    "Fatalities",
    breaks = seq(0, 150000, length.out = 6),
    labels = function(x) {
      # format() fügt Tausender‐Punkte und Komma‐Dezimal an
      format(x,
             big.mark    = ".",
             decimal.mark= ",",
             scientific  = FALSE,
             trim        = TRUE)
    }
  ) + 
  scale_x_continuous(
    breaks = seq(min(fatalities_months), max(fatalities_months), by = 6),
    labels = month_lookup_vec[as.character(seq(min(fatalities_months), max(fatalities_months), by = 6))]
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.ticks       = element_line(color = "black")
  )

ggsave("final_plots/fatalities_worldwide.png",
       plot = fatalities_worldwide_plot, width = 1.2 * 2822, height = 1.2 * 1322, dpi = 300, units = "px",
       bg="white")


## -----------------------------------------------------------------------------
## Number of Countries experiencing > 0 fatalities per month for the test window
## -----------------------------------------------------------------------------
greater_zero_fatalities_plot_data <- observations_18_23 %>%
  filter(outcome > 0) %>%
  group_by(month_id) %>%
  summarise(
    n_countries = n_distinct(country_id),
    .groups     = "drop"
  )

greater_zero_fatalities_plot <- ggplot(greater_zero_fatalities_plot_data, aes(y = n_countries, x = month_id)) +
  geom_bar(position="stack", stat="identity",fill  = "orange3", alpha = 0.7, width = 0.7, color = "black") +
  labs(title = "Countries Experiencing Conflict Fatalities by Month",
       x = "Month") + 
  scale_x_continuous(
    breaks = seq(min(fatalities_months), max(fatalities_months), by = 6),
    labels = month_lookup_vec[as.character(seq(min(fatalities_months), max(fatalities_months), by = 6))]
  ) +
  scale_y_continuous(
    "Countries"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.ticks       = element_line(color = "black")
  )

ggsave("final_plots/greater_zero_fatalities.png",
       plot = greater_zero_fatalities_plot, width = 1.2 * 2822, height = 1.2 * 1322, dpi = 300, units = "px",
       bg="white")

## -----------------------------------------------------------------------------
## Number of Onset events per month for the test window
## -----------------------------------------------------------------------------
conflict_situations <- conflict_situations %>%
  select(-situation_year)
  
onsets_per_month_plot_data <- conflict_situations %>%
  group_by(month_id) %>%
  summarise(
    n_onset = sum(situation_month == "onset", na.rm = TRUE),
    .groups  = "drop"
  )


onsets_per_month_plot <- ggplot(onsets_per_month_plot_data, aes(y = n_onset, x = month_id)) + 
  geom_bar(position="stack", stat="identity",fill  = "orange3", alpha = 0.7, width = 0.7, color = "black") +
  labs(title = "Onset Events by Month for all Countries",
       x = "Month") + 
  scale_x_continuous(
    breaks = seq(min(fatalities_months), max(fatalities_months), by = 6),
    labels = month_lookup_vec[as.character(seq(min(fatalities_months), max(fatalities_months), by = 6))]
  ) +
  scale_y_continuous(
    "Onset Events",
    breaks = seq(0,10, by = 2)
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.ticks       = element_line(color = "black")
  )

ggsave("final_plots/onsets_per_month.png",
       plot = onsets_per_month_plot, width = 1.2 * 2822, height = 1.2 * 1322, dpi = 300, units = "px",
       bg="white")

## -----------------------------------------------------------------------------
## Maybe: CRPS per Country: highlighting 5-10 most important ones (and other) 18 - 23
## -----------------------------------------------------------------------------

average_CRPSscores_over_all_models_18_23 <- prob_models_all_situation_long %>%
  select(
    month_id,
    country_id,
    crps
  ) %>%
  group_by(
    month_id,
    country_id
  ) %>%
  summarise(
    mean_crps = mean(crps, na.rm = TRUE),
    .groups   = "drop"
  )


# 3. Top-Länder bestimmen 
top_country_crps <- average_CRPSscores_over_all_models_18_23 %>%
  group_by(country_id) %>%
  summarise(
    total_crps = sum(mean_crps, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  slice_max(total_crps, n = 5) 

top_countries_crps <- top_country_crps %>%
  pull(country_id)


crps_worldwide_plot_data <- average_CRPSscores_over_all_models_18_23 %>%
  mutate(
    # if_else liefert immer einen Vektor gleicher Länge
    country_category = if_else(
      country_id %in% top_countries_crps,
      # TRUE-Zweig: ID als Zeichenkette, damit "others" passt
      as.character(country_id),
      # FALSE-Zweig
      "others"
    )
  ) %>%
  group_by(month_id, country_category) %>%
  summarise(
    sum_avg_crps = sum(mean_crps, na.rm = TRUE),
    .groups      = "drop"
  )

# change country_id's to names
crps_worldwide_plot_data <- crps_worldwide_plot_data %>%
  mutate(
    country_category = if_else(
      country_category == "others",
      "others",
      # match liefert für jede ID die Position in country_id_list
      country_id_list$name[
        match(
          as.character(country_category),
          as.character(country_id_list$country_id)
        )
      ]
    )
  )

crps_months <- sort(unique(crps_worldwide_plot_data$month_id))

# to determin maximum y range
crps_worldwide_y_range <- crps_worldwide_plot_data %>%
  group_by(month_id) %>%
  summarise(
    avg_crps = sum(sum_avg_crps)
  )
# max(crps_worldwide_y_range$avg_crps)


crps_worldwide_plot <- ggplot(crps_worldwide_plot_data, aes(fill=country_category, y=sum_avg_crps, x=month_id)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(
    values = c("Ukraine" = "#E69F00", "Yemen" = "#0072B2",
               "Afghanistan"="#D55E00","Israel" = "#51103C",
               "Ethiopia" = "grey30", "others" = "#009E73")
  ) +
  labs(title = "Mean CRPS Across Models for Top 5 Countries and All Others",
       x = "Month") + 
  scale_y_continuous(
    "Mean CRPS",
    breaks = seq(0, 150000, length.out = 6),
    labels = function(x) {
      # format() fügt Tausender‐Punkte und Komma‐Dezimal an
      format(x,
             big.mark    = ".",
             decimal.mark= ",",
             scientific  = FALSE,
             trim        = TRUE)
    }
  ) + 
  scale_x_continuous(
    breaks = seq(min(crps_months), max(crps_months), by = 6),
    labels = month_lookup_vec[as.character(seq(min(crps_months), max(crps_months), by = 6))]
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.ticks       = element_line(color = "black")
  )

ggsave("final_plots/crps_worldwide.png",
       plot = crps_worldwide_plot, width = 1.2 * 2822, height = 1.2 * 1322, dpi = 300, units = "px",
       bg="white")


## -----------------------------------------------------------------------------
## Maybe: BRIER for onset cases (previous peace) per Country: 
## highlighting 5-10 most important ones (and other) 18 - 23
## -----------------------------------------------------------------------------
average_BRIERscores_over_all_models_18_23 <- prev_peace_prob_month_long %>%
  select(
    month_id,
    country_id,
    brier_onset
  ) %>%
  group_by(
    month_id,
    country_id
  ) %>%
  summarise(
    mean_brier = mean(brier_onset, na.rm = TRUE),
    .groups   = "drop"
  )


# 3. Top-Länder bestimmen 
top_country_brier <- average_BRIERscores_over_all_models_18_23 %>%
  group_by(country_id) %>%
  summarise(
    total_brier = sum(mean_brier, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  slice_max(total_brier, n = 5) 

top_countries_brier <- top_country_brier %>%
  pull(country_id)


brier_worldwide_plot_data <- average_BRIERscores_over_all_models_18_23 %>%
  mutate(
    # if_else liefert immer einen Vektor gleicher Länge
    country_category = if_else(
      country_id %in% top_countries_brier,
      # TRUE-Zweig: ID als Zeichenkette, damit "others" passt
      as.character(country_id),
      # FALSE-Zweig
      "others"
    )
  ) %>%
  group_by(month_id, country_category) %>%
  summarise(
    sum_avg_brier = sum(mean_brier, na.rm = TRUE),
    .groups      = "drop"
  )

# change country_id's to names
brier_worldwide_plot_data <- brier_worldwide_plot_data %>%
  mutate(
    country_category = if_else(
      country_category == "others",
      "others",
      # match liefert für jede ID die Position in country_id_list
      country_id_list$name[
        match(
          as.character(country_category),
          as.character(country_id_list$country_id)
        )
      ]
    )
  )

brier_months <- sort(unique(brier_worldwide_plot_data$month_id))


brier_worldwide_plot <- ggplot(brier_worldwide_plot_data, aes(fill=country_category, y=sum_avg_brier, x=month_id)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(
    values = c("Libya" = "#E69F00", "Bangladesh" = "#0072B2",
               "Burundi"="#D55E00","Peru" = "#51103C",
               "Rwanda" = "grey30", "others" = "#009E73")
  ) +
  labs(title = "Mean BRIER-Score (Onset) Across Models for Top 5 Countries and All Others",
       x = "Month") + 
  scale_y_continuous(
    "Mean BRIER Score",
    labels = function(x) {
      # format() fügt Tausender‐Punkte und Komma‐Dezimal an
      format(x,
             big.mark    = ".",
             decimal.mark= ",",
             scientific  = FALSE,
             trim        = TRUE)
    }
  ) + 
  scale_x_continuous(
    breaks = seq(min(brier_months), max(brier_months), by = 6),
    labels = month_lookup_vec[as.character(seq(min(brier_months), max(brier_months), by = 6))]
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 13),
    axis.ticks       = element_line(color = "black")
  )

ggsave("final_plots/brier_worldwide.png",
       plot = brier_worldwide_plot, width = 1.2 * 2822, height = 1.2 * 1322, dpi = 300, units = "px",
       bg="white")


## -----------------------------------------------------------------------------
## Selected Onset Events
## -----------------------------------------------------------------------------

onset_event_plot <- function(data, selected_model_bool, country){
  
  onset_bars <- data %>%
    filter(situation_month == "onset") %>%
    select(month_id, actual) %>%
    distinct()
  
  max_onset <- max(onset_bars$actual)
  
  months <- sort(unique(data$month_id))
  
  if(selected_model_bool == TRUE){
    scale_color_own <- scale_color_manual(
      values = selected_model_colors,
      breaks = names(selected_model_labels),
      labels = selected_model_labels
    )
  } else{
    scale_color_own <- scale_color_manual(
      values = remaining_model_colors,
      breaks = names(remaining_model_labels),
      labels = remaining_model_labels
    )
    
  }
  
  title_string <- paste0("Estimated Onset Probabilities – ", country, " (Subset)")
  
  p <- ggplot(data, 
              aes(x = month_id, y = onset_prob_pred, 
                  group = model, color = model)) +
    geom_line(size = 0.7, alpha = 0.5) +
    geom_point(size = 2,    # Punktgröße
               shape = 21,  # gefüllter Kreis
               fill = "white",
               stroke = 0.8) + 
    geom_col(
      data       = onset_bars,
      inherit.aes = FALSE,
      aes(
        x = month_id,
        y = actual / max_onset
      ),
      fill  = "grey40",
      width = 0.2,
      alpha = 0.9,
      linewidth = 0.6
    ) +
    labs(title = title_string,
         x = "Month") + 
    theme_minimal() +
    scale_y_continuous(
      "Predicted Onset Probability", 
      sec.axis = sec_axis(~ . * max_onset, name = "Onset Magnitude")
    ) +
    scale_x_continuous(
      breaks = seq(min(months), max(months), by = 1),
      labels = month_lookup_vec[as.character(months)],
      expand = expansion(add = c(0.3, 0.3))
    ) +
    scale_color_own +
    theme(
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 13),
      axis.ticks       = element_line(color = "black"),
      axis.line.y.right = element_line(color = "grey40"),
      axis.ticks.y.right = element_line(color = "grey40"),
      axis.text.y.right = element_text(color = "grey40"), 
      axis.title.y.right = element_text(color = "grey40")
    )
  return(p)
  
}





## ---
## In-depth: 8 models
## ---

##### onsets per year > 1
# data_onset_events_month_notyear <- prev_peace_prob_month_long %>%
#   filter(situation_month == "onset" &
#            situation_year == "conflict")

## Multiple Small Onsets in a short period of time
# CAR: country 70 
# ab onset monat 467 bis 477 interessant 
# hat in test periode immer wieder ausbrüche von sehr 
# niedriger höhe ~2 und einen höheren 46
CAR_onset_predictions_dataset_selected <- prob_models_all_situation_long %>%
  filter(country_id == 70 & 
           month_id >= 466 &
           month_id <= 477 &
           model %in% selected_models)


car_onset_plot_selected <- onset_event_plot(CAR_onset_predictions_dataset_selected, TRUE, "Central African Republic")

ggsave("final_plots/car_onset_selected.png",
       plot = car_onset_plot_selected, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")


## highest number of fatalities of any onset event
# Ethopia: country 57 
# onset month 491 mit max überhaupt von 1474 aber 
# auch onset month 485 mit 16 fatalities 
# (modelle täuschen sich stark: großer CRPS)
Ethiopia_onset_predictions_dataset_selected <- prob_models_all_situation_long %>%
  filter(country_id == 57 & 
           month_id >= 481 &
           month_id <= 492 &
           model %in% selected_models)


ethiopia_onset_plot_selected <- onset_event_plot(Ethiopia_onset_predictions_dataset_selected, TRUE, "Ethiopia")

ggsave("final_plots/ethiopia_onset_selected.png",
       plot = ethiopia_onset_plot_selected, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")


##### onsets per year = 1
# data_onset_events_year <- prev_peace_prob_month_long %>%
#   filter(situation_month == "onset" &
#            situation_year == "onset")

## Israel: country 218 
# in month 497 maximum aus onset onset mit 252
Israel_onset_predictions_dataset_selected <- prob_models_all_situation_long %>%
  filter(country_id == 218 & 
           month_id >= 486 &
           month_id <= 497 &
           model %in% selected_models)


israel_onset_plot_selected <- onset_event_plot(Israel_onset_predictions_dataset_selected, TRUE, "Israel")

ggsave("final_plots/israel_onset_selected.png",
       plot = israel_onset_plot_selected, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")

##### onsets per year = 0
#Poland 101
# month 457 - 468
Poland_onset_predictions_dataset_selected <- prob_models_all_situation_long %>%
  filter(country_id == 101 & 
           month_id >= 517 &
           month_id <= 528 &
           model %in% selected_models)


poland_onset_plot_selected <- onset_event_plot(Poland_onset_predictions_dataset_selected, TRUE, "Poland")

ggsave("final_plots/poland_onset_selected.png",
       plot = poland_onset_plot_selected, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")




## ---
## Remaining: 9 models
## ---

##### onsets per year > 1
# data_onset_events_month_notyear <- prev_peace_prob_month_long %>%
#   filter(situation_month == "onset" &
#            situation_year == "conflict")

## Multiple Small Onsets in a short period of time
# CAR: country 70 
# ab onset monat 467 bis 477 interessant 
# hat in test periode immer wieder ausbrüche von sehr 
# niedriger höhe ~2 und einen höheren 46
CAR_onset_predictions_dataset_remaining <- prob_models_all_situation_long %>%
  filter(country_id == 70 & 
           month_id >= 466 &
           month_id <= 477 &
           model %in% remaining_models)


car_onset_plot_remaining <- onset_event_plot(CAR_onset_predictions_dataset_remaining, FALSE, "Central African Republic")

ggsave("final_plots/car_onset_remaining.png",
       plot = car_onset_plot_remaining, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")


## highest number of fatalities of any onset event
# Ethopia: country 57 
# onset month 491 mit max überhaupt von 1474 aber 
# auch onset month 485 mit 16 fatalities 
# (modelle täuschen sich stark: großer CRPS)
Ethiopia_onset_predictions_dataset_remaining <- prob_models_all_situation_long %>%
  filter(country_id == 57 & 
           month_id >= 481 &
           month_id <= 492 &
           model %in% remaining_models)


ethiopia_onset_plot_remaining <- onset_event_plot(Ethiopia_onset_predictions_dataset_remaining, FALSE, "Ethiopia")

ggsave("final_plots/ethiopia_onset_remaining.png",
       plot = ethiopia_onset_plot_remaining, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")


##### onsets per year = 1
# data_onset_events_year <- prev_peace_prob_month_long %>%
#   filter(situation_month == "onset" &
#            situation_year == "onset")

## Israel: country 218 
# in month 497 maximum aus onset onset mit 252
Israel_onset_predictions_dataset_remaining <- prob_models_all_situation_long %>%
  filter(country_id == 218 & 
           month_id >= 486 &
           month_id <= 497 &
           model %in% remaining_models)


israel_onset_plot_remaining <- onset_event_plot(Israel_onset_predictions_dataset_remaining, FALSE, "Israel")

ggsave("final_plots/israel_onset_remaining.png",
       plot = israel_onset_plot_remaining, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")

##### onsets per year = 0
#Poland 101
# month 457 - 468
Poland_onset_predictions_dataset_remaining <- prob_models_all_situation_long %>%
  filter(country_id == 101 & 
           month_id >= 517 &
           month_id <= 528 &
           model %in% remaining_models)


poland_onset_plot_remaining <- onset_event_plot(Poland_onset_predictions_dataset_remaining, FALSE, "Poland")


ggsave("final_plots/poland_onset_remaining.png",
       plot = poland_onset_plot_remaining, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
       bg="white")


## -----------------------------------------------------------------------------
## distribution of onset probabilities >= 0
## -----------------------------------------------------------------------------

# ideal forecast would have everything green for p close to 0 and 
# everything yellow for high probabilities (p close to 1).

# kernel density estimation parameter
adjust_value <- 0.4


# used for >= and > 0
create_distribution_plot <- function(data, model_name){
  
  filter_data <- data %>%
    filter(model == model_name)
  
  p <- ggplot(data = filter_data, aes(x = onset_prob_pred)) +
    # "peace" and "onset" density (i.e. "previous peace")
    geom_density(aes(y = after_stat(density), fill = "previous peace"), 
                 alpha = 0.9, 
                 size = 0.8,
                 adjust = adjust_value) +
    
    # "onset" density weighted relative to overall dist of prev_peace_prob_month_long
    geom_density(data = filter(filter_data, situation_month == "onset"),
                 aes(y = after_stat(density) * nrow(filter(filter_data, situation_month == "onset")) / nrow(filter_data),
                     fill = "onset"),
                 alpha = 0.8,
                 size = 0.8,
                 adjust = adjust_value) +
    
    scale_y_sqrt() +
    
    scale_fill_manual(name="conflict situation",
                      values=c("previous peace"="#556B2F","onset"="#df9b1b")) + 
    theme_minimal_grid() +
    labs( 
      title = model_labels[model_name]
    ) +
    coord_fixed() +
    # x-Achse: Breaks bei 0, .25, .5, .75, 1
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25)
    ) +
    theme(
      panel.background = element_blank(),
      plot.title = element_text(size = 15, hjust = 0.5, face = "plain"),
      legend.position = "none",
      aspect.ratio = 1,
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  
  if(model_name == "zero"){
    
    
    p <- p +
      annotate(
        "rect",
        xmin   = 0,
        ymin   = 0,
        xmax   = 0.1,
        ymax   = 2,
        fill   = "transparent",
        colour = "black",
        size   = 0.7
      )
  }
  
  if(model_name == "submission_final_gpcmm"){
    
    
    p <- p +
      annotate(
        "rect",
        xmin   = 0.9,
        ymin   = 0,
        xmax   = 1,
        ymax   = 2,
        fill   = "transparent",
        colour = "black",
        size   = 0.7
      )
  }
  
  
  
  
  return(p)
  
}


#create common x and y labels
# used for >= and > 0
y_density.grob <- textGrob("Density", 
                           gp=gpar(col="black", fontsize=15), rot=90)

x_density.grob <- textGrob("Forecast value", 
                           gp=gpar(col="black", fontsize=15))


## ---
## In-depth: 8 models
## ---
# list to store the density plots
density_greqZero_list_selected <- list()

for (model_name in selected_models) {
  
  plot <- create_distribution_plot(prev_peace_prob_month_long, model_name)
  
  density_greqZero_list_selected[[model_name]] <- list(plot = plot)
  
}

# modify plot list
density_greqZero_grid_modified_list_selected <- density_greqZero_list_selected

density_greqZero_model_names_selected <- names(density_greqZero_grid_modified_list_selected)


for(i in seq_along(density_greqZero_list_selected)){
  
  model <- density_greqZero_model_names_selected[i]
  
  p <- NULL
  
  if(i == 1){
    p <- density_greqZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i > 1 && i < 5){
    p <- density_greqZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i == 5){
    p <- density_greqZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else {
    p <- density_greqZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  }
  
  
  density_greqZero_grid_modified_list_selected[[model]]$plot <- p
  
}

plots_density_greqZero_gird_modified_selected <- lapply(density_greqZero_grid_modified_list_selected, function(x) x$plot)

# Kombiniere in einem 2x4 Layout
p_density_greqZero_selected <- plot_grid(plotlist = plots_density_greqZero_gird_modified_selected, nrow = 2, ncol = 4, align = "hv", scale = 0.99)


#add to plot
density_greqZero_selected_models_plot <- grid.arrange(arrangeGrob(p_density_greqZero_selected, 
                                                                  left = y_density.grob, 
                                                                  bottom = x_density.grob))

ggsave("final_plots/density_greqZero_selected_models.png",
       plot = density_greqZero_selected_models_plot, width = 1.4 * 3200, height = 1.4 * 1700, dpi = 300, units = "px")

## ---
## Remaining: 9 models
## ---
# list to store the density plots
density_greqZero_list_remaining <- list()

for (model_name in remaining_models) {
  
  plot <- create_distribution_plot(prev_peace_prob_month_long, model_name)
  
  density_greqZero_list_remaining[[model_name]] <- list(plot = plot)
  
}

# modify plot list
density_greqZero_grid_modified_list_remaining <- density_greqZero_list_remaining

density_greqZero_model_names_remaining <- names(density_greqZero_grid_modified_list_remaining)



for(i in seq_along(density_greqZero_list_remaining)){
  
  model <- density_greqZero_model_names_remaining[i]
  
  p <- NULL
  
  if(i == 1 || i == 4){
    p <- density_greqZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if((i > 1 && i < 4) || (i > 4 && i < 7)){
    p <- density_greqZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i == 7){
    p <- density_greqZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else {
    p <- density_greqZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  }
  
  
  density_greqZero_grid_modified_list_remaining[[model]]$plot <- p
  
}


plots_density_greqZero_gird_modified_remaining <- lapply(density_greqZero_grid_modified_list_remaining, function(x) x$plot)

# Kombiniere in einem 2x4 Layout
p_density_greqZero_remaining <- plot_grid(plotlist = plots_density_greqZero_gird_modified_remaining, nrow = 3, ncol = 3, align = "hv", scale = 0.99)

#add to plot
density_greqZero_remaining_models_plot <- grid.arrange(arrangeGrob(p_density_greqZero_remaining, 
                                                                    left = y_density.grob, 
                                                                    bottom = x_density.grob))

ggsave("final_plots/density_greqZero_remaining_models.png",
       plot = density_greqZero_remaining_models_plot, width = 1.5 * 2200, height = 1.5 * 2000, dpi = 300, units = "px")



## -----------------------------------------------------------------------------
## distribution of onset probabilities > 0
## -----------------------------------------------------------------------------
# filter dataset
decomp_prev_peace_prob_month_long <- prev_peace_prob_month_long %>%
  filter(onset_prob_pred > 0)


## ---
## In-depth: 8 models
## ---
# list to store the density plots
density_greatZero_list_selected <- list()

for (model_name in selected_models) {
  
  plot <- create_distribution_plot(decomp_prev_peace_prob_month_long, model_name)
  
  density_greatZero_list_selected[[model_name]] <- list(plot = plot)
  
}

# modify plot list
density_greatZero_grid_modified_list_selected <- density_greatZero_list_selected

density_greatZero_model_names_selected <- names(density_greatZero_grid_modified_list_selected)


for(i in seq_along(density_greatZero_list_selected)){
  
  model <- density_greatZero_model_names_selected[i]
  
  p <- NULL
  
  if(i == 1){
    p <- density_greatZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i > 1 && i < 5){
    p <- density_greatZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i == 5){
    p <- density_greatZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else {
    p <- density_greatZero_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  }
  
  
  density_greatZero_grid_modified_list_selected[[model]]$plot <- p
  
}

plots_density_greatZero_gird_modified_selected <- lapply(density_greatZero_grid_modified_list_selected, function(x) x$plot)

# Kombiniere in einem 2x4 Layout
p_density_greatZero_selected <- plot_grid(plotlist = plots_density_greatZero_gird_modified_selected, nrow = 2, ncol = 4, align = "hv", scale = 0.99)

#create common x and y labels
y_density_greatZero_.grob <- textGrob("Density", 
                                      gp=gpar(col="black", fontsize=15), rot=90)

x_density_greatZero_.grob <- textGrob("Forecast value", 
                                      gp=gpar(col="black", fontsize=15))

#add to plot
density_greatZero_selected_models_plot <- grid.arrange(arrangeGrob(p_density_greatZero_selected, 
                                                                   left = y_density.grob, 
                                                                   bottom = x_density.grob))

ggsave("final_plots/density_greatZero_selected_models.png",
       plot = density_greatZero_selected_models_plot, width = 1.4 * 3200, height = 1.4 * 1700, dpi = 300, units = "px")



## ---
## Remaining: 9 models
## ---
# list to store the density plots
density_greatZero_list_remaining <- list()

for (model_name in remaining_models) {
  
  plot <- create_distribution_plot(decomp_prev_peace_prob_month_long, model_name)
  
  density_greatZero_list_remaining[[model_name]] <- list(plot = plot)
  
}

# modify plot list
density_greatZero_grid_modified_list_remaining <- density_greatZero_list_remaining

density_greatZero_model_names_remaining <- names(density_greatZero_grid_modified_list_remaining)



for(i in seq_along(density_greatZero_list_remaining)){
  
  model <- density_greatZero_model_names_remaining[i]
  
  p <- NULL
  
  if(i == 1 || i == 4){
    p <- density_greatZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if((i > 1 && i < 4) || (i > 4 && i < 7)){
    p <- density_greatZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i == 7){
    p <- density_greatZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else {
    p <- density_greatZero_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  }
  
  
  density_greatZero_grid_modified_list_remaining[[model]]$plot <- p
  
}


plots_density_greatZero_gird_modified_remaining <- lapply(density_greatZero_grid_modified_list_remaining, function(x) x$plot)

# Kombiniere in einem 2x4 Layout
p_density_greatZero_remaining <- plot_grid(plotlist = plots_density_greatZero_gird_modified_remaining, nrow = 3, ncol = 3, align = "hv", scale = 0.99)

#add to plot
density_greatZero_remaining_models_plot <- grid.arrange(arrangeGrob(p_density_greatZero_remaining, 
                                                                   left = y_density.grob, 
                                                                   bottom = x_density.grob))

ggsave("final_plots/density_greatZero_remaining_models.png",
       plot = density_greatZero_remaining_models_plot, width = 1.5 * 2200, height = 1.5 * 2000, dpi = 300, units = "px")




## -----------------------------------------------------------------------------
## CRPS decomposition
## -----------------------------------------------------------------------------

## ---
## One model per team: 13 models
## ---
crps_month_selected_models <- crps_month

excluded_secondary_models <- c("submission_final_gpcmm", "submission_final_hpmm", "TW_GLMM",  "Neg_Bin_GLMM")

crps_month_selected_models <- crps_month_selected_models[-which(crps_month_selected_models$Model %in% excluded_secondary_models), ]

# Rename models
crps_month_selected_models$Model <- recode(
  crps_month_selected_models$Model,
  !!!model_labels
)

# determine color of y-axis model names (black if selected for in-depth analysis)
one_model_per_team_names <- unique(crps_month_selected_models$Model)[order(unique(crps_month_selected_models$Model))]
selected_models_bool <- one_model_per_team_names %in% selected_model_labels 

col_vec <- ifelse(selected_models_bool, "black", "grey40")

# plot
crps_conflict_situation_plot <- ggplot(crps_month_selected_models, 
                                       aes(fill = Situation, y = Model, x = .data[["CRPS"]])) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean CRPS per Conflict Situation",
       x = "Mean CRPS") +
  scale_fill_manual("Situation", values = c("conflict" = "#a22223", "deescalation" = "#009682", "onset" = "#df9b1b", "peace" = "#4664aa"),
                    labels = c("conflict" = "Conflict", "deescalation" = "Deescalation", "onset" = "Onset", "peace" = "Peace")) +
  theme_setup +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(colour = col_vec))

crps_conflict_situation_plot

ggsave("final_plots/crps_conflict_situation.png",
         plot = crps_conflict_situation_plot, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px")


## -----------------------------------------------------------------------------
## murphy diagram
## -----------------------------------------------------------------------------
murphy_data <- prev_peace_prob_month_long_binary_actual %>%
  select(1:6,8) %>%
  pivot_wider(names_from = model, values_from = onset_prob_pred) %>%
  select(5:ncol(.))


## ---
## In-depth: 8 models
## ---
mr_conflict_selected <- murphy(subset(murphy_data, select = c(selected_models,"actual")), 
                      y_var = "actual")

df_est_conflict_selected <- estimates(mr_conflict_selected)

murphy_diagram_selected_models <- ggplot(df_est_conflict_selected) + 
  geom_path(aes(x = knot, y = mean_score, col = forecast), linewidth = 1, alpha = 0.8) + 
  scale_color_manual(
    values = selected_model_colors,
    breaks = names(selected_model_labels),
    labels = selected_model_labels
  ) +
  theme_bw() +
  labs(
    title = "Murphy Diagram",
    x = expression(paste("Parameter ", theta)),
    y = "Mean elementary score",
    color = ""
  ) + 
  theme_setup + 
  theme(legend.position = "none",#"bottom
        aspect.ratio = 1)

murphy_diagram_selected_models


## ---
## Remaining: 9 models
## ---
mr_conflict_remaining <- murphy(subset(murphy_data, select = c(remaining_models,"actual")), 
                               y_var = "actual")

df_est_conflict_remaining <- estimates(mr_conflict_remaining)

murphy_diagram_remaining_models <- ggplot(df_est_conflict_remaining) + 
  geom_path(aes(x = knot, y = mean_score, col = forecast), linewidth = 1, alpha = 0.8) + 
  scale_color_manual(
    values = remaining_model_colors,
    breaks = names(remaining_model_labels),
    labels = remaining_model_labels
  ) +
  theme_bw() +
  labs(
    title = "Murphy Diagram",
    x = expression(paste("Parameter ", theta)),
    y = "Mean elementary score",
    color = ""
  ) + 
  theme_setup + 
  theme(legend.position = "none",#"bottom
        aspect.ratio = 1)

murphy_diagram_remaining_models


## -----------------------------------------------------------------------------
## ROC curves
## -----------------------------------------------------------------------------

## 
## https://cran.r-project.org/web/packages/precrec/vignettes/introduction.html
## 
roc_data <- prev_peace_prob_month_long_binary_actual


## ---
## In-depth: 8 models
## ---
# filter relevant models
roc_data_selected <- roc_data %>% 
  filter(model %in% selected_models)

# create list of scores
score_list_selected <- lapply(selected_models, function(m) {
  roc_data_selected %>%
    filter(model == m) %>%
    pull(onset_prob_pred)
})

# create list of labels
label_list_selected <- lapply(selected_models, function(m) {
  roc_data_selected %>%
    filter(model == m) %>%
    pull(actual)
})

# join lists
scores_joined_selected <- join_scores(score_list_selected)
labels_joined_selected <- join_labels(label_list_selected)

# roc-curve data
mm_selected_models <- mmdata(
  scores   = scores_joined_selected,
  labels   = labels_joined_selected,
  modnames = selected_models
)


roc_curve_selected_models <- autoplot(evalmod(mm_selected_models), curvetype = "ROC") + 
  ggplot2::geom_line(size = 1, alpha = 0.8) + 
  ggplot2::scale_color_manual(
    values = selected_model_colors,
    breaks = names(selected_model_labels),
    labels = selected_model_labels
  ) +
  ggplot2::labs(
    title = "ROC Curve",
    x = "FAR",
    y = "HR",
    color = ""
  ) +
  theme_setup +
  ggplot2::theme(legend.position = "none") #"bottom

roc_curve_selected_models


## ---
## Remaining: 9 models
## ---
# filter relevant models
roc_data_remaining <- roc_data %>% 
  filter(model %in% remaining_models)

# create list of scores
score_list_remaining <- lapply(remaining_models, function(m) {
  roc_data_remaining %>%
    filter(model == m) %>%
    pull(onset_prob_pred)
})

# create list of labels
label_list_remaining <- lapply(remaining_models, function(m) {
  roc_data_remaining %>%
    filter(model == m) %>%
    pull(actual)
})

# join lists
scores_joined_remaining <- join_scores(score_list_remaining)
labels_joined_remaining <- join_labels(label_list_remaining)

# roc-curve data
mm_remaining_models <- mmdata(
  scores   = scores_joined_remaining,
  labels   = labels_joined_remaining,
  modnames = remaining_models
)


roc_curve_remaining_models <- autoplot(evalmod(mm_remaining_models), curvetype = "ROC") + 
  ggplot2::geom_line(size = 1, alpha = 0.8) + 
  ggplot2::scale_color_manual(
    values = remaining_model_colors,
    breaks = names(remaining_model_labels),
    labels = remaining_model_labels
  ) +
  ggplot2::labs(
    title = "ROC Curve",
    x = "FAR",
    y = "HR",
    color = ""
  ) +
  theme_setup +
  ggplot2::theme(legend.position = "none") #"bottom

roc_curve_remaining_models

## -----------------------------------------------------------------------------
## Reliability diagram
## -----------------------------------------------------------------------------

## function to create custom reliability diagram
create_reliability_diag <- function(data, forecast_model) {
  
  reliability_data_selected <- data %>%
    filter(model == forecast_model)
  
  r_selected <- reliabilitydiag(x = reliability_data_selected$onset_prob_pred, y = reliability_data_selected$actual)
  
  
  
  # strip the geomSegment out (red and overwritten)
  reliability_plot <- autoplot(r_selected)
  
  # is_seg <- sapply(reliability_plot$layers, function(layer) {
  #   inherits(layer$geom, "GeomSegment")
  # })
  # print(is_seg[5])
  
  is_seg <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  
  # 4) strip them out
  reliability_plot$layers <- reliability_plot$layers[!is_seg]
  
  
  
  
  
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
  point_draw <- NULL
  
  if(length(df_clean$CEP) > 1){
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
  } else{
    point_draw <- ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$CEP,
      ),
      data = data_estim,
      colour = path_color,
      shape = 19,
      size = 2 #3
    )
  }
  
  
  
  p <- reliability_plot + ggplot2::labs( #autoplot(r_selected) + ggplot2::labs(
    title = model_labels[model_name],
    x = "Forecast value",
    y = "CEP") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$x, 
        y = .data$CEP, 
        xend = .data$x_end, 
        yend = .data$CEP_end),
      data = df_segments,
      linewidth = 1.4, #1.8
      colour = path_color
    ) +
    ggplot2::geom_path(
      mapping = ggplot2::aes(
        x = .data$x,
        y = .data$CEP
      ),
      data = data_estim,
      linewidth = 1, #1
      colour = path_color
    ) + point_draw
  
  return(p)
  
}




## ---
## In-depth: 8 models
## ---
# list to store the CORP
corp_plots_list_selected <- list()

for (model_name in selected_models) {
  
  plot <- create_reliability_diag(prev_peace_prob_month_long_binary_actual, model_name)
  
  corp_plots_list_selected[[model_name]] <- list(plot = plot)
  
}

# modify plot list
corp_plots_grid_modified_list_selected <- corp_plots_list_selected

corp_plot_model_names_selected <- names(corp_plots_grid_modified_list_selected)


for(i in seq_along(corp_plots_list_selected)){
  
  model <- corp_plot_model_names_selected[i]
  
  p <- NULL
  
  if(i == 1){
    p <- corp_plots_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i > 1 && i < 5){
    p <- corp_plots_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i == 5){
    p <- corp_plots_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else {
    p <- corp_plots_grid_modified_list_selected[[model]]$plot + 
      ggplot2::theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  }
  
  
  corp_plots_grid_modified_list_selected[[model]]$plot <- p
  
}

plots_corp_gird_modified_selected <- lapply(corp_plots_grid_modified_list_selected, function(x) x$plot)

# Kombiniere in einem 2x4 Layout
p_all_selected <- plot_grid(plotlist = plots_corp_gird_modified_selected, nrow = 2, ncol = 4, align = "hv", scale = 0.99)

#create common x and y labels
y.grob <- textGrob("CEP", 
                   gp=gpar(col="black", fontsize=15), rot=90)

x.grob <- textGrob("Forecast value", 
                   gp=gpar(col="black", fontsize=15))

#add to plot
reliability_diag_selected_models_plot <- grid.arrange(arrangeGrob(p_all_selected, left = y.grob, bottom = x.grob))

## ---
## Remaining: 9 models
## ---
# list to store the CORP
corp_plots_list_remaining <- list()

for (model_name in remaining_models) {
  
  plot <- create_reliability_diag(prev_peace_prob_month_long_binary_actual, model_name)
  
  corp_plots_list_remaining[[model_name]] <- list(plot = plot)
  
}

# modify plot list
corp_plots_grid_modified_list_remaining <- corp_plots_list_remaining

corp_plot_model_names_remaining <- names(corp_plots_grid_modified_list_remaining)


for(i in seq_along(corp_plots_list_remaining)){
  
  model <- corp_plot_model_names_remaining[i]
  
  p <- NULL
  
  if(i == 1 || i == 4){
    p <- corp_plots_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if((i > 1 && i < 4) || (i > 4 && i < 7)){
    p <- corp_plots_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else if(i == 7){
    p <- corp_plots_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  } else {
    p <- corp_plots_grid_modified_list_remaining[[model]]$plot + 
      ggplot2::theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     plot.margin = unit(c(0,0,0,0), "cm"))
  }
  
  
  corp_plots_grid_modified_list_remaining[[model]]$plot <- p
  
}

plots_corp_gird_modified_remaining <- lapply(corp_plots_grid_modified_list_remaining, function(x) x$plot)

# Kombiniere in einem 2x4 Layout
p_all_remaining <- plot_grid(plotlist = plots_corp_gird_modified_remaining, nrow = 3, ncol = 3, align = "hv", scale = 0.99)

#create common x and y labels
y_remaining.grob <- textGrob("CEP", 
                   gp=gpar(col="black", fontsize=15), rot=90)

x_remaining.grob <- textGrob("Forecast value", 
                   gp=gpar(col="black", fontsize=15))

#add to plot
reliability_diag_remaining_models_plot <- grid.arrange(arrangeGrob(p_all_remaining, left = y_remaining.grob, bottom = x_remaining.grob))



## -----------------------------------------------------------------------------
## Merge Murpy Diagram and ROC curve
## -----------------------------------------------------------------------------
# method get legend in cowplot not longer working so method to get legend from
# https://stackoverflow.com/questions/78163631/r-get-legend-from-cowplot-package-no-longer-work-for-ggplot2-version-3-5-0
get_legend2 <- function(plot, legend = NULL) {
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


## ---
## In-depth: 8 models
## ---

# get legend from roc_curve_plot
## ACHTUNG aktuell bei roc curv keine legend vorhanden!!!
#legend_roc <- get_legend2(roc_curve_selected_models + theme(legend.position="bottom"))

# arrange the three plots in a single row
murphy_roc_selected_models_plot <- plot_grid(murphy_diagram_selected_models + theme(legend.position="none"),
                                             NULL,
                  roc_curve_selected_models + theme(legend.position="none"),
                  align = 'vh',
                  hjust = -1,
                  nrow = 1,
                  rel_widths = c(1, 0, 1)
)


## ---
## Remaining: 9 models
## ---
# arrange the three plots in a single row
murphy_roc_remaining_models_plot <- plot_grid(murphy_diagram_remaining_models + theme(legend.position="none"),
                                             NULL,
                                             roc_curve_remaining_models + theme(legend.position="none"),
                                             align = 'vh',
                                             hjust = -1,
                                             nrow = 1,
                                             rel_widths = c(1, 0, 1)
)

## -----------------------------------------------------------------------------
## Merge Murpy Diagram, ROC curve and Reliability Diagram
## -----------------------------------------------------------------------------

## ---
## In-depth: 8 models
## ---

new_triptych_selected <- grid.arrange(
  murphy_roc_selected_models_plot,
  reliability_diag_selected_models_plot,
  nrow = 2,
  heights = c(1, 0.895)
)

# ggsave("final_plots/tryptich_selected_models.png",
#        plot = new_triptych_selected, width = 1.5 * 2100, height = 1.5 * 2200, dpi = 300, units = "px")

## ---
## Remaining: 9 models
## ---

new_triptych_remaining <- grid.arrange(
  murphy_roc_remaining_models_plot,
  reliability_diag_remaining_models_plot,
  nrow = 2,
  heights = c(1, 1.65)
)

# ggsave("final_plots/tryptich_remaining_models.png",
#        plot = new_triptych_remaining, width = 1.4 * 1875, height = 1.4 * 2890.6, dpi = 300, units = "px")


## -----------------------------------------------------------------------------
## BRIER: MSC-DSC-plots
## -----------------------------------------------------------------------------

## IMPORTANT: approach is not elegnat. possible add version with geom_rect later


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


## ---
## In-depth: 8 models
## ---
# filter for selected models
brier_decomposition_results_selected <- brier_decomposition_results[selected_models]

## important remark for BRIER score:
# UNC = variance of X~Ber(p=E(y)=mean(y)) -> p(1-p)
#p <- models_scoring_rules$boot_240$actual_conflict
#mean(p)*(1-mean(p))

# combine to one dataframe
brier_score_decomposition_selected <- bind_rows(brier_decomposition_results_selected) %>%
  select(model, mean_score, MCB, DSC, UNC) %>%
  rename(MeanScore = mean_score)

brier_score_decomposition_barplot_selected <- brier_score_decomposition_selected %>%
  mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)

# data frame with format for the barchart
brier_score_decomposition_barplot_selected <- brier_score_decomposition_barplot_selected %>%
  arrange(MeanScore) %>%  # sort by MeanScore
  pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap, gap1, gap2),
               names_to = "component",
               values_to = "value") %>%
  mutate(class = case_when(
    component %in% c("MeanScore") ~ "SCORE",
    component %in% c("MCB", "UNC") ~ "MCB_UNC",
    component %in% c("DSC", "score_invisible") ~ "DSC_score", 
    component %in% c("gap") ~ "GAP",
    component %in% c("gap1") ~ "GAPmeanscoreDSC",
    component %in% c("gap2") ~ "GAPdscUNC",
  )) %>%
  select(model, class, component, value)

brier_score_decomposition_barplot_selected <- brier_score_decomposition_barplot_selected %>%
  mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))


# dataset for model orderbased on MeanScore
brier_levs_selected <- brier_score_decomposition_barplot_selected %>%
  filter(component == "MeanScore") %>%
  arrange(desc(value)) %>%
  pull(model)

# data for the plot
brier_score_decomposition_barplot_selected <- brier_score_decomposition_barplot_selected %>%
  mutate(
    class_group = case_when(
      class == "MCB_UNC"   ~ 1,
      class == "DSC_score"       ~ 2,
      class == "SCORE" ~ 3,
      class == "GAP" ~ 4,
      class == "GAPmeanscoreDSC" ~ 5,
      class == "GAPdscUNC" ~ 6
      
    ),  
    model = factor(model, levels = brier_levs_selected, ordered = TRUE)
  )

brier_score_decomposition_barplot_selected$model <- recode(
  brier_score_decomposition_barplot_selected$model,
  !!!model_labels
)


# set order of subbars within group
brier_score_decomposition_barplot_selected$component <- factor(
  brier_score_decomposition_barplot_selected$component,
  levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
)

# set order of sub bars within group
brier_score_decomposition_barplot_selected$class_group <- factor(
  brier_score_decomposition_barplot_selected$class_group,
  levels = c(4,1,6,2,3,5)    
)

# set width of mean score bar
mean_score_bar_selected <- 4.5
# set width of mcb,dcs,unc bars
msc_dsc_unc_bar_selected <- 2

# width column
brier_score_decomposition_barplot_selected <- brier_score_decomposition_barplot_selected %>%
  mutate(
    wdth = case_when(
      # MeanScore
      class_group == 3 ~ mean_score_bar_selected, #4
      # MCB DSC
      class_group %in% c(1, 2) ~ msc_dsc_unc_bar_selected, #1
      # gapswithin groups
      class_group == 6 ~ msc_dsc_unc_bar_selected,
      # gap between groups top
      class_group == 5 ~ 15,
      #gap between groups
      class_group == 4 ~ 8 #10
      
    )
  )


brier_score_decomposition_selected <- ggplot(brier_score_decomposition_barplot_selected, aes(x = class_group, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack", width = brier_score_decomposition_barplot_selected$wdth) +
  facet_grid(model ~ ., switch = "y") +
  coord_flip() +
  scale_fill_manual(
    values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
    ##
    breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
    ##
    name = ""
  ) +
  labs(title = "MSC-DSC Mean Brier Score Decomposition") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 18),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "points"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
    plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = ggplot2::element_text(size = 14)
  )

brier_score_decomposition_selected


# ggsave("final_plots/brier_score_decomposition_selected.png",
#        plot = brier_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px")



## ---
## All models
## ---


# combine to one dataframe
brier_score_decomposition_all <- bind_rows(brier_decomposition_results) %>%
  select(model, mean_score, MCB, DSC, UNC) %>%
  rename(MeanScore = mean_score)

brier_score_decomposition_barplot_all <- brier_score_decomposition_all %>%
  mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)

# data frame with format for the barchart
brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
  arrange(MeanScore) %>%  # sort by MeanScore
  pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap, gap1, gap2),
               names_to = "component",
               values_to = "value") %>%
  mutate(class = case_when(
    component %in% c("MeanScore") ~ "SCORE",
    component %in% c("MCB", "UNC") ~ "MCB_UNC",
    component %in% c("DSC", "score_invisible") ~ "DSC_score", 
    component %in% c("gap") ~ "GAP",
    component %in% c("gap1") ~ "GAPmeanscoreDSC",
    component %in% c("gap2") ~ "GAPdscUNC",
  )) %>%
  select(model, class, component, value)

brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
  mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))


# dataset for model orderbased on MeanScore
brier_levs_all <- brier_score_decomposition_barplot_all %>%
  filter(component == "MeanScore") %>%
  arrange(desc(value)) %>%
  pull(model)

# data for the plot
brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
  mutate(
    class_group = case_when(
      class == "MCB_UNC"   ~ 1,
      class == "DSC_score"       ~ 2,
      class == "SCORE" ~ 3,
      class == "GAP" ~ 4,
      class == "GAPmeanscoreDSC" ~ 5,
      class == "GAPdscUNC" ~ 6
      
    ),  
    model = factor(model, levels = brier_levs_all, ordered = TRUE)
  )

brier_score_decomposition_barplot_all$model <- recode(
  brier_score_decomposition_barplot_all$model,
  !!!model_labels
)


# set order of subbars within group
brier_score_decomposition_barplot_all$component <- factor(
  brier_score_decomposition_barplot_all$component,
  levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
)

# set order of sub bars within group
brier_score_decomposition_barplot_all$class_group <- factor(
  brier_score_decomposition_barplot_all$class_group,
  levels = c(4,1,6,2,3,5)
)

# set width of mean score bar
mean_score_bar_all <- 7
# set width of mcb,dcs,unc bars
msc_dsc_unc_bar_all <- 2

# width column
brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
  mutate(
    wdth = case_when(
      # MeanScore
      class_group == 3 ~ mean_score_bar_all, #4
      # MCB DSC
      class_group %in% c(1, 2) ~ msc_dsc_unc_bar_all,#msc_dsc_unc_bar_all, #1
      # gap between groups top
      class_group == 5 ~ 15,
      # gapswithin groups
      class_group == 6 ~ 2,
      # gap between groups bottom
      class_group == 4 ~ 10 #10
      
    )
  )


brier_score_decomposition_all <- ggplot(brier_score_decomposition_barplot_all, aes(x = class_group, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack", width = brier_score_decomposition_barplot_all$wdth) +
  facet_grid(model ~ ., switch = "y") +
  coord_flip() +
  scale_fill_manual(
    values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
    ##
    breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
    ##
    name = ""
  ) +
  labs(title = "MSC-DSC Mean Brier Score Decomposition") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 15),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "points"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
    plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = ggplot2::element_text(size = 14)
  )

brier_score_decomposition_all

ggsave("final_plots/brier_decomposition_all_models.png",
       plot = brier_score_decomposition_all, width = 1.2 * 2297, height = 1.4 * 2181, dpi = 300, units = "px")


## -----------------------------------------------------------------------------
## LOG: MSC-DSC-plots
## -----------------------------------------------------------------------------



log_decomposition_results <- map2(
  .x = models_scoring_rules,
  .y = names(models_scoring_rules),
  .f = function(df, model_name) {
    # caclulate MSC-DSC-UNC 
    #@param score 
    # A string specifying the score function.
    #'   One of: `"Brier_score"` (default), `"log_score"`, `"MR_score"`.
    result <- mcbdsc(df %>% select(onset_prob_pred), y = df$actual_conflict, score="log_score")
    
    # Extrahiere estimates + Modellname
    estimates(result) %>%
      mutate(model = model_name)
  }
)


## ---
## In-depth: 8 models
## ---
# filter for selected models
log_decomposition_results_selected <- log_decomposition_results[selected_models]

## important remark for BRIER score:
# UNC = variance of X~Ber(p=E(y)=mean(y)) -> p(1-p)
#p <- models_scoring_rules$boot_240$actual_conflict
#mean(p)*(1-mean(p))

# combine to one dataframe
log_score_decomposition_selected <- bind_rows(log_decomposition_results_selected) %>%
  select(model, mean_score, MCB, DSC, UNC) %>%
  rename(MeanScore = mean_score)

log_score_decomposition_barplot_selected <- log_score_decomposition_selected %>%
  mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)

# data frame with format for the barchart
log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
  arrange(MeanScore) %>%  # sort by MeanScore
  pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap, gap1, gap2),
               names_to = "component",
               values_to = "value") %>%
  mutate(class = case_when(
    component %in% c("MeanScore") ~ "SCORE",
    component %in% c("MCB", "UNC") ~ "MCB_UNC",
    component %in% c("DSC", "score_invisible") ~ "DSC_score", 
    component %in% c("gap") ~ "GAP",
    component %in% c("gap1") ~ "GAPmeanscoreDSC",
    component %in% c("gap2") ~ "GAPdscUNC",
  )) %>%
  select(model, class, component, value)

log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
  mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))


# dataset for model orderbased on MeanScore
log_levs_selected <- log_score_decomposition_barplot_selected %>%
  filter(component == "MeanScore") %>%
  arrange(desc(value)) %>%
  pull(model)

# data for the plot
log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
  mutate(
    class_group = case_when(
      class == "MCB_UNC"   ~ 1,
      class == "DSC_score"       ~ 2,
      class == "SCORE" ~ 3,
      class == "GAP" ~ 4,
      class == "GAPmeanscoreDSC" ~ 5,
      class == "GAPdscUNC" ~ 6
      
    ),  
    model = factor(model, levels = log_levs_selected, ordered = TRUE)
  )

log_score_decomposition_barplot_selected$model <- recode(
  log_score_decomposition_barplot_selected$model,
  !!!model_labels
)


# set order of subbars within group
log_score_decomposition_barplot_selected$component <- factor(
  log_score_decomposition_barplot_selected$component,
  levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
)

# set order of sub bars within group
log_score_decomposition_barplot_selected$class_group <- factor(
  log_score_decomposition_barplot_selected$class_group,
  levels = c(4,1,6,2,3,5)    
)

# set width of mean score bar
mean_score_bar_selected <- 4.5
# set width of mcb,dcs,unc bars
msc_dsc_unc_bar_selected <- 2

# width column
log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
  mutate(
    wdth = case_when(
      # MeanScore
      class_group == 3 ~ mean_score_bar_selected, #4
      # MCB DSC
      class_group %in% c(1, 2) ~ msc_dsc_unc_bar_selected, #1
      # gapswithin groups
      class_group == 6 ~ msc_dsc_unc_bar_selected,
      # gap between groups top
      class_group == 5 ~ 15,
      #gap between groups
      class_group == 4 ~ 8 #10
      
    )
  )


log_score_decomposition_selected <- ggplot(log_score_decomposition_barplot_selected, aes(x = class_group, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack", width = log_score_decomposition_barplot_selected$wdth) +
  facet_grid(model ~ ., switch = "y") +
  coord_flip() +
  scale_fill_manual(
    values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
    ##
    breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
    ##
    name = ""
  ) +
  labs(title = "MSC-DSC Mean log Score Decomposition") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 18),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "points"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
    plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = ggplot2::element_text(size = 14)
  )

log_score_decomposition_selected


# ggsave("final_plots/log_score_decomposition_selected.png",
#        plot = log_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px")




## ---
## All models
## ---


# combine to one dataframe
log_score_decomposition_all <- bind_rows(log_decomposition_results) %>%
  select(model, mean_score, MCB, DSC, UNC) %>%
  rename(MeanScore = mean_score)

log_score_decomposition_barplot_all <- log_score_decomposition_all %>%
  mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)

# data frame with format for the barchart
log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
  arrange(MeanScore) %>%  # sort by MeanScore
  pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap, gap1, gap2),
               names_to = "component",
               values_to = "value") %>%
  mutate(class = case_when(
    component %in% c("MeanScore") ~ "SCORE",
    component %in% c("MCB", "UNC") ~ "MCB_UNC",
    component %in% c("DSC", "score_invisible") ~ "DSC_score", 
    component %in% c("gap") ~ "GAP",
    component %in% c("gap1") ~ "GAPmeanscoreDSC",
    component %in% c("gap2") ~ "GAPdscUNC",
  )) %>%
  select(model, class, component, value)

log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
  mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))


# dataset for model orderbased on MeanScore
log_levs_all <- log_score_decomposition_barplot_all %>%
  filter(component == "MeanScore") %>%
  arrange(desc(value)) %>%
  pull(model)

# data for the plot
log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
  mutate(
    class_group = case_when(
      class == "MCB_UNC"   ~ 1,
      class == "DSC_score"       ~ 2,
      class == "SCORE" ~ 3,
      class == "GAP" ~ 4,
      class == "GAPmeanscoreDSC" ~ 5,
      class == "GAPdscUNC" ~ 6
      
    ),  
    model = factor(model, levels = log_levs_all, ordered = TRUE)
  )

log_score_decomposition_barplot_all$model <- recode(
  log_score_decomposition_barplot_all$model,
  !!!model_labels
)


# set order of subbars within group
log_score_decomposition_barplot_all$component <- factor(
  log_score_decomposition_barplot_all$component,
  levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
)

# set order of sub bars within group
log_score_decomposition_barplot_all$class_group <- factor(
  log_score_decomposition_barplot_all$class_group,
  levels = c(4,1,6,2,3,5)
)

# set width of mean score bar
mean_score_bar_all <- 7
# set width of mcb,dcs,unc bars
msc_dsc_unc_bar_all <- 2

# width column
log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
  mutate(
    wdth = case_when(
      # MeanScore
      class_group == 3 ~ mean_score_bar_all, #4
      # MCB DSC
      class_group %in% c(1, 2) ~ msc_dsc_unc_bar_all,#msc_dsc_unc_bar_all, #1
      # gap between groups top
      class_group == 5 ~ 15,
      # gapswithin groups
      class_group == 6 ~ 2,
      # gap between groups bottom
      class_group == 4 ~ 10 #10
      
    )
  )


log_score_decomposition_all <- ggplot(log_score_decomposition_barplot_all, aes(x = class_group, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "stack", width = log_score_decomposition_barplot_all$wdth) +
  facet_grid(model ~ ., switch = "y") +
  coord_flip() +
  scale_fill_manual(
    values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
    ##
    breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
    ##
    name = ""
  ) +
  labs(title = "MSC-DSC Mean log Score Decomposition") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 15),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "points"),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
    plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
    legend.title = ggplot2::element_text(size = 15),
    legend.text = ggplot2::element_text(size = 14),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = ggplot2::element_text(size = 14)
  )

log_score_decomposition_all

# ggsave("final_plots/log_score_decomposition_all_models.png",
#        plot = log_score_decomposition_all, width = 1.2 * 2297, height = 1.4 * 2181, dpi = 300, units = "px")
