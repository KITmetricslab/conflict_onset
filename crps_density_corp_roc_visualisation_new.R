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
# data_path <- "../Data/"
data_path <- ifelse(os == "Windows",
                    "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/",  # Windows
                    "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/")  # macOS/Linux

# round predictive samples to integers or not
round_samples <- TRUE

# actual data from 2018-2023 in the directory
files_actuals_from18 <- list.files(data_path, pattern = "cm_actuals_\\d{4}\\.parquet", full.names = TRUE)
# read all files and bind them into a single data frame
observations_18_23 <- do.call(rbind, lapply(files_actuals_from18, arrow::read_parquet)) %>% # observations from 2018 - 2023
  arrange(country_id, month_id)

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

appendix_names <- paste0("_cm_Y20", c(18, 19, 20, 21, 22, 23))
subfolder_names <- c(18, 19, 20, 21, 22, 23)
predictive_samples <- list()

for (m in 1:length(model_names)) {
  model_files <- paste0("../Data/", "all_available_cm_predictions/", model_names[m], "/cm/window=Y20", subfolder_names, "/", model_names[m], appendix_names, ".parquet")
  print(m)
  # model_files <- paste0(data_path, "all_available_cm_predictions/", model_names[m], "/cm/window=Y20", subfolder_names, "/", model_names[m], appendix_names, ".parquet")
  predictive_samples[[m]] <- lapply(model_files, arrow::read_parquet)
  predictive_samples[[m]] <- bind_rows(predictive_samples[[m]]) %>%
    arrange(country_id, month_id) %>%
    filter(country_id %in% unique(observations_18_23$country_id))
  predictive_samples[[m]]$member <- rep(1:1000, nrow(observations_18_23)) # relevant because some models count samples 0:999 and others 1:1000

  if("draw" %in% names(predictive_samples[[m]])) { # relevant for one model
    predictive_samples[[m]] <- predictive_samples[[m]] %>% select(-draw)
  }
}

names(predictive_samples) <- model_names

if(round_samples == TRUE) predictive_samples <- lapply(predictive_samples, function(df) round(df, 0))

predictive_samples_wide <- lapply(predictive_samples, function(df) df %>%
                                    pivot_wider( # reshape so that there's one sample per row (instead of stacked columns)
                                      id_cols = c(country_id, month_id),
                                      names_from = "member",
                                      values_from = "outcome",
                                      names_prefix = "sample"
                                    ))

# merge with outcome data
complete_data <- lapply(predictive_samples_wide, function(df) merge(observations_18_23, df) %>%
                          arrange(country_id, month_id))

# set country ids and observations etc
country_ids <- unique(observations_18_23$country_id)
actuals_ids <- unique(observations_18_23$month_id)

date_seq <- seq(
  from = as.Date("2018-01-01"),
  to   = as.Date("2023-12-01"),
  by   = "month"
)

# format as "MM-YYYY"
month_labels <- format(date_seq, "%m-%Y")

# Named character vector: names are month_ids
month_lookup_vec <- setNames(month_labels, actuals_ids)


### CRPS - BRIER - LOG Score Calculation -------------------------------------------------------------------------

## -----
# compute crps on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----

models_crps <- lapply(complete_data, function(df) apply(df, 1, function(df_row) {
  crps_sample(y = df_row[3],
              dat = df_row[4:1003])
}))

models_crps <- lapply(models_crps, function(df) data.frame(observations_18_23[,c("country_id", "month_id")], "crps" = df))

# save(models_crps, file = "output/models_crps.RData")
# load("output/models_crps.RData")



## ------------------------------------------------------------------------------------------------------------
## define lower threshold for binary event of a present conflict
## -----
# conflict = monthly_fatalities > fatality_thresh
lower_fatalitiy_thresh = 0
#lower_fatalitiy_thresh = 24
## ------------------------------------------------------------------------------------------------------------

# compute empirical probabilities for the binary onset event
models_predictive_probabilities <- lapply(predictive_samples, function(pred_sample) {
  pred_sample %>%
    mutate("predicted_conflict" = outcome > lower_fatalitiy_thresh) %>%
    group_by(country_id, month_id) %>%
    summarise(predictive_probability = mean(predicted_conflict),
              predictive_probability_log_nplustwo = (sum(predicted_conflict)+1)/(length(outcome)+2))
})


# merge models_predictive_probabilities and models_crps into new list "models_scoring_rules"
models_scoring_rules <- list()
for (m in 1:n_models) {
  models_scoring_rules[[m]] <- models_crps[[m]] %>%
    left_join(models_predictive_probabilities[[m]], by = c("country_id", "month_id")) %>%
    rename(onset_prob_pred = predictive_probability,
           onset_prob_pred_nplustwo = predictive_probability_log_nplustwo)
}
names(models_scoring_rules) <- model_names

# add the acutal observations to list
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    left_join(observations_18_23 %>%
                select(country_id, month_id, actual = outcome),
              by = c("country_id", "month_id"))
})


## -----
# compute Brier score on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----
# compute summands of brier score for each model, month and country
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    mutate(
      actual_conflict = actual > lower_fatalitiy_thresh,
      brier_onset = (actual_conflict - onset_prob_pred)^2
    )
})

# remove lists that are not longer needed
rm(models_predictive_probabilities, models_crps)


## -----
# compute log-score on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----
# compute log-scores
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    mutate(
      log_score_onset = - (actual_conflict * log(onset_prob_pred_nplustwo) +
                             (1 - actual_conflict) * log(1 - onset_prob_pred_nplustwo))
    )
})


## -----
# perform model checks (are there instances where crps < brier?)
## -----
## get combinations (model, country_id, month_id) where crps < brier_onset
crps_less_brier <- lapply(models_scoring_rules, function(df) df[which(round(df$crps,5) < round(df$brier_onset,5)), c("country_id", "month_id")])
sum_crps_less_brier <- unlist(lapply(crps_less_brier, nrow))

## check whether submitted forecasts per model consist of only integer values or not
only_integer <- unlist(lapply(predictive_samples, function(df) ifelse(sum(df$outcome%%1!=0)>0, FALSE, TRUE)))

## create combined overview
model_check <- data.frame(only_integer, sum_crps_less_brier); model_check


## check whether the forecasts for only-integer models and (country_id, month_id)
## with crps < brier_onset are only zeros or also other values
problem_model_ids <- which(model_check$only_integer & model_check$sum_crps_less_brier>0)
integer_crps_less_brier_data <- lapply(problem_model_ids,
                                       function(model) {
                                         df1 <- merge(crps_less_brier[[model]], models_scoring_rules[[model]])
                                         merge(df1, complete_data[[model]])
                                       })

names(integer_crps_less_brier_data) <- rownames(model_check)[problem_model_ids]


## -----
## REMOVE excluded models
## -----
# exclude some models for all plots and analysis
# excluded_models <- c("Neg_Bin_GAM", "P_GAM", "TW_GAM")
#
# # delete models
# models_scoring_rules <- models_scoring_rules[
#   ! names(models_scoring_rules) %in% excluded_models
# ]
#
# model_names <- names(models_scoring_rules)
# n_models <- length(model_names)





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
  arrange(country_id, month_id) %>%
  mutate(conflict = outcome>lower_fatalitiy_thresh)


actual_conflict <- observations_17_23 %>%
  filter(country_id %in% country_ids) %>%
  filter(month_id %in% actuals_ids) %>%
  mutate(prev_oct_id = (actuals_ids[1]-3) + 12 * ((month_id - 457) %/% 12))


actual_conflict <- actual_conflict %>%
  left_join(
    observations_17_23 %>%
      select(country_id,
             month_id,
             conflict_prev_oct = conflict),
    by = c("country_id",
           "prev_oct_id" = "month_id")
  )

situation_month <- ifelse(!actual_conflict$conflict & !actual_conflict$conflict_prev_oct, "peace", # no conflict, no conflict
                          ifelse(actual_conflict$conflict & !actual_conflict$conflict_prev_oct, "onset", # no conflict, conflict
                                 ifelse(actual_conflict$conflict & actual_conflict$conflict_prev_oct, "conflict", # conflict, conflict
                                        "deescalation"))) # conflict, no conflict


conflict_situations <- data.frame(actual_conflict[,c(2,3)], "situation_month" = as.vector(situation_month))
# write.csv(conflict_situations, "conflict_situations.csv")

## -----
## Plot Data: CRPS values by conflict situation
## -----
models_crps_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "crps")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_crps_conflict) <- c("country_id", "month_id", model_names)
models_crps_conflict <- list(models_crps_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

models_crps_conflict_month <- models_crps_conflict %>% select(!c("country_id", "month_id")) %>% group_by(situation_month) %>% summarise_all(sum)

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
models_brier_conflict_month <- models_brier_conflict %>% select(!c("country_id", "month_id")) %>% group_by(situation_month) %>% summarise_all(sum)

models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] <- models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] / nrow(models_brier_conflict) # compute contributions to average brier

brier_month <- data.frame("Brier" = unlist(c(models_brier_conflict_month[,2:ncol(models_brier_conflict_month)])),
                          "Situation" = rep(models_brier_conflict_month$situation_month, ncol(models_brier_conflict_month)-1),
                          "Model" = rep(names(models_brier_conflict_month)[2:ncol(models_brier_conflict_month)], each = 4))



## ---
## Plot-Data: log-score values by conflict situation for the onset problem (y \in {0,1})
## ---
models_logscore_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "log_score_onset")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_logscore_conflict) <- c("country_id", "month_id", model_names)
models_logscore_conflict <- list(models_logscore_conflict, conflict_situations) %>% reduce(left_join, c("country_id", "month_id"))

models_logscore_conflict_month <- models_logscore_conflict %>% select(!c("country_id", "month_id")) %>% group_by(situation_month) %>% summarise_all(sum)

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
  situation_month = character()
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


## -----
## save plots in folders
## ----
store_plot <- TRUE

## -----
## labels, colors, textsize etc.
## ----

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

model_labels_df <- data.frame("name_original" = names(model_labels), "name_paper" = model_labels)

## -----
## Selected models
## -----

## --
## change this if needed. everything else is dynamic!
# selected_models <- c("zero","last",
#                      "conflictology","conflictforecast_v2",
#                      "submission_muchlinski_thornhill", "submission_final_omm",
#                      "unito_transformer", "P_GLMM")
selected_models <- c("bodentien_rueter_negbin", "conflictforecast_v2", # replaced P_GLMM with bodentien_rueter_negbin
                     "submission_muchlinski_thornhill", "ShapeFinder",
                     "submission_final_omm", "unito_transformer",
                     "conflictology", "zero")
##
## --



selected_colors <- c(
  "#009682",  # green
  "#4664aa",  # blue
  "#23a1e0",  # maygreen
  "black",    # grey
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

# remaining_models <- setdiff(model_names, selected_models)
#
# remaining_colors <- c(
#   "#009682",  # green
#   "#4664aa",  # blue
#   "#23a1e0",  # maygreen
#   "black", # grey
#   "#df9b1b",  # orange
#   "#8cb63c",  # yellow
#   "#a22223",  # red
#   "#a3107c",   # purple
#   "#51103C"    # additional color
# )
#
#
# remaining_model_colors <- setNames(remaining_colors, remaining_models)
#
# remaining_model_labels <- model_labels[remaining_models]
#
# model_colors <- c(
#   selected_model_colors,
#   remaining_model_colors
# )




theme_fontsize <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
  axis.title = ggplot2::element_text(size = 13),
  axis.text = ggplot2::element_text(size = 12),
  legend.text = element_text(size = 12),
)

theme_fontsize_large <- ggplot2::theme(
  plot.title = ggplot2::element_text(size = 22, hjust = 0.5),
  axis.title = ggplot2::element_text(size = 16),
  axis.text = ggplot2::element_text(size = 16),
  legend.text = element_text(size = 16),
)

# margin and subtitle for grid plots
sg <- textGrob('', gp = gpar(fontsize = 4))
margin <- unit(0.5, "line")


## -----------------------------------------------------------------------------
## Figure 1a Number of total fatalities worldwide per month, stacked bar plot
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

# to determine maximum y range
fatalities_worldwide_y_range <- fatalities_worldwide_plot_data %>%
  group_by(month_id) %>%
  summarise(
    total_fatalities = sum(n_fatalities)
  )
max(fatalities_worldwide_y_range$total_fatalities)



fatalities_worldwide_plot <- ggplot(fatalities_worldwide_plot_data, aes(fill=country_category, y=n_fatalities, x=month_id)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(
    values = c("Ukraine" = "#E69F00", "Yemen" = "#0072B2",
               "Afghanistan"="#D55E00","Syria" = "#51103C",
               "Ethiopia" = "#009E73", "others" = "grey30"),
    breaks = c("Ethiopia", "Ukraine", "Afghanistan", "Yemen", "Syria", "others"),
    ##
    labels = c("Ethiopia", "Ukraine", "Afghanistan", "Yemen", "Syria", "Others")
  ) +
  labs(title = "Aggregated Conflict Fatalities by Month and Country",
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
  ) + theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  theme_fontsize

# --- find December→January transitions ---
month_boundaries <- fatalities_worldwide_plot_data %>%
  mutate(month_num = (month_id - 1) %% 12 + 1) %>%
  filter(month_num == 12) %>%
  pull(month_id) %>%
  unique() %>%
  sort()

# --- get October bar tops ---
october_tops <- fatalities_worldwide_plot_data %>%
  mutate(month_num = (month_id - 1) %% 12 + 1) %>%
  filter(month_num == 10) %>%
  group_by(month_id) %>%
  summarise(y_pos = sum(n_fatalities), .groups = "drop")

# --- add to plot ---
fatalities_worldwide_plot <- fatalities_worldwide_plot +
  geom_vline(xintercept = month_boundaries + 0.5,   # dashed lines after each December
             linetype = "dashed", color = "black") +
  geom_point(data = october_tops,                   # triangles above October bars
             aes(x = month_id, y = y_pos + 3000),   # offset above bar
             shape = 6, size = 2.5, fill = "black")

fatalities_worldwide_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/fatalities_worldwide_1.png",
           plot = fatalities_worldwide_plot, width = 1.2 * 3000, height = 1.2 * 1322, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/fatalities_worldwide_25.png",
           plot = fatalities_worldwide_plot, width = 1.2 * 3000, height = 1.2 * 1322, dpi = 300, units = "px",
           bg="white")
  }
}



## -----------------------------------------------------------------------------
## Figure 1b: Number of Onset events per month for the test window
## -----------------------------------------------------------------------------
onsets_per_month_plot_data <- conflict_situations %>%
  group_by(month_id) %>%
  summarise(
    n_onset = sum(situation_month == "onset", na.rm = TRUE),
    .groups  = "drop"
  )

# total number of onsets
sum(onsets_per_month_plot_data$n_onset)

# unique countries that experienced onset
nrow(unique(conflict_situations %>% filter(situation_month == "onset") %>% select(country_id)))

onset_country_year <- merge(conflict_situations %>% filter(situation_month == "onset"),
                            data.frame("month_id" = names(month_lookup_vec), "month_year" = month_lookup_vec, "year" = sub('.*(\\d{4}).*', '\\1', month_lookup_vec)))

# write.csv(onset_country_year, "onset_country_year.csv")
unique(onset_country_year %>% select(country_id, year)) %>% arrange(country_id)

onsets_per_month_plot <- ggplot(onsets_per_month_plot_data, aes(y = n_onset, x = month_id)) +
  geom_bar(position="stack", stat="identity",fill  = "darkgrey", alpha = 0.7, width = 1, color = "black") +
  labs(title = "Onset Events by Month for all Countries",
       x = "Month") +
  scale_x_continuous(
    breaks = seq(min(fatalities_months), max(fatalities_months), by = 6),
    labels = month_lookup_vec[as.character(seq(min(fatalities_months), max(fatalities_months), by = 6))]
  ) +
  scale_y_continuous(
    "Onset Events",
    breaks = seq(0,14, by = 2)
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  theme_fontsize +
  geom_vline(xintercept = month_boundaries + 0.5,   # dashed lines after each December
             linetype = "dashed", color = "black")

onsets_per_month_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/onsets_per_month_1.png",
           plot = onsets_per_month_plot, width = 1.2 * 3000, height = 1.2 * 1322, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/onsets_per_month_25.png",
           plot = onsets_per_month_plot, width = 1.2 * 3000, height = 1.2 * 1322, dpi = 300, units = "px",
           bg="white")
  }
}



## -----------------------------------------------------------------------------
## Figure 1: 1a & 1b combined
## -----------------------------------------------------------------------------

library(patchwork)
fatalities_onsets_plot <- fatalities_worldwide_plot / onsets_per_month_plot +
  plot_layout(heights = c(4, 3))  # adjust ratio of top/bottom

fatalities_onsets_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/fatalities_onset_greq1.png",
           plot = fatalities_onsets_plot, width = 1.2 * 3000, height = 1.2 * 2000, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/fatalities_onset_greq25.png",
           plot = fatalities_onsets_plot, width = 1.2 * 3000, height = 1.2 * 2000, dpi = 300, units = "px",
           bg="white")
  }
}



## -----------------------------------------------------------------------------
## Figure 2a: CRPS per Country: highlighting 5-10 most important ones (and other)
# for the years 2018 - 2023
## -----------------------------------------------------------------------------






## -----------------------------------------------------------------------------
## Figure 2b: CRPS decomposition
## -----------------------------------------------------------------------------
## ---
## Selected models
## ---
crps_month_selected_models <- crps_month %>% filter(Model %in% selected_models) %>%
  mutate("Model_orig" = Model)

# Rename models
crps_month_selected_models$Model <- recode(
  crps_month_selected_models$Model,
  !!!model_labels
)

# plot
library(forcats)

crps_conflict_situation_plot <- crps_month_selected_models %>%
  group_by(Model) %>%
  summarise(total_crps = sum(CRPS), .groups = "drop") %>%
  right_join(crps_month_selected_models, by = "Model") %>%
  mutate(Model = fct_reorder(Model, total_crps, .desc = TRUE)) %>%
  ggplot(aes(fill = Situation, y = Model, x = CRPS)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean CRPS per Conflict Situation",
       x = "Mean CRPS") +
  scale_fill_manual("Situation",
                    values = c("conflict" = "#a22223",
                               "deescalation" = "#d09191",
                               "onset" = "#4664aa",
                               "peace" = "#a2b2d4"),
                    labels = c("conflict" = "Continued conflict",
                               "deescalation" = "End / interruption of conflict",
                               "onset" = "Onset",
                               "peace" = "Continued peace")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  theme_fontsize

crps_conflict_situation_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/crps_conflict_situation_1.png",
           plot = crps_conflict_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/crps_conflict_situation_25.png",
           plot = crps_conflict_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  }
}



## -----------------------------------------------------------------------------
## Figure 2c: CRPS decomposition for previous peace with Brier score
## -----------------------------------------------------------------------------
## ---
## Selected models
## ---
crps_selected_models_prev_peace <- crps_month %>% filter(Model %in% selected_models) %>%
  filter(Situation %in% c("peace", "onset")) %>%
  mutate(Model_name = Model)

model_order <- unlist(crps_selected_models_prev_peace %>%
  group_by(Model) %>%
  summarise(total_crps = sum(CRPS), .groups = "drop") %>%
  arrange(total_crps) %>%
  select(Model))

model_order <- rev(model_order)

brier_selected_models_prev_peace <- brier_month %>% filter(Model %in% selected_models) %>%
  filter(Situation %in% c("peace", "onset"))

# Compute CRPS remainder
crps_brier_selected_models_prev_peace <- merge(crps_selected_models_prev_peace, brier_selected_models_prev_peace) %>%
  mutate("CRPS_Remainder" = CRPS-Brier)

# Create CRPS remainder part
crps_brier_plot_data_Rem <- crps_brier_selected_models_prev_peace %>%
  mutate(Situation = paste0(Situation, "-Remainder")) %>%
  rename(Score = CRPS_Remainder) %>%
  select(-c("CRPS", "Brier"))
crps_brier_plot_data_BS <- crps_brier_selected_models_prev_peace %>%
  mutate(Situation = paste0(Situation, "-Brier")) %>%
  rename(Score = Brier) %>%
  select(-c("CRPS", "CRPS_Remainder"))

crps_brier_plot_data <- rbind(crps_brier_plot_data_Rem, crps_brier_plot_data_BS) %>%
  arrange(factor(Model, levels = model_order))

crps_brier_plot_data <- crps_brier_plot_data %>%
  mutate(Model = factor(Model, levels = model_order))

# Rename models
crps_brier_plot_data$Model <- recode(
  crps_brier_plot_data$Model,
  !!!model_labels
)

# Plot
crps_brier_prev_peace_plot <- crps_brier_plot_data  %>%
  ggplot(aes(fill = Situation, y = Model, x = Score)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean CRPS for Previous Peace with Brier Score",
       x = "Mean CRPS") +
  scale_fill_manual("Situation",
                    values = c("onset-Remainder" = "#4664aa",
                               "peace-Remainder" = "#a2b2d4",
                               "onset-Brier" = "#555555",
                               "peace-Brier" = "lightgrey"),
                    labels = c("onset-Remainder" = "Onset (CRPS remainder)",
                               "peace-Remainder" = "Continued peace (CRPS remainder)",
                               "onset-Brier" = "Onset (Brier score)",
                               "peace-Brier" = "Continued peace (Brier score)")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  theme_fontsize

crps_brier_prev_peace_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/crps_brier_prev_peace_1.png",
           plot = crps_brier_prev_peace_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/crps_brier_prev_peace_25.png",
           plot = crps_brier_prev_peace_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  }
}




## -----------------------------------------------------------------------------
## Figure 2d: Brier score decomposition for previous peace
## -----------------------------------------------------------------------------
## ---
## Selected models
## ---
brier_selected_models_prev_peace <- brier_month %>% filter(Model %in% selected_models) %>%
  filter(Situation %in% c("peace", "onset"))

model_order_brier <- unlist(brier_selected_models_prev_peace %>%
                              group_by(Model) %>%
                              summarise(total_brier = sum(Brier), .groups = "drop") %>%
                              arrange(total_brier) %>%
                              select(Model))

model_order_brier <- rev(model_order_brier)

brier_selected_models_prev_peace <- brier_selected_models_prev_peace %>%
  arrange(factor(Model, levels = model_order_brier))

brier_selected_models_prev_peace <- brier_selected_models_prev_peace %>%
  mutate(Model = factor(Model, levels = model_order_brier))

# Rename models
brier_selected_models_prev_peace$Model <- recode(
  brier_selected_models_prev_peace$Model,
  !!!model_labels
)


# Plot
brier_prev_peace_plot <- brier_selected_models_prev_peace  %>%
  ggplot(aes(fill = Situation, y = Model, x = Brier)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean Brier Score for Previous Peace",
       x = "Mean CRPS") +
  scale_fill_manual("Situation",
                    values = c("onset" = "#555555",
                               "peace" = "lightgrey"),
                    labels = c("onset" = "Onset (Brier score)",
                               "peace" = "Continued peace (Brier score)")) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  theme_fontsize

brier_prev_peace_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/brier_prev_peace_1.png",
           plot = brier_prev_peace_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/brier_prev_peace_25.png",
           plot = brier_prev_peace_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  }
}


crps_brier_prev_peace_comb_plot <- crps_brier_prev_peace_plot / brier_prev_peace_plot +
  plot_layout(heights = c(1, 1), guides = "keep") &
  theme(
    # same x for both legends -> horizontally aligned
    legend.position = c(1, 0.55),        # tweak 1.02 -> 1.15 if you want it further right
    legend.justification = c(0, 0.5),     # anchor left-middle of the legend box
    legend.spacing.y = unit(0.2, "cm"),   # spacing within legend if stacked items
    plot.margin = margin(5, 115, 5, 5)     # increase right margin so legends have room
  ) &
  coord_cartesian(clip = "off")             # prevent clipping of legend


if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/crps_brier_comb_prev_peace_1.png",
           plot = crps_brier_prev_peace_comb_plot, width = 1.0 * 4222, height = 1.0 * 2600, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/crps_brier_comb_prev_peace_25.png",
           plot = crps_brier_prev_peace_comb_plot, width = 1.0 * 4222, height = 1.0 * 2600, dpi = 300, units = "px",
           bg="white")
  }
}


crps_brier_comb_plot <- crps_conflict_situation_plot / crps_brier_prev_peace_plot / brier_prev_peace_plot +
  plot_layout(heights = c(1, 1, 1), guides = "keep") &
  theme(
    # same x for both legends -> horizontally aligned
    legend.position = c(1, 0.55),        # tweak 1.02 -> 1.15 if you want it further right
    legend.justification = c(0, 0.5),     # anchor left-middle of the legend box
    legend.spacing.y = unit(0.2, "cm"),   # spacing within legend if stacked items
    plot.margin = margin(5, 115, 5, 5)     # increase right margin so legends have room
  ) &
  coord_cartesian(clip = "off")             # prevent clipping of legend

crps_brier_comb_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/crps_brier_comb_1.png",
           plot = crps_brier_comb_plot, width = 1.0 * 4222, height = 1.0 * 3900, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/crps_brier_comb_25.png",
           plot = crps_brier_comb_plot, width = 1.0 * 4222, height = 1.0 * 3900, dpi = 300, units = "px",
           bg="white")
  }
}

# library(magick)
#
# if(store_plot == TRUE){
#   if(lower_fatalitiy_thresh == 0){
#     img <- image_read("plots_fatalities_greq1/crps_brier_comb_1.png")
#   } else {
#     img <- image_read("plots_fatalities_greq25/crps_brier_comb_25.png")
#   }
#   info <- image_info(img)
#   width <- info$width
#   height <- info$height
#   part_height <- floor(height / 3)
#
#   for (i in 0:2) {
#     top <- i * part_height
#     bottom <- if (i < 2) part_height else height - top  # include remainder in last part
#     part <- image_crop(img, geometry = geometry_area(width, bottom, 0, top))
#     if(lower_fatalitiy_thresh == 0){
#       image_write(part, path = paste0("plots_fatalities_greq1/crps_brier_comb_part_", i + 1, "_1.png"))
#     } else {
#       image_write(part, path = paste0("plots_fatalities_greq25/crps_brier_comb_part_", i + 1, "_25.png"))
#     }
#   }
# }

## -----------------------------------------------------------------------------
## Figure 3: ROC curves
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

mm_eval <- evalmod(mm_selected_models)
auc_selected_models <- cbind(auc(mm_eval) %>% filter(curvetypes == "ROC"), selected_model_labels)


labels_AUC <- paste0(auc_selected_models$selected_model_labels, ", AUC = ", round(auc_selected_models$aucs, 2))

roc_curve_selected_models <- autoplot(mm_eval, curvetype = "ROC") +
  ggplot2::geom_line(size = 1, alpha = 0.8) +
  ggplot2::scale_color_manual(
    values = selected_model_colors,
    breaks = names(selected_model_labels),
    labels = labels_AUC
  ) +
  ggplot2::labs(
    title = "ROC Curve",
    x = "FAR",
    y = "HR",
    color = ""
  ) +
  # ggplot2::theme(legend.position = "none") +
  theme_fontsize

roc_curve_selected_models

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/roc_curve_selected_models_1.png",
           plot = roc_curve_selected_models, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/roc_curve_selected_models_25.png",
           plot = roc_curve_selected_models, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
           bg="white")
  }
}


## -----------------------------------------------------------------------------
## Figure 4: Reliability diagrams
## -----------------------------------------------------------------------------
# -------------------------------------------------------------------
# Define model labels and colors (must exist for your function)
# -------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(reliabilitydiag)

# -------------------------------------------------------------------
# Title grob for the combined plot
# -------------------------------------------------------------------
tg <- textGrob("Reliability Diagrams", gp = gpar(fontsize = 18, hjust = 0.5))
# theme_fontsize <- ggplot2::theme(
#   plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
#   axis.title = ggplot2::element_text(size = 13),
#   axis.text = ggplot2::element_text(size = 12),
#   legend.text = element_text(size = 12),
# )


# -------------------------------------------------------------------
# Define model labels and colors
# -------------------------------------------------------------------
model_labels_selected_models <- model_labels_df %>%
  filter(name_original %in% selected_models) %>%
  arrange(factor(name_original, levels = selected_models))

model_labels <- setNames(model_labels_selected_models$name_paper, model_labels_selected_models$name_original)

# Pick a palette with at least as many colors as models
# selected_colors <- RColorBrewer::brewer.pal(n = length(selected_models), "Set2")
names(selected_colors) <- selected_models

# -------------------------------------------------------------------
# Function: Reliability Diagram
# -------------------------------------------------------------------
create_reliability_diag <- function(data, forecast_model) {

  # subset for selected model
  reliability_data_selected <- data %>%
    dplyr::filter(model == forecast_model)

  # compute reliability diagram
  r_selected <- reliabilitydiag(
    x = reliability_data_selected$onset_prob_pred,
    y = reliability_data_selected$actual
  )

  # base plot
  reliability_plot <- autoplot(r_selected)

  # strip out unwanted geom_segment layers
  is_seg <- sapply(reliability_plot$layers, function(layer) {
    inherits(layer$geom, "GeomSegment")
  })
  reliability_plot$layers <- reliability_plot$layers[!is_seg]

  # CEP estimates
  data_estim <- estimates(
    reliability(
      x = reliability_data_selected$onset_prob_pred,
      y = reliability_data_selected$actual
    )
  ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(CEP)   # ensure order

  # build horizontal segment df
  df_segments <- data.frame()
  if (nrow(data_estim) > 1) {
    for (i in 1:(nrow(data_estim) - 1)) {
      if (data_estim$CEP[i] == data_estim$CEP[i + 1]) {
        df_segments <- rbind(
          df_segments,
          data.frame(
            x = min(data_estim$x[i], data_estim$x[i + 1]),
            x_end = max(data_estim$x[i], data_estim$x[i + 1]),
            CEP = data_estim$CEP[i],
            CEP_end = data_estim$CEP[i + 1]
          )
        )
      }
    }
  }

  # color for this model
  path_color <- selected_colors[forecast_model]

  # final plot
  p <- reliability_plot +
    ggplot2::labs(
      title = model_labels[forecast_model],
      x = "Forecast value",
      y = "CEP"
    ) +
    # theme_fontsize +
    # ggplot2::theme(
    #   legend.position = "none") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::geom_path(
      mapping = ggplot2::aes(x = .data$x, y = .data$CEP),
      data = data_estim,
      linewidth = 1,
      colour = path_color
    )

  # add flat horizontal segments if present
  if (nrow(df_segments) > 0) {
    p <- p + ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$x, y = .data$CEP,
        xend = .data$x_end, yend = .data$CEP_end
      ),
      data = df_segments,
      linewidth = 1.4,
      colour = path_color
    )
  } else {
    # single-point case
    p <- p + ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$x, y = .data$CEP),
      data = data_estim,
      colour = path_color,
      shape = 19,
      size = 2
    )
  }

  return(p)
}

# -------------------------------------------------------------------
# Generate plots for all selected models
# -------------------------------------------------------------------
corp_plots_list_selected <- list()
for (model_name in selected_models) {
  corp_plots_list_selected[[model_name]] <- create_reliability_diag(
    prev_peace_prob_month_long_binary_actual,
    model_name
  )
}

# -------------------------------------------------------------------
# Arrange them in a grid with a title
# -------------------------------------------------------------------
grid.arrange(
  tg,
  do.call(arrangeGrob, c(corp_plots_list_selected, ncol = 4)),
  ncol = 1,
  heights = c(0.1, 1)
)

if (store_plot == TRUE) {

  # build the arranged plot as a grob
  reliability_grid <- gridExtra::arrangeGrob(
    tg,
    do.call(arrangeGrob, c(corp_plots_list_selected, ncol = 4)),
    ncol = 1,
    heights = c(0.1, 1)
  )

  if (lower_fatalitiy_thresh == 0) {
    ggsave(
      "plots_fatalities_greq1/reliability_diagram_selected_models_1.png",
      plot = reliability_grid,
      width = 1.0 * 3000, height = 1.0 * 1800, dpi = 300, units = "px",
      bg = "white"
    )
  } else {
    ggsave(
      "plots_fatalities_greq25/reliability_diagram_selected_models_25.png",
      plot = reliability_grid,
      width = 1.0 * 3000, height = 1.0 * 1800, dpi = 300, units = "px",
      bg = "white"
    )
  }
}


# -------------------------------------------------------------------
# Align ROC curve an select RV O MM (submission_final_omm)
# -------------------------------------------------------------------

roc_reliability_plot <-
  (roc_curve_selected_models + theme(legend.position = "none")) /
  (corp_plots_list_selected$submission_final_omm +
     labs(title = "Reliability Diagram") +
     theme_fontsize) +
  plot_layout(heights = c(1, 1))


if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plots_fatalities_greq1/roc_reliability_greq1.png",
           plot = roc_reliability_plot, width = 1.2 * 1000, height = 1.2 * 2000, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plots_fatalities_greq25/roc_reliability_greq25.png",
           plot = roc_reliability_plot, width = 1.2 * 1000, height = 1.2 * 2000, dpi = 300, units = "px",
           bg="white")
  }
}

