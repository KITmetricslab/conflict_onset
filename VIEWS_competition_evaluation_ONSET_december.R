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
library(pracma)
library(forcats)

#devtools::install_github("aijordan/reliabilitydiag")

## ---------
# load observational data (actuals)
## ---------

# path to directory on "share-alle"
data_path <- "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/" # Windows version
#data_path <- "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/" # MacOS version


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
  #model_files <- paste0("../Data/", "all_available_cm_predictions/", model_names[m], "/cm/window=Y20", subfolder_names, "/", model_names[m], appendix_names, ".parquet")
  
  print(m)
  # model_files <- paste0(data_path, "all_available_cm_predictions/", model_names[m], "/cm/window=Y20", subfolder_names, "/", model_names[m], appendix_names, ".parquet")
  model_files <- paste0(data_path, "all_available_cm_predictions/", model_names[m], "/cm/window=Y20", subfolder_names, "/", model_names[m], appendix_names, ".parquet")
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


## ------------------------------------------------------------------------------------------------------------
## define lower threshold for binary event of a present conflict
## -----
# conflict = monthly_fatalities > fatality_thresh
lower_fatalitiy_thresh = 24 #0
## ------------------------------------------------------------------------------------------------------------

# compute empirical probabilities for the binary onset event
models_predictive_probabilities <- lapply(predictive_samples, function(pred_sample) {
  pred_sample %>%
    mutate("predicted_conflict" = outcome > lower_fatalitiy_thresh) %>%
    group_by(country_id, month_id) %>%
    summarise(predictive_probability = mean(predicted_conflict),
              predictive_probability_log_nplustwo = (sum(predicted_conflict)+1)/(length(outcome)+2),
              .groups = "drop")
})

# merge models_predictive_probabilities and models_crps into new list "models_scoring_rules"
models_scoring_rules <- models_predictive_probabilities
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
      brier_onset = (actual_conflict - predictive_probability)^2,
      brier_onset_log_target = (1/(lower_fatalitiy_thresh + 1)) * (actual_conflict - predictive_probability)^2
    )
})

# remove lists that are not longer needed
rm(models_predictive_probabilities)


## -----
# compute log-score on benchmark models and submitted forecasts for each country-month pair for 2018 to 2023
## -----
# compute log-scores
models_scoring_rules <- lapply(models_scoring_rules, function(df) {
  df %>%
    mutate(
      log_score_onset = - (actual_conflict * log(predictive_probability_log_nplustwo) +
                             (1 - actual_conflict) * log(1 - predictive_probability_log_nplustwo))
    )
})






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


## ---
## Data for Onset Prediction
## ---

# dataframe to store onset probabilities in long format
prob_models_all_situation_long <- data.frame(
  model = character(),
  month_id = integer(),
  country_id = integer(),
  predictive_probability = numeric(),
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

# onset actuals binary target
prev_peace_prob_month_long_binary_actual <- prev_peace_prob_month_long %>%
  mutate(actual = ifelse(actual > lower_fatalitiy_thresh, 1, 0))




## ---------
## Filtering for only December
## ---------
prev_peace_prob_month_long_binary_actual <- prev_peace_prob_month_long_binary_actual %>%
  filter(month_id %% 12 == 0)








################################################################################
## PLOTS
################################################################################


## ---
## save plots in folders
## ---
store_plot <- TRUE

## ---
## labels, colors, textsize etc.
## --

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

## ---
## Selected models
## ---

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
## Figure 2d: Brier Score decomposition
## -----------------------------------------------------------------------------
## ---
## Prepare plot data: Brier values by conflict situation
## ---
models_brier_conflict <- lapply(models_scoring_rules, function(m) m %>% select("country_id", "month_id", "brier_onset")) %>%
  reduce(left_join, c("country_id", "month_id"))


names(models_brier_conflict) <- c("country_id", "month_id", model_names)
models_brier_conflict <- list(models_brier_conflict, conflict_situations) %>% 
  reduce(left_join, c("country_id", "month_id"))  %>%
  filter(month_id %% 12 == 0)

models_brier_conflict_month <- models_brier_conflict %>% 
  select(!c("country_id", "month_id")) %>% 
  group_by(situation_month) %>% 
  summarise_all(sum)

models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] <- models_brier_conflict_month[,2:ncol(models_brier_conflict_month)] / nrow(models_brier_conflict) # compute contributions to average brier


# create ggplot data frames
brier_month <- data.frame("Brier" = unlist(c(models_brier_conflict_month[,2:ncol(models_brier_conflict_month)])),
                          "Situation" = rep(models_brier_conflict_month$situation_month, ncol(models_brier_conflict_month)-1),
                          "Model" = rep(names(models_brier_conflict_month)[2:ncol(models_brier_conflict_month)], each = 4))



## ---
## Selected models
## ---
brier_month_selected_models <- brier_month %>% filter(Model %in% selected_models) %>%
  mutate("Model_orig" = Model)

# Rename models
brier_month_selected_models$Model <- recode(
  brier_month_selected_models$Model,
  !!!model_labels
)


brier_conflict_situation_plot <- brier_month_selected_models %>%
  group_by(Model) %>%
  summarise(total_brier = sum(Brier), .groups = "drop") %>%
  right_join(brier_month_selected_models, by = "Model") %>%
  mutate(Model = fct_reorder(Model, total_brier, .desc = TRUE)) %>%
  ggplot(aes(fill = Situation, y = Model, x = Brier)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean Brier Score per Conflict Situation",
       x = "Mean Brier Score") +
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

brier_conflict_situation_plot

if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plotsOnsetDecember_fatalities_greq1/DecOnset_brier_conflict_situation_1.png",
           plot = brier_conflict_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plotsOnsetDecember_fatalities_greq25/DecOnset_brier_conflict_situation_25.png",
           plot = brier_conflict_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
           bg="white")
  }
}


## -----------------------------------------------------------------------------
## Figure 3: ROC curves for previous peace
## -----------------------------------------------------------------------------
##
## https://cran.r-project.org/web/packages/precrec/vignettes/introduction.html
##
roc_data_prev_peace <- prev_peace_prob_month_long_binary_actual


## ---
## In-depth: 8 models
## ---
# filter relevant models
roc_data_prev_peace_selected <- roc_data_prev_peace %>%
  filter(model %in% selected_models)

# create list of scores
score_list_selected <- lapply(selected_models, function(m) {
  roc_data_prev_peace_selected %>%
    filter(model == m) %>%
    pull(predictive_probability)
})

# create list of labels
label_list_selected <- lapply(selected_models, function(m) {
  roc_data_prev_peace_selected %>%
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
    ggsave("plotsOnsetDecember_fatalities_greq1/DecOnset_roc_curve_selected_models_prev_peace_greq1.png",
           plot = roc_curve_selected_models, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plotsOnsetDecember_fatalities_greq25/DecOnset_roc_curve_selected_models_prev_peace_greq25.png",
           plot = roc_curve_selected_models, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
           bg="white")
  }
}


## -----------------------------------------------------------------------------
## Figure 4: Reliability diagrams for previous peace
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
tg <- textGrob("CORP Reliability Diagrams", gp = gpar(fontsize = 18, hjust = 0.5))
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
reliabilitydiag.custom <- function(fcst, obs, pathclr = "red", confnveau = 0.9, bndtype="diagonal", unc_mthd = "resampling", annt_score_decom = "large") {
  
  # compute reliability diagram
  r_selected <- reliabilitydiag(
    x = fcst,
    y = obs,
    region.level = confnveau,
    region.method = unc_mthd,
    region.position = bndtype
  )
  
  # base plot
  reliability_plot <- autoplot(r_selected)
  
  # # strip out unwanted geom_segment layers
  # is_seg <- sapply(reliability_plot$layers, function(layer) {
  #   inherits(layer$geom, "GeomSegment")
  # })
  # reliability_plot$layers <- reliability_plot$layers[!is_seg]
  
  
  
  
  is_seg <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  reliability_plot$layers <- reliability_plot$layers[!is_seg]
  
  
  
  # CEP estimates
  data_estim <- estimates(
    reliability(
      x = fcst,
      y = obs
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
  
  # final plot
  p <- reliability_plot +
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
      linewidth = 0.9,
      colour = pathclr
    ) + ggplot2::geom_segment(
      aes(
        x = 0,
        y = 0,
        xend = 1,
        yend = 1
      ),
      inherit.aes = FALSE,
      colour = "grey50",
      linewidth = 0.4
    )
  
  # add flat horizontal segments if present
  if (nrow(df_segments) > 0) {
    p <- p + ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = .data$x, y = .data$CEP,
        xend = .data$x_end, yend = .data$CEP_end
      ),
      data = df_segments,
      linewidth = 1.3,
      colour = pathclr
    )
  } else {
    # single-point case
    p <- p + ggplot2::geom_point(
      mapping = ggplot2::aes(x = .data$x, y = .data$CEP),
      data = data_estim,
      colour = pathclr,
      shape = 19,
      size = 1.7
    ) +
      ggplot2::coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
  }
  
  if (annt_score_decom == "small") {
    p <- p +
      annotate(
        "text",
        x = .125,
        y = 1.01,
        label = sprintf("MCB = .%03d",
                        round(summary(r_selected)$miscalibration * 1000)),
        color = "red",
        size = 2
      ) +
      annotate(
        "text",
        x = .125,
        y = .95,
        label = sprintf("DSC = .%03d",
                        round(summary(r_selected)$discrimination * 1000)),
        size = 2
      ) +
      annotate(
        "text",
        x = .125,
        y = .89,
        label = sprintf("UNC = .%03d",
                        round(summary(r_selected)$uncertainty * 1000)),
        size = 2
      )
    
  } else if (annt_score_decom == "large"){
    p <- p +
      annotate(
        "text",
        x = .125,
        y = .96,
        label = sprintf("MCB = .%03d",
                        round(summary(r_selected)$miscalibration * 1000)),
        color = "red",
        size = 4
      ) +
      annotate(
        "text",
        x = .125,
        y = .90,
        label = sprintf("DSC = .%03d",
                        round(summary(r_selected)$discrimination * 1000)),
        size = 4
      ) +
      annotate(
        "text",
        x = .125,
        y = .84,
        label = sprintf("UNC = .%03d",
                        round(summary(r_selected)$uncertainty * 1000)),
        size = 4
      )
  } else if (annt_score_decom == "medium"){
    p <- p +
      annotate(
        "text",
        x = .125,
        y = .96,
        label = sprintf("MCB = .%03d",
                        round(summary(r_selected)$miscalibration * 1000)),
        color = "red",
        size = 3
      ) +
      annotate(
        "text",
        x = .125,
        y = .90,
        label = sprintf("DSC = .%03d",
                        round(summary(r_selected)$discrimination * 1000)),
        size = 3
      ) +
      annotate(
        "text",
        x = .125,
        y = .84,
        label = sprintf("UNC = .%03d",
                        round(summary(r_selected)$uncertainty * 1000)),
        size = 3
      )
  }
  
  
  return(p)
}

# -------------------------------------------------------------------
# Generate plots for selected models
# -------------------------------------------------------------------
# corp_plots_list_selected <- list()
# for (model_name in selected_models) {
#   corp_plots_list_selected[[model_name]] <- create_reliability_diag(
#     prev_peace_prob_month_long_binary_actual,
#     model_name
#   )
# }
corp_plots_list_selected_prev_peace <- list()
for (model_name in selected_models) {
  
  reliability_df <- prev_peace_prob_month_long_binary_actual %>%
    dplyr::filter(model == model_name)
  
  corp_plots_list_selected_prev_peace[[model_name]] <- reliabilitydiag.custom(
    reliability_df$predictive_probability,
    reliability_df$actual,
    selected_colors[model_name],
    0.9,
    "diagonal",  #"estimate" or "diagonal"
    "resampling",
    "small"
  ) + ggplot2::labs(
    title = model_labels[model_name],
    x = "Forecast value",
    y = "CEP"
  )
}


# -------------------------------------------------------------------
# Arrange them in a grid with a title
# -------------------------------------------------------------------
grid.arrange(
  tg,
  do.call(arrangeGrob, c(corp_plots_list_selected_prev_peace, ncol = 4)),
  ncol = 1,
  heights = c(0.1, 1)
)

if (store_plot == TRUE) {
  
  # build the arranged plot as a grob
  reliability_grid <- gridExtra::arrangeGrob(
    tg,
    do.call(arrangeGrob, c(corp_plots_list_selected_prev_peace, ncol = 4)),
    ncol = 1,
    heights = c(0.1, 1)
  )
  
  if (lower_fatalitiy_thresh == 0) {
    ggsave(
      "plotsOnsetDecember_fatalities_greq1/DecOnset_reliability_diagram_selected_models_prev_peace_1.png",
      plot = reliability_grid,
      width = 1.0 * 3000, height = 1.0 * 1800, dpi = 300, units = "px",
      bg = "white"
    )
  } else {
    ggsave(
      "plotsOnsetDecember_fatalities_greq25/DecOnset_reliability_diagram_selected_models_prev_peace_25.png",
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
  (corp_plots_list_selected_prev_peace$submission_final_omm +
     labs(title = "Reliability Diagram") +
     theme_fontsize) +
  plot_layout(heights = c(1, 1))


if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plotsOnsetDecember_fatalities_greq1/DecOnset_roc_reliability_selected_models_prev_peace_greq1.png",
           plot = roc_reliability_plot, width = 1.2 * 1000, height = 1.2 * 2000, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plotsOnsetDecember_fatalities_greq25/DecOnset_roc_reliability_selected_models_prev_peace_greq25.png",
           plot = roc_reliability_plot, width = 1.2 * 1000, height = 1.2 * 2000, dpi = 300, units = "px",
           bg="white")
  }
}



## -----------------------------------------------------------------------------
## BRIER: MSC-DSC-plots for previous peace
## -----------------------------------------------------------------------------

brier_decomposition_results <- prev_peace_prob_month_long_binary_actual %>%
  group_by(model) %>%
  group_split(.keep = TRUE) %>%
  set_names(map_chr(., ~ unique(.x$model))) %>%
  map(function(df) {
    res <- mcbdsc(df %>% select(predictive_probability),
                  y = df$actual,
                  score = "Brier_score") #'   One of: `"Brier_score"` (default), `"log_score"`, `"MR_score"`.
    
    estimates(res) %>%
      mutate(model = unique(df$model))
  })

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
  labs(title = "MSC-DSC Mean Brier Score Decomposition for Previous Peace") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.01)) +
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


if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plotsOnsetDecember_fatalities_greq1/DecOnset_brier_score_decomposition_selected_models_prev_peace_greq1.png",
           plot = brier_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plotsOnsetDecember_fatalities_greq25/DecOnset_brier_score_decomposition_selected_models_prev_peace_greq25.png",
           plot = brier_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px",
           bg="white")
  }
}

## ---
## All models
## ---

# 
# # combine to one dataframe
# brier_score_decomposition_all <- bind_rows(brier_decomposition_results) %>%
#   select(model, mean_score, MCB, DSC, UNC) %>%
#   rename(MeanScore = mean_score)
# 
# brier_score_decomposition_barplot_all <- brier_score_decomposition_all %>%
#   mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)
# 
# # data frame with format for the barchart
# brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
#   arrange(MeanScore) %>%  # sort by MeanScore
#   pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap, gap1, gap2),
#                names_to = "component",
#                values_to = "value") %>%
#   mutate(class = case_when(
#     component %in% c("MeanScore") ~ "SCORE",
#     component %in% c("MCB", "UNC") ~ "MCB_UNC",
#     component %in% c("DSC", "score_invisible") ~ "DSC_score", 
#     component %in% c("gap") ~ "GAP",
#     component %in% c("gap1") ~ "GAPmeanscoreDSC",
#     component %in% c("gap2") ~ "GAPdscUNC",
#   )) %>%
#   select(model, class, component, value)
# 
# brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
#   mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))
# 
# 
# # dataset for model orderbased on MeanScore
# brier_levs_all <- brier_score_decomposition_barplot_all %>%
#   filter(component == "MeanScore") %>%
#   arrange(desc(value)) %>%
#   pull(model)
# 
# # data for the plot
# brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
#   mutate(
#     class_group = case_when(
#       class == "MCB_UNC"   ~ 1,
#       class == "DSC_score"       ~ 2,
#       class == "SCORE" ~ 3,
#       class == "GAP" ~ 4,
#       class == "GAPmeanscoreDSC" ~ 5,
#       class == "GAPdscUNC" ~ 6
#       
#     ),  
#     model = factor(model, levels = brier_levs_all, ordered = TRUE)
#   )
# 
# brier_score_decomposition_barplot_all$model <- recode(
#   brier_score_decomposition_barplot_all$model,
#   !!!model_labels
# )
# 
# 
# # set order of subbars within group
# brier_score_decomposition_barplot_all$component <- factor(
#   brier_score_decomposition_barplot_all$component,
#   levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
# )
# 
# # set order of sub bars within group
# brier_score_decomposition_barplot_all$class_group <- factor(
#   brier_score_decomposition_barplot_all$class_group,
#   levels = c(4,1,6,2,3,5)
# )
# 
# # set width of mean score bar
# mean_score_bar_all <- 7
# # set width of mcb,dcs,unc bars
# msc_dsc_unc_bar_all <- 2
# 
# # width column
# brier_score_decomposition_barplot_all <- brier_score_decomposition_barplot_all %>%
#   mutate(
#     wdth = case_when(
#       # MeanScore
#       class_group == 3 ~ mean_score_bar_all, #4
#       # MCB DSC
#       class_group %in% c(1, 2) ~ msc_dsc_unc_bar_all,#msc_dsc_unc_bar_all, #1
#       # gap between groups top
#       class_group == 5 ~ 15,
#       # gapswithin groups
#       class_group == 6 ~ 2,
#       # gap between groups bottom
#       class_group == 4 ~ 10 #10
#       
#     )
#   )
# 
# 
# brier_score_decomposition_all <- ggplot(brier_score_decomposition_barplot_all, aes(x = class_group, y = value, fill = component)) +
#   geom_bar(stat = "identity", position = "stack", width = brier_score_decomposition_barplot_all$wdth) +
#   facet_grid(model ~ ., switch = "y") +
#   coord_flip() +
#   scale_fill_manual(
#     values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
#     ##
#     breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
#     ##
#     name = ""
#   ) +
#   labs(title = "MSC-DSC Mean Brier Score Decomposition") + 
#   theme_classic() +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   theme(
#     panel.spacing = unit(0, "points"),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text.y.left = element_text(angle = 0, hjust = 1, size = 15),
#     axis.text.y = element_blank(),
#     axis.ticks.length.y = unit(0, "points"),
#     axis.ticks.x = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
#     plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
#     legend.title = ggplot2::element_text(size = 15),
#     legend.text = ggplot2::element_text(size = 14),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.text.x = ggplot2::element_text(size = 14)
#   )
# 
# brier_score_decomposition_all

# ggsave("final_plots/brier_decomposition_all_models.png",
#        plot = brier_score_decomposition_all, width = 1.2 * 2297, height = 1.4 * 2181, dpi = 300, units = "px")


## -----------------------------------------------------------------------------
## LOG: MSC-DSC-plots
## -----------------------------------------------------------------------------
log_decomposition_results <- prev_peace_prob_month_long_binary_actual %>%
  group_by(model) %>%
  group_split(.keep = TRUE) %>%
  set_names(map_chr(., ~ unique(.x$model))) %>%
  map(function(df) {
    res <- mcbdsc(df %>% select(predictive_probability_log_nplustwo),
                  y = df$actual,
                  score = "log_score") #'   One of: `"Brier_score"` (default), `"log_score"`, `"MR_score"`.
    
    estimates(res) %>%
      mutate(model = unique(df$model))
  })


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


if(store_plot == TRUE){
  if(lower_fatalitiy_thresh == 0){
    ggsave("plotsOnsetDecember_fatalities_greq1/DecOnset_log_score_decomposition_selected_models_prev_peace_greq1.png",
           plot = log_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px",
           bg="white")
  } else {
    ggsave("plotsOnsetDecember_fatalities_greq25/DecOnset_log_score_decomposition_selected_models_prev_peace_greq25.png",
           plot = log_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px",
           bg="white")
  }
}

## ---
## All models
## ---
# 
# 
# # combine to one dataframe
# log_score_decomposition_all <- bind_rows(log_decomposition_results) %>%
#   select(model, mean_score, MCB, DSC, UNC) %>%
#   rename(MeanScore = mean_score)
# 
# log_score_decomposition_barplot_all <- log_score_decomposition_all %>%
#   mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)
# 
# # data frame with format for the barchart
# log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
#   arrange(MeanScore) %>%  # sort by MeanScore
#   pivot_longer(cols = c(MeanScore, MCB, DSC, UNC, score_invisible, gap, gap1, gap2),
#                names_to = "component",
#                values_to = "value") %>%
#   mutate(class = case_when(
#     component %in% c("MeanScore") ~ "SCORE",
#     component %in% c("MCB", "UNC") ~ "MCB_UNC",
#     component %in% c("DSC", "score_invisible") ~ "DSC_score", 
#     component %in% c("gap") ~ "GAP",
#     component %in% c("gap1") ~ "GAPmeanscoreDSC",
#     component %in% c("gap2") ~ "GAPdscUNC",
#   )) %>%
#   select(model, class, component, value)
# 
# log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
#   mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))
# 
# 
# # dataset for model orderbased on MeanScore
# log_levs_all <- log_score_decomposition_barplot_all %>%
#   filter(component == "MeanScore") %>%
#   arrange(desc(value)) %>%
#   pull(model)
# 
# # data for the plot
# log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
#   mutate(
#     class_group = case_when(
#       class == "MCB_UNC"   ~ 1,
#       class == "DSC_score"       ~ 2,
#       class == "SCORE" ~ 3,
#       class == "GAP" ~ 4,
#       class == "GAPmeanscoreDSC" ~ 5,
#       class == "GAPdscUNC" ~ 6
#       
#     ),  
#     model = factor(model, levels = log_levs_all, ordered = TRUE)
#   )
# 
# log_score_decomposition_barplot_all$model <- recode(
#   log_score_decomposition_barplot_all$model,
#   !!!model_labels
# )
# 
# 
# # set order of subbars within group
# log_score_decomposition_barplot_all$component <- factor(
#   log_score_decomposition_barplot_all$component,
#   levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
# )
# 
# # set order of sub bars within group
# log_score_decomposition_barplot_all$class_group <- factor(
#   log_score_decomposition_barplot_all$class_group,
#   levels = c(4,1,6,2,3,5)
# )
# 
# # set width of mean score bar
# mean_score_bar_all <- 7
# # set width of mcb,dcs,unc bars
# msc_dsc_unc_bar_all <- 2
# 
# # width column
# log_score_decomposition_barplot_all <- log_score_decomposition_barplot_all %>%
#   mutate(
#     wdth = case_when(
#       # MeanScore
#       class_group == 3 ~ mean_score_bar_all, #4
#       # MCB DSC
#       class_group %in% c(1, 2) ~ msc_dsc_unc_bar_all,#msc_dsc_unc_bar_all, #1
#       # gap between groups top
#       class_group == 5 ~ 15,
#       # gapswithin groups
#       class_group == 6 ~ 2,
#       # gap between groups bottom
#       class_group == 4 ~ 10 #10
#       
#     )
#   )
# 
# 
# log_score_decomposition_all <- ggplot(log_score_decomposition_barplot_all, aes(x = class_group, y = value, fill = component)) +
#   geom_bar(stat = "identity", position = "stack", width = log_score_decomposition_barplot_all$wdth) +
#   facet_grid(model ~ ., switch = "y") +
#   coord_flip() +
#   scale_fill_manual(
#     values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
#     ##
#     breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
#     ##
#     name = ""
#   ) +
#   labs(title = "MSC-DSC Mean log Score Decomposition") + 
#   theme_classic() +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   theme(
#     panel.spacing = unit(0, "points"),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text.y.left = element_text(angle = 0, hjust = 1, size = 15),
#     axis.text.y = element_blank(),
#     axis.ticks.length.y = unit(0, "points"),
#     axis.ticks.x = element_blank(),
#     axis.title = element_blank(),
#     legend.position = "bottom",
#     panel.grid.major.x = element_line(color = "#D9D9D9", size = 0.1),
#     plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
#     legend.title = ggplot2::element_text(size = 15),
#     legend.text = ggplot2::element_text(size = 14),
#     axis.line.x = element_blank(),
#     axis.line.y = element_blank(),
#     axis.text.x = ggplot2::element_text(size = 14)
#   )
# 
# log_score_decomposition_all

# ggsave("final_plots/log_score_decomposition_all_models.png",
#        plot = log_score_decomposition_all, width = 1.2 * 2297, height = 1.4 * 2181, dpi = 300, units = "px")



