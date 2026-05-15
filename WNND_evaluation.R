## clear environment
rm(list = ls())

## set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))



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
library(distfromq)
library(readxl)


# Erkenne das Betriebssystem
os <- Sys.info()["sysname"]

# Setze den Pfad abhängig vom Betriebssystem
# data_path <- "../Data/"
data_path <- ifelse(os == "Windows",
                    "//stat-meth-file1.stat.kit.edu/share-alle/Data/WNV forecasting challenge 2020/",  # Windows
                    "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/WNV forecasting challenge 2020/")  # macOS/Linux



## -----
## prediction dataset
## -----
# all files that begin with "2020-04-30"
prediction_files <- list.files(path = data_path, 
                               pattern = "^2020-04-30.*\\.csv$", 
                               full.names = TRUE)

# read all files and bind them into a single data frame
df_predictions <- do.call(rbind, lapply(prediction_files, read.csv))

df_predictions <- df_predictions %>%
  filter(type == "Bin")%>%
  select(-type, -unit)

## -----
## WNVND prediction competition participants
## -----
model_names <- c("ARS", "LANL", "MHC", "MSSM", "NCSU", "NYSW", 
                 "NYSW-CVD", "Rutgers", "Stanford", "UA", 
                 "UCD", "UI", "UI-NCSA", "WDH")

n_models <- length(model_names)


## -----
## FIPS county codes for easy identification and alignement of location names
## to the ones in the prediction dataframe
## -----

# FIPS codes from
# https://www.census.gov/geographies/reference-files/2020/demo/popest/2020-fips.html
geocodes <- read_excel(paste0(data_path, "all-geocodes-v2020.xlsx"), skip = 4)

geocodes <- geocodes %>%
  filter(`Summary Level` %in% c("040", "050"))


state_lookup <- geocodes %>%
  filter(`Summary Level` == "040") %>%
  select(`State Code (FIPS)`, location = `Area Name (including legal/statistical area description)`)

geocodes <- geocodes %>%
  left_join(state_lookup, by = "State Code (FIPS)") %>%
  filter(`Summary Level` == "050") %>%
  mutate(location = paste0(location, "-", `Area Name (including legal/statistical area description)`))

geocodes <- geocodes %>%
  mutate(location = str_replace_all(location, "city", "City"))

county_endings_to_delete <- "\\s+(County|Parish|Borough|Census Area|Municipality)$"

# delete endings in location that are in county_endings_to_delete
geocodes <- geocodes %>%
  mutate(location = str_remove(location, county_endings_to_delete))

# standardize naming convention from geocodes codebook to the one used in the prediction data
geocodes <- geocodes %>%
  mutate(location = str_replace(location, "Illinois-LaSalle", "Illinois-La Salle"),
         location = str_replace(location, "Indiana-DeKalb", "Indiana-De Kalb"),
         location = str_replace(location, "Indiana-LaPorte", "Indiana-La Porte"),
         location = str_replace(location, "Indiana-LaGrange", "Indiana-Lagrange"),
         location = str_replace(location, "Iowa-O'Brien", "Iowa-OBrien"),
         location = str_replace(location, "Louisiana-LaSalle", "Louisiana-La Salle"),
         location = str_replace(location, "Maryland-Prince George's", "Maryland-Prince Georges"),
         location = str_replace(location, "Maryland-Queen Anne's", "Maryland-Queen Annes"),
         location = str_replace(location, "Maryland-St. Mary's", "Maryland-St. Marys"),
         location = str_replace(location, "New Mexico-Doña Ana", "New Mexico-Dona Ana"),
         location = str_replace(location, "New Mexico-De Baca", "New Mexico-DeBaca"))

# uique list of the counties that predictions have been made for
all_locations_predictions <- unique(df_predictions$location)

# # test if all counties from the prediction data are in the geocodes dataset
# all_locations_geocodes <- unique(geocodes$location)
# missing_counties <- setdiff(all_locations_predictions, all_locations_geocodes)
# str_subset(all_locations_geocodes, "^Maryland-Q")


geocodes <- geocodes %>%
  filter(location %in% all_locations_predictions) %>%
  mutate(FIPS = paste0(`State Code (FIPS)`, `County Code (FIPS)`)) %>%
  mutate(FIPS = as.numeric(FIPS)) %>%
  select(location, FIPS)


## -----
## actual dataset
## -----
df_actual <- read.csv(paste0(data_path, "West Nile virus human and non-human activity by area of residence for year selected below_.csv"))

### Surveillance data are reported by county of residence, not the location (county or state) of exposure.
# https://www.cdc.gov/west-nile-virus/data-maps/historic-data.html

### achtung nur 480 counties im Datensatz: ist es üblich die zeros einfach nicht zu reporten?
df_actual_2020 <- df_actual %>%
  filter(Year == "2020") %>%
  select(-Activity, 
         -Total.human.disease.cases,
         -X..Presumptive.viremic.blood.donors, 
         -Notes)

df_actual_2020 <- df_actual_2020 %>%
  rename(FIPS = Location,
         actual = Neuroinvasive.disease.cases) %>%
  select(-FullGeoName)

df_actual_2020 <- left_join(df_actual_2020, geocodes, by = "FIPS")

new_actual_rows_for_counties_no_data_available <- geocodes %>%
  anti_join(df_actual_2020, by = "FIPS") %>%
  
  mutate(
    Year = 2020,
    actual = 0
  ) %>%
  
  select(Year, FIPS, actual, location)

df_actual_2020 <- bind_rows(df_actual_2020, new_actual_rows_for_counties_no_data_available)

# sort by FIPS
df_actual_2020 <- df_actual_2020 %>%
  arrange(FIPS)




## -----
## create list of predictive probabilites for all models
## -----

## join actual and predictions dataset
df_predictions <- left_join(df_predictions, df_actual_2020, by = "location") %>%
  select(-target)%>%
  relocate(Year, FIPS, .after = 1)


# list to store probabilistic predictions (in terms of CDFs)
predictive_probs <- list()

# upper bounds of the bins (quantiles)
qs <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 150, 200, 999)
# count values of the PDF
x_counts <- 0:999

# function to estimate PMF and CDF from binned forecasts
expand_prediction <- function(bin_probs) {
  
  if (bin_probs[1] >= 0.9999) {
    cdf_at_counts <- rep(1, length(x_counts))       # CDF is 1 for all counts
    prob_at_counts <- c(1, rep(0, length(x_counts) - 1)) # PMF is 1 at 0, 0 otherwise
    
    return(tibble(
      count = x_counts,
      PMF = prob_at_counts,
      CDF = cdf_at_counts
    ))
  }
  
  ### PRÜFEN OB DAS SINNVOLL IST #######################
  
  
  # 2. NEU: "Epsilon Smoothing" gegen die NaNs
  # Wir addieren eine extrem kleine Zahl, damit kein Wert jemals EXAKT 0 ist.
  bin_probs_safe <- bin_probs + 1e-10
  
  # Dadurch summiert sich alles auf etwas mehr als 1. 
  # Wir teilen durch die neue Summe, um es wieder perfekt auf 1.0 zu normieren.
  bin_probs_safe <- bin_probs_safe / sum(bin_probs_safe)
  
  ps <- cumsum(bin_probs_safe)
  
  ###############################################
  
  #ps <- cumsum(bin_probs)
  ps[ps > 1] <- 1 # CDF has to be <= 1, rounding errors
  
  # fit cdf via distfromq
  p_lognormal_approx <- make_p_fn(ps = ps,
                                  qs = qs,
                                  tail_dist = "lnorm")
  
  # get CDF at the count values
  cdf_at_counts <- p_lognormal_approx(x_counts)
  
  # PMF out of the step function CDF
  prob_at_counts <- c(cdf_at_counts[1], diff(cdf_at_counts))
  
  # return dataframe
  tibble(
    count = x_counts,
    PMF = prob_at_counts,
    CDF = cdf_at_counts
  )
}

# fill 
for (m in 1:length(model_names)) { 
  
  current_team <- model_names[m]
  
  df_model <- df_predictions %>%
    filter(team == current_team)
  
  # group by (location, FIPS, Year, actual)
  # isolate 15 rows (predictions per bin), 
  # apply 'expand_prediction' to column 'value' 
  df_expanded <- df_model %>%
    group_by(location, FIPS, Year, actual) %>%
    reframe(expand_prediction(value)) %>%
    ungroup()
  
  # order df
  df_expanded <- df_expanded %>%
    select(location, FIPS, count, Year, PMF, CDF, actual)
  
  predictive_probs[[current_team]] <- as.data.frame(df_expanded)
}




## -----
## scoring of the predictions via RPS, twRPS, Brier-Score
## -----
# define lower threshold for binary event of a present outbreak
lower_infections_thresh = 0



models_scoring_rules <- list()

for (m in names(predictive_probs)) {
  
  models_scoring_rules[[m]] <- predictive_probs[[m]] %>%
    group_by(FIPS, location, Year, actual) %>%
    # all new (score) calculations
    summarise(
      
      rps = sum((CDF - as.numeric(count >= actual))^2),
      
      twrps = sum((1 / (count + 1)) * (CDF - as.numeric(count >= actual))^2),
      
      onset_prob_pred = sum(PMF[count > 0]),
      
      actual_infection = first(actual) > lower_infections_thresh,
      
      brier_outbreak = (onset_prob_pred - actual_infection)^2,
      
      #########brier_outbreak_log_target,
      
      onset_prob_pred_nplustwo = (onset_prob_pred * 1000 + 1) / 1002,
      
      log_score_outbreak = ifelse(actual_infection == 1, 
                               log(onset_prob_pred_nplustwo), 
                               log(1 - onset_prob_pred_nplustwo)),
      
      .groups = "drop"
    ) %>%
    as.data.frame()
}


## -----
# label observations as either no infections ("none"), 
# WNVND outbreak ("onset"), ongoing infections ("ongoing") or 
# end of infections ("resolved")
# based on the WNVND cases from the previous year 2019
## -----
df_actual_2019 <- df_actual %>%
  filter(Year == "2019") %>%
  select(-Activity, 
         -Total.human.disease.cases,
         -X..Presumptive.viremic.blood.donors, 
         -Notes)

df_actual_2019 <- df_actual_2019 %>%
  rename(FIPS = Location,
         actual = Neuroinvasive.disease.cases) %>%
  select(-FullGeoName)

df_actual_2019 <- left_join(df_actual_2019, geocodes, by = "FIPS")

new_actual_rows_for_counties_no_data_available_2019 <- geocodes %>%
  anti_join(df_actual_2019, by = "FIPS") %>%
  
  mutate(
    Year = 2019,
    actual = 0
  ) %>%
  
  select(Year, FIPS, actual, location)

df_actual_2019 <- bind_rows(df_actual_2019, new_actual_rows_for_counties_no_data_available_2019)

# sort by FIPS
df_actual_2019 <- df_actual_2019 %>%
  arrange(FIPS)


wnv_situations <- df_actual_2020 %>%
  # rename actual to avoid confusion with the 2019 actuals
  select(FIPS, location, Year, actual_2020 = actual) %>%
  
  # join df_actual_2020 with .._2019 by FIPS
  left_join(
    df_actual_2019 %>% select(FIPS, actual_2019 = actual),
    by = "FIPS"
  ) %>%
  
  # if county wasnt present in 2019 set the actuals to 0
  mutate(actual_2019 = replace_na(actual_2019, 0)) %>%
  
  # assign situation labels
  mutate(
    outbreak_2020 = actual_2020 > 0,
    outbreak_2019 = actual_2019 > 0,
    
    situation = case_when(
      !outbreak_2020 & !outbreak_2019 ~ "none",
      outbreak_2020 & !outbreak_2019  ~ "onset",
      outbreak_2020 & outbreak_2019   ~ "ongoing",
      !outbreak_2020 & outbreak_2019  ~ "resolved"
    )
  ) %>%
  select(FIPS, location, Year, situation)









## -----
## Plot Data: RPS values by infectious situation
## -----
models_rps_infection <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "rps")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(models_rps_infection) <- c("FIPS", "location", "Year", model_names)
models_rps_infection <- list(models_rps_infection, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

n_total_counties <- nrow(models_rps_infection)

models_rps_infection <- models_rps_infection %>% select(!c("FIPS", "location", "Year")) %>% group_by(situation) %>% summarise_all(sum)

models_rps_infection[,2:ncol(models_rps_infection)] <- models_rps_infection[,2:ncol(models_rps_infection)] / n_total_counties # compute contributions to average CRPS


# create ggplot data frames
rps_year <- data.frame("RPS" = unlist(c(models_rps_infection[,2:ncol(models_rps_infection)])),
                       "Situation" = rep(models_rps_infection$situation, ncol(models_rps_infection)-1),
                       "Model" = rep(names(models_rps_infection)[2:ncol(models_rps_infection)], each = 4))

# # print RPS contribution of "none"
# print(colSums(models_rps_infection[1,2:ncol(models_rps_infection)] ))
# 
# # print overall CRPS per model
# print(colSums(models_rps_infection[1:nrow(models_rps_infection),2:ncol(models_rps_infection)] ))





## -----
## Plot Data: twRPS values by infectious situation
## -----
models_twrps_infection <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "twrps")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(models_twrps_infection) <- c("FIPS", "location", "Year", model_names)
models_twrps_infection <- list(models_twrps_infection, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

n_total_counties <- nrow(models_twrps_infection)

models_twrps_infection <- models_twrps_infection %>% select(!c("FIPS", "location", "Year")) %>% group_by(situation) %>% summarise_all(sum)

models_twrps_infection[,2:ncol(models_twrps_infection)] <- models_twrps_infection[,2:ncol(models_twrps_infection)] / n_total_counties # compute contributions to average CRPS


# create ggplot data frames
twrps_year <- data.frame("twRPS" = unlist(c(models_twrps_infection[,2:ncol(models_twrps_infection)])),
                       "Situation" = rep(models_twrps_infection$situation, ncol(models_twrps_infection)-1),
                       "Model" = rep(names(models_twrps_infection)[2:ncol(models_twrps_infection)], each = 4))



## -----
## Plot Data: Brier values by infectious situation
## -----
models_brier_infection <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "brier_outbreak")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(models_brier_infection) <- c("FIPS", "location", "Year", model_names)
models_brier_infection <- list(models_brier_infection, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

n_total_counties <- nrow(models_brier_infection)

models_brier_infection <- models_brier_infection %>% select(!c("FIPS", "location", "Year")) %>% group_by(situation) %>% summarise_all(sum)

models_brier_infection[,2:ncol(models_brier_infection)] <- models_brier_infection[,2:ncol(models_brier_infection)] / n_total_counties # compute contributions to average CRPS


# create ggplot data frames
brier_year <- data.frame("Brier" = unlist(c(models_brier_infection[,2:ncol(models_brier_infection)])),
                         "Situation" = rep(models_brier_infection$situation, ncol(models_brier_infection)-1),
                         "Model" = rep(names(models_brier_infection)[2:ncol(models_brier_infection)], each = 4))





## ---
## Data for onset prediction
## ---

# dataframe to store onset probabilities in long format
prob_models_all_situation_long <- data.frame(
  model = character(),
  Year = integer(),
  FIPS = integer(),
  onset_prob_pred = numeric(),
  situation = character()
)

for (model_name in names(models_scoring_rules)) {
  
  joined_data <- models_scoring_rules[[model_name]] %>%
    inner_join(wnv_situations, by = c("FIPS", "location", "Year")) %>%
    mutate(model = model_name) # add model name
  
  # add to dataframe
  prob_models_all_situation_long <- bind_rows(prob_models_all_situation_long, joined_data)
}



# filter for "onset" or "none" month (conflict this is previous peace)
prev_none_prob_month_long <- prob_models_all_situation_long %>%
  filter(situation == "none" | situation == "onset")

# onset actuals binary target
prev_none_prob_month_long_binary_actual <- prev_none_prob_month_long %>%
  mutate(actual = ifelse(actual > lower_infections_thresh, 1, 0))

















################################################################################
## PLOTS
################################################################################


## -----
## save plots in folders
## ----
store_plot <- FALSE

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





































































###################################################################################################################
### END
#################################################################################################################################
# 
# 
# 
# 
# # example
# bin_probs <- df_predictions$value[1:15]
# 
# # cdf (ps)
# ps <- cumsum(bin_probs)
# ps[ps > 1] <- 1 # CDF has to be <= 1, rounding errors
# 
# # fit cdf via distfromq
# p_lognormal_approx <- make_p_fn(ps = ps,
#                                 qs = qs,
#                                 tail_dist = "lnorm")
# 
# 
# # get CDF at the count values (equals step function CDF)
# cdf_at_counts <- p_lognormal_approx(x_counts)
# 
# # PMF out of the step function CDF
# prob_at_counts <- c(cdf_at_counts[1], diff(cdf_at_counts))
# 
# # save in dataframe
# df_counts <- data.frame(
#   Count = x_counts,
#   CDF = cdf_at_counts,       # P(X <= Count)
#   Probability = prob_at_counts # P(X == Count)
# )



















# 
# x <- seq(from = -0.5, to = 1000, length = 1000)
# 
# p_lognormal_approx <- make_p_fn(ps = ps,
#                                 qs = qs,
#                                 tail_dist = "lnorm")
# cdf_lognormal_approx <- p_lognormal_approx(x)
# 
# data.frame(
#   x = x,
#   y = cdf_lognormal_approx
# ) %>%
#   ggplot() +
#   geom_line(
#     mapping = aes(x = x, y = y),
#     size = 0.8
#   ) +
#   geom_point(
#     data = data.frame(q = qs, p = ps),
#     mapping = aes(x = q, y = p),
#     size = 1.2
#   ) +
#   scale_x_log10() + 
#   ylim(0, 1) +
#   ylab("Probability") +
#   xlab("") +
#   theme_bw() + 
#   xlim(0,20)
# 
# 
# 












