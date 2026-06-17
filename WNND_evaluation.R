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
library(patchwork)
library(sf)
library(tigris)
library(scales)
library(ggbreak)


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
  ##########################
  ## ebenfalls prüfen
  ##ps[length(ps)] <- 1 ------ sind alle ps davor kleiner 1??
  ##ps <- ps / ps[length(ps)]------ alternativ?
  
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
# for (m in 1:length(model_names)) { 
#   
#   current_team <- model_names[m]
#   
#   # Wir lassen R mit dir sprechen, damit du weißt, wo es steht!
#   cat(sprintf("[%s] Verarbeite Modell %d von %d: %s ...\n", 
#               Sys.time(), m, length(model_names), current_team))
#   
#   
#   df_model <- df_predictions %>%
#     filter(team == current_team)
#   
#   if (nrow(df_model) == 0) {
#     cat("-> Keine Daten für dieses Modell. Überspringe...\n")
#     next
#   }
#   
#   # group by (location, FIPS, Year, actual)
#   # isolate 15 rows (predictions per bin), 
#   # apply 'expand_prediction' to column 'value' 
#   df_expanded <- df_model %>%
#     group_by(location, FIPS, Year, actual) %>%
#     reframe(expand_prediction(value)) %>%
#     ungroup()
#   
#   # order df
#   df_expanded <- df_expanded %>%
#     select(location, FIPS, count, Year, PMF, CDF, actual)
#   
#   predictive_probs[[current_team]] <- as.data.frame(df_expanded)
# }
# save(predictive_probs, file = "output/predictive_provs_WNVND.RData")
load("output/predictive_provs_WNVND.RData")




# df_all_predictions <- bind_rows(predictive_probs, .id = "model")
# write_parquet(df_all_predictions, paste0("output/", "all_predictive_probs_expanded.parquet"))
# df_all_predictions <- read_parquet(paste0("output/", "all_predictive_probs_expanded.parquet"))
# predictive_probs <- split(df_all_predictions, df_all_predictions$model)
# predictive_probs <- lapply(predictive_probs, function(df) select(df, -model))


## -----
## scoring of the predictions via crps, twcrps, Brier-Score
## -----
# define lower threshold for binary event of a present outbreak
lower_infections_thresh = 0



models_scoring_rules <- list()

for (m in names(predictive_probs)) {
  
  models_scoring_rules[[m]] <- predictive_probs[[m]] %>%
    group_by(FIPS, location, Year, actual) %>%
    # all new (score) calculations
    summarise(
      
      crps = sum((CDF - as.numeric(count >= actual))^2),
      
      twcrps = sum((1 / (count + 1)) * (CDF - as.numeric(count >= actual))^2),
      
      onset_prob_pred = sum(PMF[count > 0]),
      
      actual_infection = first(actual) > lower_infections_thresh,
      
      brier_outbreak = (onset_prob_pred - actual_infection)^2,
      
      #### überprüfen
      onset_prob_pred_nplustwo = (onset_prob_pred * 1000 + 1) / 1002,
      #######
      ### clipping?
      # # Wahrscheinlichkeit auf [1e-15, 1 - 1e-15] begrenzen
      # onset_prob_pred_safe = pmax(pmin(onset_prob_pred, 1 - 1e-15), 1e-15),
      # 
      # log_score_outbreak = ifelse(actual_infection == 1, 
      #                             log(onset_prob_pred_safe), 
      #                             log(1 - onset_prob_pred_safe))
      
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
## Plot Data: crps values by infectious situation
## -----
models_crps_infection <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "crps")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(models_crps_infection) <- c("FIPS", "location", "Year", model_names)
models_crps_infection <- list(models_crps_infection, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

n_total_counties <- nrow(models_crps_infection)

models_crps_infection <- models_crps_infection %>% select(!c("FIPS", "location", "Year")) %>% group_by(situation) %>% summarise_all(sum)

models_crps_infection[,2:ncol(models_crps_infection)] <- models_crps_infection[,2:ncol(models_crps_infection)] / n_total_counties # compute contributions to average crps


# create ggplot data frames
crps_year <- data.frame("CRPS" = unlist(c(models_crps_infection[,2:ncol(models_crps_infection)])),
                       "Situation" = rep(models_crps_infection$situation, ncol(models_crps_infection)-1),
                       "Model" = rep(names(models_crps_infection)[2:ncol(models_crps_infection)], each = 4))

# # print crps contribution of "none"
# print(colSums(models_crps_infection[1,2:ncol(models_crps_infection)] ))
# 
# # print overall crps per model
# print(colSums(models_crps_infection[1:nrow(models_crps_infection),2:ncol(models_crps_infection)] ))





## -----
## Plot Data: twcrps values by infectious situation
## -----
models_twcrps_infection <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "twcrps")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(models_twcrps_infection) <- c("FIPS", "location", "Year", model_names)
models_twcrps_infection <- list(models_twcrps_infection, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

n_total_counties <- nrow(models_twcrps_infection)

models_twcrps_infection <- models_twcrps_infection %>% select(!c("FIPS", "location", "Year")) %>% group_by(situation) %>% summarise_all(sum)

models_twcrps_infection[,2:ncol(models_twcrps_infection)] <- models_twcrps_infection[,2:ncol(models_twcrps_infection)] / n_total_counties # compute contributions to average crps


# create ggplot data frames
twcrps_year <- data.frame("twCRPS" = unlist(c(models_twcrps_infection[,2:ncol(models_twcrps_infection)])),
                       "Situation" = rep(models_twcrps_infection$situation, ncol(models_twcrps_infection)-1),
                       "Model" = rep(names(models_twcrps_infection)[2:ncol(models_twcrps_infection)], each = 4))



## -----
## Plot Data: Brier values by infectious situation
## -----
models_brier_infection <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "brier_outbreak")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(models_brier_infection) <- c("FIPS", "location", "Year", model_names)
models_brier_infection <- list(models_brier_infection, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

n_total_counties <- nrow(models_brier_infection)

models_brier_infection <- models_brier_infection %>% select(!c("FIPS", "location", "Year")) %>% group_by(situation) %>% summarise_all(sum)

models_brier_infection[,2:ncol(models_brier_infection)] <- models_brier_infection[,2:ncol(models_brier_infection)] / n_total_counties # compute contributions to average crps


# create ggplot data frames
brier_year <- data.frame("BRIER" = unlist(c(models_brier_infection[,2:ncol(models_brier_infection)])),
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

# actuals binary target
prob_models_all_situation_long_binary_actual <- prob_models_all_situation_long %>%
  mutate(actual = ifelse(actual > lower_infections_thresh, 1, 0))








################################################################################
## PLOTS
################################################################################


## -----
## save plots in folders
## ----
store_plot <- TRUE
# folder to store plots
folder <- "plots_infections"

## -----
## labels, colors, textsize etc.
## ----

model_labels <- c(
  "ARS" = "ARS",
  "LANL" = "LANL",
  "MHC" = "MHC",
  "MSSM" = "MSSM",
  "NCSU" = "NCSU",
  "NYSW" = "NYSW",
  "NYSW-CVD" = "NYSW-CVD",
  "Rutgers" = "Rutgers",
  "Stanford" = "Stanford",
  "UA" = "UA",
  "UCD" = "UCD",
  "UI" = "Naive",
  "UI-NCSA" = "UI-NCSA",
  "WDH" = "WDH"
  
)

model_labels_df <- data.frame("name_original" = names(model_labels), "name_paper" = model_labels)

## -----
## Selected models
## -----

## --
## change this if needed. everything else is dynamic!
selected_models <- c("ARS", "LANL",
                     "MHC", "MSSM",
                     "UI", "NYSW",
                     "Stanford", "WDH")
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
  axis.title = ggplot2::element_text(size = 14),
  axis.text = ggplot2::element_text(size = 12),
  legend.title = element_text(size = 14),
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



## labels for situations
onset_label <- "Introduction"
ongoing_label <- "Ongoing Infections"
resolved_label <- "End of Infections"
none_label <- "Ongoing Absence"


##################
## Figure 1a: cases per county and outbreak situation 2020
##################
# df with actuals and outbreak situations
actual_2020_situations <- wnv_situations %>%
  left_join(df_actual_2020, by = c("FIPS","location","Year"))

# load us counties: cartographic boundaries true (shoreline), resolution 1:20 million
us_counties <- counties(cb = TRUE, resolution = "20m", year = 2020) %>%
  # filter s.t. contiguous US is obtained (Alaska, Hawaii, Samoa, Guam, Mariana Islands, Puerto Rico, Virgin Islands)
  filter(!(STATEFP %in% c("02", "15", "60", "66", "69", "72", "78")))

# 2. df for plot
map_data <- actual_2020_situations %>%
  # add leading zero if missing
  mutate(GEOID = sprintf("%05d", as.numeric(FIPS))) %>%
  right_join(us_counties, by = "GEOID") %>%
  st_as_sf()# geographic map

hotspot_centroids <- suppressWarnings(
  map_data %>%
    filter(situation %in% c("onset", "ongoing", "resolved")) %>%
    st_centroid() # calculate the centorid of each county that is not "none"
)
hotspot_points <- cbind(hotspot_centroids, st_coordinates(hotspot_centroids)) %>%
  st_drop_geometry()


wnvnd_map_plot <- ggplot() +
  
  # map of US
  geom_sf(data = map_data, aes(fill = actual), color = "white", linewidth = 0.1) +
  
  # shape and color of the points
  geom_point(data = hotspot_points,
             aes(x = X, y = Y, shape = situation, color = situation),
             size = 1.8) +
  
  # Shapes of the midpoints of the counties
  scale_shape_manual(
    name = "Disease Situation", 
    values = c(
      "onset" = 17,       # triangle
      "ongoing" = 15,     # square
      "resolved" = 16     # circle
    ),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label, 
      "resolved" = resolved_label
    )
  ) +
  
  # Colors of the midpoints of the counties
  scale_color_manual(
    name = "Disease Situation",
    values = c(
      "onset" = "black",
      "ongoing" = "grey30",
      "resolved" = "grey"
    ),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label, 
      "resolved" = resolved_label
    )
  ) +
  
  scale_fill_gradientn(
    colors = c("#e0f3f8", "#fdae61", "#f46d43", "#d73027", "#a50026"), 
    trans = "log1p",
    name = "Cases",
    breaks = c(0, 1, 5, 10, 50, 100),
    na.value = "grey85"
  ) +
  
  theme_void() + 
  labs(
    title = "WNND Cases 2020"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 10)),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    
    legend.spacing.y = unit(0.5, "cm")
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  )

plot(wnvnd_map_plot)


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","WNND_map.png"),
         plot = wnvnd_map_plot, width = 1.0 * 4222, height = 1.5 * 1300, dpi = 300, units = "px",
         bg="white")
}




# check wether actuals coincide with the analysis of the paper on page 5
# 559 WNND cases?
sum(actual_2020_situations$actual)
test <- actual_2020_situations %>%
  filter(actual > 0)

# 181 counties
length(test$actual)


##################
## Figure 1b: distribution of the actual cases
##################
# Main-Plot: only counties cases>0
hist_main <- actual_2020_situations %>%
  filter(situation %in% c("onset", "ongoing")) %>%
  ggplot(aes(x = actual, fill = situation)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.9) +
  scale_fill_manual(
    name = "Outbreak Status",
    values = c("onset" = "#4664aa", 
               "ongoing" = "#a22223"),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label
    )
  ) +
  labs(
    title = "Distribution of WNND Cases 2020",
    x = "WNND Cases",
    y = "Counties"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# 2. Inset-plot: all counties
hist_inset <- ggplot(actual_2020_situations, aes(x = actual)) +
  geom_histogram(binwidth = 1, fill = "#74add1", color = "white", alpha = 0.9) +
  labs(
    title = "All US Counties",
    x = "Cases",
    y = "Counties (log)"
  ) +
  scale_y_continuous(
    breaks = c(0, 10, 100, 1000, 3000),
    labels = label_comma()
  ) +
  coord_transform(y = "log1p") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  )

# Picture in picture
hist_combined <- hist_main + 
  inset_element(hist_inset, left = 0.5, bottom = 0.5, right = 0.98, top = 0.98)

plot(hist_combined)

if(store_plot == TRUE){
  ggsave(paste0(folder,"/","WNND_cases_distribution.png"),
         plot = hist_combined, width = 0.9 * 3500, height = 0.8 * 2000, dpi = 300, units = "px",
         bg="white")
}




## -----------------------------------------------------------------------------
## Figure 2a: Crps per county and outbreak situation
## -----------------------------------------------------------------------------
crps_per_county <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "crps")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(crps_per_county) <- c("FIPS", "location", "Year", model_names)
crps_per_county <- list(crps_per_county, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

crps_per_county <- crps_per_county %>%
  mutate(sum_CRPS = rowSums(across(4:17), na.rm = TRUE)) %>%
  select(-c(4:17))
  
map_data_crps <- crps_per_county %>%
  # add leading zero if missing
  mutate(GEOID = sprintf("%05d", as.numeric(FIPS))) %>%
  right_join(us_counties, by = "GEOID") %>%
  st_as_sf()# geographic map

hotspot_centroids_crps <- suppressWarnings(
  map_data_crps %>%
    filter(situation %in% c("onset", "ongoing", "resolved")) %>%
    st_centroid() # calculate the centorid of each county that is not "none"
)
hotspot_points_crps <- cbind(hotspot_centroids_crps, st_coordinates(hotspot_centroids_crps)) %>%
  st_drop_geometry()


wnvnd_crps_map_plot <- ggplot() +
  
  # map of US
  geom_sf(data = map_data_crps, aes(fill = sum_CRPS), color = "white", linewidth = 0.1) +
  
  # shape and color of the points
  geom_point(data = hotspot_points_crps,
             aes(x = X, y = Y, shape = situation, color = situation),
             size = 1.8) +
  
  # Shapes of the midpoints of the counties
  scale_shape_manual(
    name = "Disease Situation", 
    values = c(
      "onset" = 17,       # triangle
      "ongoing" = 15,     # square
      "resolved" = 16     # circle
    ),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label, 
      "resolved" = resolved_label
    )
  ) +
  
  # Colors of the midpoints of the counties
  scale_color_manual(
    name = "Disease Situation",
    values = c(
      "onset" = "black",
      "ongoing" = "grey30",
      "resolved" = "grey"
    ),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label, 
      "resolved" = resolved_label
    )
  ) +
  
  scale_fill_gradientn(
    colors = c("#e0f3f8", "#fdae61", "#f46d43", "#d73027", "#a50026"), 
    trans = "log1p",
    name = "CRPS",
    breaks = c(30, 50, 200, 400, 600),
    na.value = "white"
  ) +
  
  theme_void() + 
  labs(
    title = "Sum of CRPS Across All Models for WNND Cases 2020"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 10)),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    
    legend.spacing.y = unit(0.5, "cm")
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  )

plot(wnvnd_crps_map_plot)


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","WNND_crps_map.png"),
         plot = wnvnd_crps_map_plot, width = 1.0 * 4222, height = 1.5 * 1300, dpi = 300, units = "px",
         bg="white")
}


## -----------------------------------------------------------------------------
## Figure 2b: Crps decomposition
## -----------------------------------------------------------------------------
crps_selected_models <- crps_year %>% filter(Model %in% selected_models) %>%
  mutate("Model_orig" = Model)

# Rename models
crps_selected_models$Model <- recode(
  crps_selected_models$Model,
  !!!model_labels
)

crps_infections_situation_plot <- crps_selected_models %>%
  group_by(Model) %>%
  summarise(total_crps = sum(CRPS), .groups = "drop") %>%
  right_join(crps_selected_models, by = "Model") %>%
  mutate(Model = fct_reorder(Model, total_crps, .desc = TRUE),
         Situation = factor(Situation, levels = c("ongoing", "resolved", "onset", "none"))) %>%
  ggplot(aes(fill = Situation, y = Model, x = CRPS)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean CRPS by Disease Situation",
       x = "Mean CRPS") +
  scale_fill_manual("Disease Situation",
                    values = c("ongoing" = "#a22223",
                               "resolved" = "#d09191",
                               "onset" = "#4664aa",
                               "none" = "#a2b2d4"),
                    labels = c("ongoing" = ongoing_label,
                               "resolved" = resolved_label,
                               "onset" = onset_label,
                               "none" = none_label)) +
  scale_x_break(c(1.5, 22.5), scales = "fixed") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  theme_fontsize

crps_infections_situation_plot

if(store_plot == TRUE){
  ggsave(paste0(folder,"/","crps_infections_situation.png"),
         plot = crps_infections_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
         bg="white")
}




## -----
## counties per category
## -----
# wichtig: CRPS dominated by high infectious counties?
# auf plot erstmal nicht erkennbar, da rote zu blaue balken ausgeglichen
# wenn die anzahl an counties mit no previous cases aber sehr gering ist
# spricht das trotzdem für die betonung des CRPS auf einige wenige Counties
# mit hohen Fallzahlen
counties_per_situation <- wnv_situations %>%
  summarise(onset = sum(situation == "onset", na.rm = TRUE),
            ongoing = sum(situation == "ongoing", na.rm = TRUE),
            resolved = sum(situation == "resolved", na.rm = TRUE),
            none = sum(situation == "none", na.rm = TRUE))
  

## -----------------------------------------------------------------------------
## Figure 2c: twCrps per county and outbreak situation
## -----------------------------------------------------------------------------
twcrps_per_county <- lapply(models_scoring_rules, function(m) m %>% select("FIPS", "location", "Year", "twcrps")) %>%
  reduce(left_join, c("FIPS", "location", "Year"))
names(twcrps_per_county) <- c("FIPS", "location", "Year", model_names)
twcrps_per_county <- list(twcrps_per_county, wnv_situations) %>% reduce(left_join, c("FIPS", "location", "Year"))

twcrps_per_county <- twcrps_per_county %>%
  mutate(sum_twCRPS = rowSums(across(4:17), na.rm = TRUE)) %>%
  select(-c(4:17))

map_data_twcrps <- twcrps_per_county %>%
  # add leading zero if missing
  mutate(GEOID = sprintf("%05d", as.numeric(FIPS))) %>%
  right_join(us_counties, by = "GEOID") %>%
  st_as_sf()# geographic map

hotspot_centroids_twcrps <- suppressWarnings(
  map_data_twcrps %>%
    filter(situation %in% c("onset", "ongoing", "resolved")) %>%
    st_centroid() # calculate the centorid of each county that is not "none"
)
hotspot_points_twcrps <- cbind(hotspot_centroids_twcrps, st_coordinates(hotspot_centroids_twcrps)) %>%
  st_drop_geometry()


wnvnd_twcrps_map_plot <- ggplot() +
  
  # map of US
  geom_sf(data = map_data_twcrps, aes(fill = sum_twCRPS), color = "white", linewidth = 0.1) +
  
  # shape and color of the points
  geom_point(data = hotspot_points_twcrps,
             aes(x = X, y = Y, shape = situation, color = situation),
             size = 1.8) +
  
  # Shapes of the midpoints of the counties
  scale_shape_manual(
    name = "Disease Situation", 
    values = c(
      "onset" = 17,       # triangle
      "ongoing" = 15,     # square
      "resolved" = 16     # circle
    ),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label, 
      "resolved" = resolved_label
    )
  ) +
  
  # Colors of the midpoints of the counties
  scale_color_manual(
    name = "Disease Situation",
    values = c(
      "onset" = "black",
      "ongoing" = "grey30",
      "resolved" = "grey"
    ),
    labels = c(
      "onset" = onset_label, 
      "ongoing" = ongoing_label, 
      "resolved" = resolved_label
    )
  ) +
  
  scale_fill_gradientn(
    colors = c("#e0f3f8", "#fdae61", "#f46d43", "#d73027", "#a50026"), 
    trans = "log1p",
    name = "twCRPS",
    breaks = c(5, 10, 20, 40, 60),
    na.value = "white"
  ) +
  
  theme_void() + 
  labs(
    title = "Sum of twCRPS Across All Models for WNND Cases 2020"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 18, hjust = 0.5, margin = margin(b = 10)),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    
    legend.spacing.y = unit(0.5, "cm")
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  )

plot(wnvnd_twcrps_map_plot)


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","WNND_twcrps_map.png"),
         plot = wnvnd_twcrps_map_plot, width = 1.0 * 4222, height = 1.5 * 1300, dpi = 300, units = "px",
         bg="white")
}


## -----------------------------------------------------------------------------
## Figure 2c: twCrps decomposition
## -----------------------------------------------------------------------------
## ---
## Selected models
## ---
twcrps_selected_models <- twcrps_year %>% filter(Model %in% selected_models) %>%
  mutate("Model_orig" = Model)

# Rename models
twcrps_selected_models$Model <- recode(
  twcrps_selected_models$Model,
  !!!model_labels
)

twcrps_infections_situation_plot <- twcrps_selected_models %>%
  group_by(Model) %>%
  summarise(total_twcrps = sum(twCRPS), .groups = "drop") %>%
  right_join(twcrps_selected_models, by = "Model") %>%
  mutate(Model = fct_reorder(Model, total_twcrps, .desc = TRUE),
         Situation = factor(Situation, levels = c("ongoing", "resolved", "onset", "none"))) %>%
  ggplot(aes(fill = Situation, y = Model, x = twCRPS)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean twCRPS by Disease Situation",
       x = "Mean twCRPS") +
  scale_fill_manual("Disease Situation",
                    values = c("ongoing" = "#a22223",
                               "resolved" = "#d09191",
                               "onset" = "#4664aa",
                               "none" = "#a2b2d4"),
                    labels = c("ongoing" = ongoing_label,
                               "resolved" = resolved_label,
                               "onset" = onset_label,
                               "none" = none_label)) +
  scale_x_break(c(0.4, 2.7), scales = "fixed") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  theme_fontsize

twcrps_infections_situation_plot

if(store_plot == TRUE){
  ggsave(paste0(folder,"/","twcrps_infections_situation.png"),
         plot = twcrps_infections_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
         bg="white")
}



## -----------------------------------------------------------------------------
## Figure 2d: Brier Score decomposition
## -----------------------------------------------------------------------------
## ---
## Selected models
## ---
brier_selected_models <- brier_year %>% filter(Model %in% selected_models) %>%
  mutate("Model_orig" = Model)

# Rename models
brier_selected_models$Model <- recode(
  brier_selected_models$Model,
  !!!model_labels
)

brier_infections_situation_plot <- brier_selected_models %>%
  group_by(Model) %>%
  summarise(total_brier = sum(BRIER), .groups = "drop") %>%
  right_join(brier_selected_models, by = "Model") %>%
  mutate(Model = fct_reorder(Model, total_brier, .desc = TRUE),
         Situation = factor(Situation, levels = c("ongoing", "resolved", "onset", "none"))) %>%
  ggplot(aes(fill = Situation, y = Model, x = BRIER)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Contribution to the Mean Brier-Score by Disease Situation",
       x = "Mean Brier-Score") +
  scale_fill_manual("Disease Situation",
                    values = c("ongoing" = "#a22223",
                               "resolved" = "#d09191",
                               "onset" = "#4664aa",
                               "none" = "#a2b2d4"),
                    labels = c("ongoing" = ongoing_label,
                               "resolved" = resolved_label,
                               "onset" = onset_label,
                               "none" = none_label)) +
  scale_x_break(c(0.2, 0.65), scales = "fixed") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  ) +
  theme_fontsize

brier_infections_situation_plot

if(store_plot == TRUE){
  ggsave(paste0(folder,"/","brier_infections_situation.png"),
         plot = brier_infections_situation_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
         bg="white")
}






## -----------------------------------------------------------------------------
## Figure 2e: Brier score per threshold => Crps
## -----------------------------------------------------------------------------
all_models_brier <- map_dfr(names(predictive_probs), function(m) {
  
  predictive_probs[[m]] %>%
    # calculate mean Brier-Score for each count and County
    mutate(
      # does event of threshold a being greq outcome happen?
      indicator = as.numeric(count >= actual),
      # Brier score at threshold (count)
      brier_score_k = (CDF - indicator)^2
    ) %>%
    
    # group per count (i.e. threshold of the CRPS)
    group_by(count) %>%
    
    # avg Brier over all counties
    summarise(
      brier_avg = mean(brier_score_k),
      .groups = "drop"
    ) %>%
    
    # rename columns for following plots
    mutate(
      model = m,
      a = count,
      log_a_plus1 = log(a + 1)
    ) %>%
    
    select(model, a, log_a_plus1, brier_avg)
})

selected_models_brier <- all_models_brier %>%
  filter(model %in% selected_models)


crps_brier_integrands_plot <- ggplot(selected_models_brier, aes(x = a, y = brier_avg, colour = model)) +
  geom_line(size = .95, alpha = 0.8) + 
  labs(
    title = "Mean Brier Score by Threshold (CRPS Integrand)",
    x = "Threshold a",
    y = "Mean Brier Score",
    colour = ""
  ) +
  xlim(0,200) + 
  scale_color_manual(
    values = selected_model_colors,
    breaks = names(selected_model_labels),
    labels = selected_model_labels
  ) +
  theme_bw() +
  theme_fontsize

crps_brier_integrands_plot


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","crps_brier_integrands_plot.png"),
         plot = crps_brier_integrands_plot, width = 1.0 * 4222, height = 1.0 * 1300, dpi = 300, units = "px",
         bg="white")
}


##
# check if the sum(meanBS) = meancrps
##
# crps_by_model <- crps_selected_models %>%
#   group_by(Model) %>%
#   summarise(total_crps = sum(CRPS), .groups = "drop") %>%
#   right_join(crps_selected_models, by = "Model") %>% 
#   distinct(Model_orig, total_crps)
# 
# crps_brierSUM_results <- selected_models_brier %>%
#   group_by(model) %>%
#   summarise(CRPSBrierapprox = sum(brier_avg))
# 
# crps_brierSUM_results <- crps_brierSUM_results %>%
#   left_join(crps_by_model,
#             by = c("model" = "Model_orig")) %>%
#   mutate(brier_proportion = CRPSBrierapprox/total_crps)
# 
# print(crps_brierSUM_results)




## -----------------------------------------------------------------------------
## Figure 2f: Brier score per log-change threshold => twCRPS
## -----------------------------------------------------------------------------
twcrps_brier_integrands_plot <- ggplot(selected_models_brier, aes(x = log_a_plus1, y = brier_avg, colour = model)) +
  geom_line(size = .95, alpha = 0.8) +
  labs(
    title = "Mean Brier Score by Log-Change Threshold (twCRPS Integrand)",
    x = "Threshold log(a+1)",
    y = "Mean Brier Score",
    colour = ""
  ) +
  scale_color_manual(
    values = selected_model_colors,
    breaks = names(selected_model_labels),
    labels = selected_model_labels
  ) +
  theme_bw() +
  theme_fontsize

twcrps_brier_integrands_plot


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","twcrps_brier_integrands_plot.png"),
         plot = twcrps_brier_integrands_plot, width = 1.0 * 3200, height = 1.0 * 1300, dpi = 300, units = "px",
         bg="white")
}


## -----------------------------------------------------------------------------
## Figure 3: ROC curves for previous peace
## -----------------------------------------------------------------------------

##
## https://cran.r-project.org/web/packages/precrec/vignettes/introduction.html
##
roc_data_prev_none <- prev_none_prob_month_long_binary_actual

# used to check wether AUC matches AUC in Appendix Fig S6 (matches!)
#roc_data_prev_none <- prob_models_all_situation_long_binary_actual

## ---
## In-depth: 8 models
## ---
# filter relevant models
roc_data_prev_none_selected <- roc_data_prev_none %>%
  filter(model %in% selected_models)

# create list of scores
score_list_selected <- lapply(selected_models, function(m) {
  roc_data_prev_none_selected %>%
    filter(model == m) %>%
    pull(onset_prob_pred)
})

# create list of labels
label_list_selected <- lapply(selected_models, function(m) {
  roc_data_prev_none_selected %>%
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
  ggsave(paste0(folder,"/","roc_curve_selected_models_prev_none.png"),
         plot = roc_curve_selected_models, width = 1.2 * 3391, height = 1.2 * 1225, dpi = 300, units = "px",
         bg="white")
}



## -----------------------------------------------------------------------------
## Figure 4: Reliability diagrams for previous peace
## -----------------------------------------------------------------------------
# Title grob for the combined plot
tg <- textGrob("CORP Reliability Diagrams", gp = gpar(fontsize = 18, hjust = 0.5))


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
corp_plots_list_selected_no_prev_cases <- list()
for (model_name in selected_models) {
  
  reliability_df <- prev_none_prob_month_long_binary_actual %>%
    dplyr::filter(model == model_name)
  
  corp_plots_list_selected_no_prev_cases[[model_name]] <- reliabilitydiag.custom(
    reliability_df$onset_prob_pred,
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
  do.call(arrangeGrob, c(corp_plots_list_selected_no_prev_cases, ncol = 4)),
  ncol = 1,
  heights = c(0.1, 1)
)

# build the arranged plot as a grob
reliability_grid <- gridExtra::arrangeGrob(
  tg,
  do.call(arrangeGrob, c(corp_plots_list_selected_no_prev_cases, ncol = 4)),
  ncol = 1,
  heights = c(0.1, 1)
)

plot(reliability_grid)

if (store_plot == TRUE) {
  ggsave(
    paste0(folder,"/","reliability_diagram_selected_models_no_prev_cases.png"),
    plot = reliability_grid,
    width = 1.0 * 3000, height = 1.0 * 1800, dpi = 300, units = "px",
    bg = "white"
  )
}

# -------------------------------------------------------------------
# Align ROC curve
# -------------------------------------------------------------------

roc_reliability_plot <-
  (roc_curve_selected_models + theme(legend.position = "none")) /
  (corp_plots_list_selected_no_prev_cases$ARS +
     labs(title = "Reliability Diagram") +
     theme_fontsize) +
  plot_layout(heights = c(1, 1))

plot(roc_reliability_plot)


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","roc_reliability_selected_models_no_prev_cases.png"),
         plot = roc_reliability_plot, width = 1.2 * 1000, height = 1.2 * 2000, dpi = 300, units = "px",
         bg="white")
}








## -----------------------------------------------------------------------------
## BRIER: MSC-DSC-plots for previous peace
## -----------------------------------------------------------------------------
brier_decomposition_results <- prev_none_prob_month_long_binary_actual %>%
  group_by(model) %>%
  group_split(.keep = TRUE) %>%
  set_names(map_chr(., ~ unique(.x$model))) %>%
  map(function(df) {
    res <- mcbdsc(df %>% select(onset_prob_pred),
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
  labs(title = "Mean Brier Score Decomposition") + 
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1)) +
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
  ggsave(paste0(folder,"/","brier_score_decomposition_selected_models_no_prev_cases.png"),
         plot = brier_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px",
         bg="white")
}



# 
# ## -----------------------------------------------------------------------------
# ## LOG: MSC-DSC-plots
# ## -----------------------------------------------------------------------------
# log_decomposition_results <- prev_none_prob_month_long_binary_actual %>%
#   group_by(model) %>%
#   group_split(.keep = TRUE) %>%
#   set_names(map_chr(., ~ unique(.x$model))) %>%
#   map(function(df) {
#     res <- mcbdsc(df %>% select(onset_prob_pred_nplustwo),
#                   y = df$actual,
#                   score = "log_score") #'   One of: `"Brier_score"` (default), `"log_score"`, `"MR_score"`.
#     
#     estimates(res) %>%
#       mutate(model = unique(df$model))
#   })
# 
# 
# ## ---
# ## In-depth: 8 models
# ## ---
# # filter for selected models
# log_decomposition_results_selected <- log_decomposition_results[selected_models]
# 
# ## important remark for BRIER score:
# # UNC = variance of X~Ber(p=E(y)=mean(y)) -> p(1-p)
# #p <- models_scoring_rules$boot_240$actual_conflict
# #mean(p)*(1-mean(p))
# 
# # combine to one dataframe
# log_score_decomposition_selected <- bind_rows(log_decomposition_results_selected) %>%
#   select(model, mean_score, MCB, DSC, UNC) %>%
#   rename(MeanScore = mean_score)
# 
# log_score_decomposition_barplot_selected <- log_score_decomposition_selected %>%
#   mutate(score_invisible = MeanScore, gap = 0, gap1=0, gap2 = 0)
# 
# # data frame with format for the barchart
# log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
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
# log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
#   mutate(model = factor(model, levels = unique(model[component == "MeanScore"])))
# 
# 
# # dataset for model orderbased on MeanScore
# log_levs_selected <- log_score_decomposition_barplot_selected %>%
#   filter(component == "MeanScore") %>%
#   arrange(desc(value)) %>%
#   pull(model)
# 
# # data for the plot
# log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
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
#     model = factor(model, levels = log_levs_selected, ordered = TRUE)
#   )
# 
# log_score_decomposition_barplot_selected$model <- recode(
#   log_score_decomposition_barplot_selected$model,
#   !!!model_labels
# )
# 
# 
# # set order of subbars within group
# log_score_decomposition_barplot_selected$component <- factor(
#   log_score_decomposition_barplot_selected$component,
#   levels = c("gap", "gap1", "gap2", "MCB", "UNC", "DSC", "score_invisible", "MeanScore")
# )
# 
# # set order of sub bars within group
# log_score_decomposition_barplot_selected$class_group <- factor(
#   log_score_decomposition_barplot_selected$class_group,
#   levels = c(4,1,6,2,3,5)    
# )
# 
# # set width of mean score bar
# mean_score_bar_selected <- 4.5
# # set width of mcb,dcs,unc bars
# msc_dsc_unc_bar_selected <- 2
# 
# # width column
# log_score_decomposition_barplot_selected <- log_score_decomposition_barplot_selected %>%
#   mutate(
#     wdth = case_when(
#       # MeanScore
#       class_group == 3 ~ mean_score_bar_selected, #4
#       # MCB DSC
#       class_group %in% c(1, 2) ~ msc_dsc_unc_bar_selected, #1
#       # gapswithin groups
#       class_group == 6 ~ msc_dsc_unc_bar_selected,
#       # gap between groups top
#       class_group == 5 ~ 15,
#       #gap between groups
#       class_group == 4 ~ 8 #10
#       
#     )
#   )
# 
# 
# log_score_decomposition_selected <- ggplot(log_score_decomposition_barplot_selected, aes(x = class_group, y = value, fill = component)) +
#   geom_bar(stat = "identity", position = "stack", width = log_score_decomposition_barplot_selected$wdth) +
#   facet_grid(model ~ ., switch = "y") +
#   coord_flip() +
#   scale_fill_manual(
#     values = c("gap1" = "darkgreen", "gap2" = "darkgreen","gap"="darkgreen","score_invisible" = "white","DSC" = "#808080", "UNC" = "#D9D9D9", "MCB" = "#A6A6A6", "MeanScore" = "#C05152"),          #"#1a1a2e"),  ##C05152 für MeanScore
#     ##
#     breaks = c("UNC", "MCB", "DSC", "MeanScore"), 
#     ##
#     name = ""
#   ) +
#   labs(title = "MSC-DSC Mean log Score Decomposition for No Previous Cases") + 
#   theme_classic() +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
#   theme(
#     panel.spacing = unit(0, "points"),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     strip.text.y.left = element_text(angle = 0, hjust = 1, size = 18),
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
# log_score_decomposition_selected
# 
# 
# if(store_plot == TRUE){
#   ggsave(paste0(folder,"/","log_score_decomposition_selected_models_no_prev_cases.png"),
#          plot = log_score_decomposition_selected, width = 1.0 * 4222, height = 1.0 * 2313, dpi = 300, units = "px",
#          bg="white")
# }
# 
# 



##------------------------------------------------------------------------------
## Plot: Ranking of the selected models for different situations
##------------------------------------------------------------------------------
crps_rankings_selected_models <- crps_selected_models %>%
  group_by(Model) %>%
  summarise(total_CRPS = sum(CRPS), .groups = "drop") %>%
  right_join(crps_selected_models, by = "Model") %>%
  select(Model, total_CRPS) %>%
  distinct() %>%
  arrange(total_CRPS) %>%
  mutate(rank_CRPS = row_number()) %>% 
  select(Model, rank_CRPS)

twcrps_rankings_selected_models <- twcrps_selected_models %>%
  group_by(Model) %>%
  summarise(total_twCRPS = sum(twCRPS), .groups = "drop") %>%
  right_join(twcrps_selected_models, by = "Model") %>%
  select(Model, total_twCRPS) %>%
  distinct() %>%
  arrange(total_twCRPS) %>%
  mutate(rank_twCRPS = row_number()) %>% 
  select(Model, rank_twCRPS)

brier_no_prev_cases_rankings_selected_models <- brier_selected_models %>%
  #filter(Situation == "onset" | Situation == "peace") %>%
  group_by(Model) %>%
  summarise(total_Brier = sum(BRIER), .groups = "drop") %>%
  right_join(brier_selected_models, by = "Model") %>%
  select(Model, total_Brier) %>%
  distinct() %>%
  arrange(total_Brier) %>%
  mutate(rank_Brier = row_number()) %>% 
  select(Model, rank_Brier)


model_rankings <- crps_rankings_selected_models %>%
  left_join(twcrps_rankings_selected_models, by = "Model") %>%
  left_join(brier_no_prev_cases_rankings_selected_models, by = "Model")

model_rankings <- model_rankings %>%
  rename(CRPS = rank_CRPS,
         twCRPS = rank_twCRPS,
         BRIER_Onset = rank_Brier)



# -----------------------------
# plateau-smoothstep function
# -----------------------------
smoothstep <- function(u) 3*u^2 - 2*u^3

plateau_step <- function(t, flat = 0.42) {
  flat <- pmin(pmax(flat, 0), 0.49)
  ifelse(
    t <= flat, 0,
    ifelse(t >= 1 - flat, 1, smoothstep((t - flat) / (1 - 2 * flat)))
  )
}

make_segment <- function(y0, y1, n = 160, flat = 0.42) {
  x <- seq(0, 1, length.out = n)
  s <- plateau_step(x, flat = flat)
  tibble(x = x, y = y0 + (y1 - y0) * s)
}

flat <- 0.33

left_curves <- model_rankings %>%
  rowwise() %>%
  do({
    seg <- make_segment(.$CRPS, .$twCRPS, flat = flat)
    seg$Model <- .$Model
    seg
  }) %>% ungroup()

right_curves <- model_rankings %>%
  rowwise() %>%
  do({
    seg <- make_segment(.$twCRPS, .$BRIER_Onset, flat = flat)
    seg$Model <- .$Model
    seg
  }) %>% ungroup()


left_pts  <- model_rankings %>% transmute(Model, x = 0, y = CRPS, x2 = 1, y2 = twCRPS)
right_pts <- model_rankings %>% transmute(Model, x = 0, y = twCRPS, x2 = 1, y2 = BRIER_Onset)

base_theme <- theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    strip.background = element_blank(),
    
    axis.line = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.y.left = element_blank(),
    axis.line.y.right = element_blank()
  )


ticks <- data.frame(y = 1:8)
tick_len <- 0.03
num_off  <- 0.055

axis_col <- "grey35"
axis_lwd <- 0.5

selected_colors_Model <- selected_colors
names(selected_colors_Model) <- model_labels[names(selected_colors)]

# ---------------- LEFT: CRPS (x=0) + twCRPS (x=1)
p_left <- ggplot() +
  geom_vline(xintercept = 0, colour = axis_col, linewidth = axis_lwd) +
  geom_vline(xintercept = 1, colour = axis_col, linewidth = axis_lwd) +
  
  geom_segment(
    data = ticks,
    aes(x = 0, xend = 0 - tick_len, y = y, yend = y),
    colour = axis_col, linewidth = 0.4
  ) +
  geom_text(
    data = ticks,
    aes(x = 0 - num_off, y = y, label = y),
    hjust = 1, colour = axis_col, size = 4
  ) +
  
  geom_segment(
    data = ticks,
    aes(x = 1, xend = 1 - tick_len, y = y, yend = y),
    colour = axis_col, linewidth = 0.35
  ) +
  
  geom_line(data = left_curves, aes(x = x, y = y, color = Model), linewidth = 1.1) +
  geom_point(data = left_pts,  aes(x = 0, y = y,  color = Model), size = 2) +
  geom_point(data = left_pts,  aes(x = 1, y = y2, color = Model), size = 2) +
  geom_text(data = model_rankings, aes(x = -0.15, y = CRPS, label = Model, color = Model),
            hjust = 1, size = 5) +
  annotate("text", x = 0 - 0.012, y = 0.6, label = "CRPS",
           fontface = "bold", hjust = 1, colour = axis_col, size = 5) +
  annotate("text", x = 1 - 0.012, y = 0.6, label = "twCRPS",
           fontface = "bold", hjust = 1, colour = axis_col, size = 5) +
  scale_color_manual(values = selected_colors_Model, guide = "none") +
  scale_x_continuous(limits = c(-0.15, 1), breaks = NULL, expand = c(0, 0)) +
  scale_y_reverse(limits = c(8.4, 0.6), breaks = NULL) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(plot.margin = margin(10, 0, 10, 160))

# ---------------- RIGHT: twCRPS (x=0) + BRIER (x=1)
p_right <- ggplot() +
  geom_vline(xintercept = 0, colour = axis_col, linewidth = axis_lwd) +
  geom_vline(xintercept = 1, colour = axis_col, linewidth = axis_lwd) +
  
  geom_segment(
    data = ticks,
    aes(x = 1, xend = 1 + tick_len, y = y, yend = y),
    colour = axis_col, linewidth = 0.4
  ) +
  geom_text(
    data = ticks,
    aes(x = 1 + num_off, y = y, label = y),
    hjust = 0, colour = axis_col, size = 4
  ) +
  
  geom_line(data = right_curves, aes(x = x, y = y, color = Model), linewidth = 1.1) +
  geom_point(data = right_pts,  aes(x = 0, y = y,  color = Model), size = 2) +
  geom_point(data = right_pts,  aes(x = 1, y = y2, color = Model), size = 2) +
  annotate("text", x = 1 - 0.012, y = 0.6, label = "Brier Score",
           fontface = "bold", hjust = 1, colour = axis_col, size = 5) +
  scale_color_manual(values = selected_colors_Model, guide = "none") +
  scale_x_continuous(limits = c(0, 1.1), breaks = NULL, expand = c(0, 0)) +
  scale_y_reverse(limits = c(8.4, 0.6), breaks = NULL) +
  base_theme +
  theme(plot.margin = margin(10, 80, 10, -90))

selected_models_ranking_plot <- (p_left | p_right) + plot_layout(widths = c(1.05, 1))

plot(selected_models_ranking_plot)


if(store_plot == TRUE){
  ggsave(paste0(folder,"/","selected_models_ranking_plot.png"),
         plot = selected_models_ranking_plot, width = 1.0 * 4222, height = 1.0 * 2000, dpi = 300, units = "px",
         bg="white")
}














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












