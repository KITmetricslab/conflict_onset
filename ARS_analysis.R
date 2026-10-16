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


##----------------------------
# dataset of reference year
##----------------------------
df_actual_2018 <- df_actual %>%
  filter(Year == "2018") %>%
  select(-Activity, 
         -Total.human.disease.cases,
         -X..Presumptive.viremic.blood.donors, 
         -Notes)

df_actual_2018 <- df_actual_2018 %>%
  rename(FIPS = Location,
         actual = Neuroinvasive.disease.cases) %>%
  select(-FullGeoName)

df_actual_2018 <- left_join(df_actual_2018, geocodes, by = "FIPS")

new_actual_rows_for_counties_no_data_available_2018 <- geocodes %>%
  anti_join(df_actual_2018, by = "FIPS") %>%
  
  mutate(
    Year = 2018,
    actual = 0
  ) %>%
  
  select(Year, FIPS, actual, location)

df_actual_2018 <- bind_rows(df_actual_2018, new_actual_rows_for_counties_no_data_available_2018)

# sort by FIPS
df_actual_2018 <- df_actual_2018 %>%
  arrange(FIPS)

# 2018 für Join vorbereiten
actual_2018_join <- df_actual_2018 %>%
  select(FIPS, actual_2018 = actual)


## -----
## create list of predictive probabilites for all models
## -----

## join actual and predictions dataset
df_predictions <- left_join(df_predictions, df_actual_2020, by = "location") %>%
  select(-target)%>%
  relocate(Year, FIPS, .after = 1)


# 2018 an Prediction-Datensatz hängen
df_predictions <- df_predictions %>%
  left_join(actual_2018_join, by = "FIPS")

















##----------------------------------------------------------------------------------
# Analysis
##----------------------------------------------------------------------------------
ars_raw_bins <- df_predictions %>%
  filter(team == "ARS") %>%
  group_by(location, FIPS, Year, actual, actual_2018) %>%
  mutate(bin = row_number()) %>%
  ungroup()

##-----
# 2018 = 0
##-----
check_zero_2018 <- ars_raw_bins %>%
  filter(actual_2018 == 0) %>%
  group_by(location, FIPS) %>%
  summarise(
    first_bin = value[bin == 1],
    other_bins_sum = sum(value[bin > 1]),
    .groups = "drop"
  ) %>%
  mutate(
    over_99_percent_zero = first_bin > 0.99
  )

check_zero_2018 %>%
  summarise(
    n_counties = n(),
    n_over_99 = sum(over_99_percent_zero),
    n_not_over_99 = sum(!over_99_percent_zero),
    proportion_over_99 = mean(over_99_percent_zero)
  )

# von insgesamt 2496 mit = 0 in 2018 haben 2494 counties > 99% Wk Masse im ersten bin

##-----
# 2018 > 0
##-----
check_nonzero_2018 <- ars_raw_bins %>%
  filter(actual_2018 > 0) %>%
  group_by(location, FIPS) %>%
  summarise(
    first_bin = value[bin == 1],
    other_bins_sum = sum(value[bin > 1]),
    .groups = "drop"
  ) %>%
  mutate(
    over_99_percent_zero = first_bin > 0.99
  )

check_nonzero_2018 %>%
  summarise(
    n_counties = n(),
    n_over_99 = sum(over_99_percent_zero),
    n_not_over_99 = sum(!over_99_percent_zero),
    proportion_over_99 = mean(over_99_percent_zero)
  )

# von insgesamt 612 mit > 0 in 2018 haben nur 7 counties > 99% Wk Masse im ersten bin


# korrelation anzahl WNND in 2018 mit skewness?
ars_simple <- ars_raw_bins %>%
  filter(actual_2018 > 0) %>%
  group_by(location, FIPS, actual_2018) %>%
  summarise(
    p_zero = value[bin == 1],
    p_over_20 = sum(value[bin >= 6]),
    p_over_50 = sum(value[bin >= 12]),
    p_over_100 = sum(value[bin >= 13]),
    .groups = "drop"
  )

correlations <- tibble(
  measure = c("p_zero", "p_over_20", "p_over_50", "p_over_100"),
  
  spearman_rho = c(
    cor(ars_simple$actual_2018, ars_simple$p_zero,
        method = "spearman", use = "complete.obs"),
    
    cor(ars_simple$actual_2018, ars_simple$p_over_20,
        method = "spearman", use = "complete.obs"),
    
    cor(ars_simple$actual_2018, ars_simple$p_over_50,
        method = "spearman", use = "complete.obs"),
    
    cor(ars_simple$actual_2018, ars_simple$p_over_100,
        method = "spearman", use = "complete.obs")
  )
)

correlations




# Among counties with prior WNVND activity, the number of cases in 2018 was negatively correlated 
# with the probability assigned to zero cases (Spearman’s ρ = −0.67) and positively correlated with 
# the probability assigned to more than 20 (ρ = 0.65), 50 (ρ = 0.64), and 100 cases (ρ = 0.57).



