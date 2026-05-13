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
                 "NYSW-CVD", "Rutgers", "Standford", "UA", 
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
df_actual <- df_actual %>%
  filter(Year == "2020") %>%
  select(-Activity, 
         -Total.human.disease.cases,
         -X..Presumptive.viremic.blood.donors, 
         -Notes)

df_actual <- df_actual %>%
  rename(FIPS = Location,
         actual = Neuroinvasive.disease.cases) %>%
  select(-FullGeoName)

df_actual <- left_join(df_actual, geocodes, by = "FIPS")

new_actual_rows_for_counties_no_data_available <- geocodes %>%
  anti_join(df_actual, by = "FIPS") %>%
  
  mutate(
    Year = 2020,
    actual = 0
  ) %>%
  
  select(Year, FIPS, actual, location)

df_actual <- bind_rows(df_actual, new_actual_rows_for_counties_no_data_available)

# sort by FIPS
df_actual <- df_actual %>%
  arrange(FIPS)




## -----
## create list of predictive probabilites for all models
## -----

## join actual and predictions dataset
df_predictions <- left_join(df_predictions, df_actual, by = "location") %>%
  select(-target)%>%
  relocate(Year, FIPS, .after = 1)









## ------------------------------------------------------------------------------------------------------------
## define lower threshold for binary event of a present outbreak
## -----
lower_infections_thresh = 0










# upper bounds of the bins (quantiles)
qs <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 150, 200, 999)











































### example

example_forecast_ARS <- df_predictions$value[17:31]
#example_forecast_ARS <- df_predictions$value[1:15]
bin_probs <- example_forecast_ARS



# cdf (ps)
ps <- cumsum(bin_probs)
# WICHTIG: Falls die Summe durch Rundungsfehler minimal über 1 liegt, auf 1 begrenzen
ps[ps > 1] <- 1


x <- seq(from = -0.5, to = 1000, length = 1000)

p_lognormal_approx <- make_p_fn(ps = ps,
                                qs = qs,
                                tail_dist = "lnorm")
cdf_lognormal_approx <- p_lognormal_approx(x)

data.frame(
  x = x,
  y = cdf_lognormal_approx
) %>%
  ggplot() +
  geom_line(
    mapping = aes(x = x, y = y),
    size = 0.8
  ) +
  geom_point(
    data = data.frame(q = qs, p = ps),
    mapping = aes(x = q, y = p),
    size = 1.2
  ) +
  scale_x_log10() + 
  ylim(0, 1) +
  ylab("Probability") +
  xlab("") +
  theme_bw() + 
  xlim(0,20)




x_counts <- 0:999  # Alternativ: seq(from = 0, to = 1000, by = 1)

# 2. Werte die CDF an genau diesen ganzzahligen Stellen aus
cdf_at_counts <- p_lognormal_approx(x_counts)

# Wenn du für weitere Berechnungen die Wahrscheinlichkeiten für JEDEN EINZELNEN Count brauchst
# (also P(X=0), P(X=1), P(X=2) etc. aus der Treppenfunktion ableiten willst):
# Das ist die Differenz zwischen dem aktuellen und dem vorherigen CDF-Wert.
prob_at_counts <- c(cdf_at_counts[1], diff(cdf_at_counts))

# 3. Speichere das Ergebnis übersichtlich in einem Data Frame
df_counts <- data.frame(
  Count = x_counts,
  CDF = cdf_at_counts,       # P(X <= Count)
  Probability = prob_at_counts # P(X == Count)
)


actual <- 0


indicator <- as.numeric(x_counts >= actual)
weight <- 1 / (x_counts + 1)

# RPS discrete version of CRPS
rps <- sum((df_counts$CDF - indicator)^2)
# twRPS discrete version of twCRPS
twrps <- sum(weight * (df_counts$CDF - indicator)^2)

# BRIER Score for Onset meaning > 0 WNVND infections




## für alle predictions aus der CDF (die für alle auf die selbe Art und Weise erzeugt wird) die WKs
# für die count-valued outcomes ablesen( 0-999 ). Das ist implizit die Treppenfunktion (und somit CRPS = RPS).
# Wie genau CDF erzeugt wird laut Johannes nicht so wichtig.





