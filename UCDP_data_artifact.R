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
library(lubridate)

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

month_lookup <- tibble(
  month_id = 457:528,
  month_date = seq(
    as.Date("2018-01-01"),
    as.Date("2023-12-01"),
    by = "month"
  )
)

# Check whether all month IDs match
stopifnot(setequal(observations_18_23$month_id, month_lookup$month_id))



## all GED events to compare with VIEWS actuals
GED_events_all <- read_csv(
  file.path(data_path, "ged251-csv", "GEDEvent_v24_1.csv"), #GEDEvent_v25_1.csv
  show_col_types = FALSE,
  guess_max = Inf
) %>%
  filter(
    type_of_violence == 1,
    year %in% 2018:2023
  ) %>%
  rename(
    country_id_ged = country_id,
    ged_date_start = date_start,
    ged_date_end = date_end
  ) %>%
  mutate(
    country_match = country %>%
      str_remove("\\s*\\([^)]*\\)") %>%
      str_trim() %>%
      recode("DR Congo" = "Congo, DRC"),
    
    ged_date_start = as.Date(ged_date_start),
    ged_date_end = as.Date(ged_date_end),
    
    # Assign each event to its ending month
    month_date = lubridate::floor_date(ged_date_end, "month")
  ) %>%
  left_join(
    country_id_list %>% select(country_id, name),
    by = c("country_match" = "name")
  ) %>%
  relocate(country_match, country_id, .after = country_id_ged)

GED_events_all <- GED_events_all %>%
  left_join(
    month_lookup %>% rename(month_id_end = month_id),
    by = "month_date"
  ) %>%
  select(-month_date)

# check if all countries are assigned
GED_events_all %>%
  filter(is.na(country_id)) %>%
  count(country, sort = TRUE)


## aggregate GED on country-month level----------------------------------------

# assign all fatalities from events that end in a month for a given country
# to the outcome variable
GED_monthly <- GED_events_all %>%
  filter(!is.na(country_id)) %>%
  mutate(
    ged_event_months =
      (lubridate::year(ged_date_end) -
         lubridate::year(ged_date_start)) * 12 +
      lubridate::month(ged_date_end) -
      lubridate::month(ged_date_start) + 1
  ) %>%
  group_by(country_id, month_id_end) %>%
  summarise(
    ged_best = sum(best, na.rm = TRUE),
    n_events = n(),  # Total number of events assigned to the ending month
    
    n_multimonth = sum(ged_event_months > 1, na.rm = TRUE), # Number of multimonth events assigned to the ending month
    
    fatalities_multimonth = sum(
      best[ged_event_months > 1],
      na.rm = TRUE
    ), # Total fatalities from multimonth events assigned to the ending month
    
    .groups = "drop"
  )

# Check whether all GED country-month pairs are present in the VIEWS actuals
GED_monthly %>%
  anti_join(
    observations_18_23,
    by = c("country_id", "month_id_end" = "month_id")
  ) 


## compare GED with VIEWS actuals ----------------------------------------------
monthly_comparison <- observations_18_23 %>%
  left_join(
    GED_monthly,
    by = c("country_id", "month_id" = "month_id_end")
  ) %>%
  mutate(
    ged_best = coalesce(ged_best, 0),
    n_events = coalesce(n_events, 0L),
    n_multimonth = coalesce(n_multimonth, 0L),
    fatalities_multimonth = coalesce(fatalities_multimonth, 0),
    difference = outcome - ged_best,
    
    share_multimonth = if_else(
      outcome > 0,
      fatalities_multimonth / outcome * 100,
      NA_real_
    )
  )


## analysis of affected country-month-pairs------------------------------------
monthly_affected <- monthly_comparison %>%
  filter(n_multimonth > 0) %>%
  arrange(desc(fatalities_multimonth))









## individual multimonth events -----------------------------------------------
GED_multimonth_details <- GED_events_all %>%
  mutate(
    ged_event_months =
      (lubridate::year(ged_date_end) -
         lubridate::year(ged_date_start)) * 12 +
      lubridate::month(ged_date_end) -
      lubridate::month(ged_date_start) + 1
  ) %>%
  filter(
    ged_event_months > 1,
    !is.na(country_id)
  ) %>%
  select(
    id, country_id, month_id_end,
    ged_date_start, ged_date_end,
    ged_event_months, best
  ) %>%
  arrange(desc(best))


GED_multimonth_details %>%
  filter(month_id_end == 514,
         country_id == 57)
