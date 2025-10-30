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
library(lubridate)

#devtools::install_github("aijordan/reliabilitydiag")

## ---------
# load ucdp data
## ---------

# Erkenne das Betriebssystem
os <- Sys.info()["sysname"]

# Setze den Pfad abhängig vom Betriebssystem
data_path <- ifelse(os == "Windows",
                    "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS",  # Windows
                    "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS")  # macOS/Linux

# ucdp ged dataset
ucdp_ged_event_data <- read.csv(file.path(data_path, "ged251-csv/GEDEvent_v25_1.csv"), header = TRUE)

ucdp_ged_event_data <- ucdp_ged_event_data %>%
  select(id, year, conflict_name, source_article, source_date,
         country, country_id, date_start, date_end, best)


ucdp_ged_event_data <- ucdp_ged_event_data %>%
  mutate(
    event_duration = as.Date(date_end) - as.Date(date_start),
    event_years = year(as.Date(date_end)) - year(as.Date(date_start)) + 1,
    event_months = month(as.Date(date_end)) - month(as.Date(date_start)) + 1
  ) %>%
  arrange(desc(event_duration))



ucdp_ged_multiple_month_event_data <- ucdp_ged_event_data %>%
  filter(event_months > 1,
         as.Date(date_start) >= as.Date("2018-01-01 00:00:00.000"),
         as.Date(date_end) <= as.Date("2023-12-31 00:00:00.000"))
  

#### greater 1 month i.e. those events that last more than one month
test_set_december_event_great1month_data <- ucdp_ged_event_data %>%
  filter(event_months > 1,
         as.Date(date_start) >= as.Date("2018-01-01 00:00:00.000"),
         as.Date(date_end) <= as.Date("2023-12-31 00:00:00.000"))


test_set_december_event_great1month_data <- test_set_december_event_great1month_data %>%
  filter(month(as.Date(date_end)) == 12)



test_set_december_event_great1month_data <- test_set_december_event_great1month_data %>%
  mutate(year = year(as.Date(date_end)),
         month = month(as.Date(date_end))) %>%
    group_by(country_id, year, month) %>%
    summarise(outcome_great1month = sum(best, na.rm = TRUE))



#### single month i.e. those events that only last 1 month

ucdp_ged_single_month_event_data <- ucdp_ged_event_data %>%
  filter(event_months == 1,
         as.Date(date_start) >= as.Date("2018-01-01 00:00:00.000"),
         as.Date(date_end) <= as.Date("2023-12-31 00:00:00.000"))

test_set_december_event_singlemonth_data <- ucdp_ged_single_month_event_data %>%
  filter(month(as.Date(date_end)) == 12)



test_set_december_event_singlemonth_data <- test_set_december_event_singlemonth_data %>%
  mutate(year = year(as.Date(date_end)),
         month = month(as.Date(date_end))) %>%
  group_by(country_id, year, month) %>%
  summarise(outcome_1month = sum(best, na.rm = TRUE))


### combine datasets
test_set_december_event_great1month_data_common <- test_set_december_event_great1month_data %>%
  inner_join(test_set_december_event_singlemonth_data, by = c("country", "year", "month")) %>%
  mutate(
    flag_true = (outcome_1month + outcome_great1month > 24) & (outcome_1month < 25),
    sum_outcomes = outcome_great1month + outcome_1month
  )
  
test_set_december_event_great1month_data_common <- test_set_december_event_great1month_data_common %>%
  filter(flag_true == TRUE)



# -----
## Do VIEWS, UCDP_GED, UCDP_GED_CANDIDATE coinside?
## Month: Dec 2023
# -----


## VIEWS data from 2018-2023
files_actuals_from18 <- list.files(data_path, pattern = "cm_actuals_\\d{4}\\.parquet", full.names = TRUE)
# read all files and bind them into a single data frame
observations_18_23 <- do.call(rbind, lapply(files_actuals_from18, arrow::read_parquet)) %>% # observations from 2018 - 2023
  arrange(country_id, month_id)

VIEWS_dec_23_data <- observations_18_23 %>%
  filter(month_id == 528)


## WICHTIG
# hier die country_ids anhand der namensliste von VIEWS mit denen von VIEWS ersetzen

## ucdp ged candicate dataset
ucdp_ged_event_candidate_dec_23_data <- read.csv(file.path(data_path, "ged251-csv/GEDEvent_v23_0_12.csv"), header = TRUE)

ucdp_ged_event_candidate_dec_23_data <- ucdp_ged_event_candidate_dec_23_data %>%
  filter(month(as.Date(date_end)) == 12) %>%
  group_by(country_id, year) %>%
  summarise(outcome = sum(best, na.rm = TRUE))

unique(ucdp_ged_event_candidate_dec_23_data$country)




ucdp_ged_event_dec_23_data <- ucdp_ged_event_data %>%
  filter(month(as.Date(date_end)) == 12,
         year(as.Date(date_end)) == 2023) %>%
  group_by(country_id) %>%
  summarise(outcome = sum(best, na.rm = TRUE))



