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

# actual data from 2018-2023 in the directory
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
         as.Date(date_end) <= as.Date("2023-12-31 00:00:00.000"),
         best > 25)
  


#length(unique(ucdp_ged_multiple_month_event_data$country))
  