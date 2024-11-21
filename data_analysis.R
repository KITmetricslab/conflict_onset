# descriptive data analysis of the actual years 2018-2023

# load packages
library(arrow)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

# path to directory
data_path <- "C:/Users/fn3307/Documents/data views/"

# find all parquet files in the directory
files <- list.files(data_path, pattern = "cm_actuals_\\d{4}\\.parquet", full.names = TRUE)

# read all files and bind them into a single data frame
data_2018_to_2023 <- do.call(rbind, lapply(files, arrow::read_parquet))

# read the feature dataset from 1990
data_from_1990 <- arrow::read_parquet("C:/Users/fn3307/Documents/data views/cm_features.parquet")
# drop all columns but "country_id", "month_id" and "ged_sb"
data_from_1990 <- data_from_1990 %>% select(country_id, month_id, ged_sb)


# create new dataset that tracks each outbreak in the period 2018-2023
# it has the columns "outbreak_level", "peace_months_prior", "conflict_months_after", "month_id", "country_id"

### do for every row (that is for every country and every month) so at time t(country)
## column "outbreak_level":
# if fatalities(t-1) = 0  and fatalities (t) > 0 then outbreak_level(t) = fat(t) - fat(t-1) = fat(t)
## column "peace_months_prior":
# if outbreak_level(t) > 0 then 
# peace_months_prior(t) = |{z| z in [z_start, t), fat(z) = 0, z_start = min{x in [z_start, t)| f(x)=0 and f(x - 1) > 0}}| 
## column "conflict_months_after":
# if outbreak_level(t) > 0 then
# conflict_months_after(t) = |{z| z in (t, z_end], fat(z) >0, z_end = min{x in (t, z_end]| f(x)>0 and f(x + 1) = 0}|

# create a new data frame that will store the outbreak data
outbreak_data <- data.frame(outbreak_level = integer(), 
                            peace_months_prior = integer(), 
                            conflict_months_after = integer(), 
                            month_id = integer(), 
                            country_id = integer())


### Achtung vorläufig und nicht fertig!

# iterate over all rows of the data_2018_to_2023 data frame
for (i in 1:nrow(data_2018_to_2023)) {
  # get the country_id and month_id of the current row
  country_id <- data_2018_to_2023$country_id[i]
  month_id <- data_2018_to_2023$month_id[i]
  
  # get the fatality count of the current row
  fatality_count <- data_2018_to_2023$outcome[i]
  
  # get the fatality count of the previous month
  if month_id = 457:
    #fatality_count_previous_month <- # data_from_1990$ged_sb[which(data_from_1990$country_id == country usw.
  else:
    fatality_count_previous_month <- data_2018_to_2023$ged_sb[i - 1]
  
  # outbreak level
  outbreak_level <- fatality_count - fatality_count_previous_month

  
  
  
  # peace months prior
  
  
  
  # conflict months after
  
  
  if outbreak_level > 0:
    # add the data to the outbreak_data data frame
    outbreak_data <- rbind(outbreak_data, data.frame(outbreak_level = outbreak_level, 
                                                     peace_months_prior = peace_months_prior, 
                                                     conflict_months_after = conflict_months_after, 
                                                     month_id = month_id)
}