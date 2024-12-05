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



## iterate over all rows of the data_2018_to_2023 data frame
for (i in 1:nrow(data_2018_to_2023)) {
  # get the country_id and month_id of the current row
  country_id <- data_2018_to_2023$country_id[i]
  month_id <- data_2018_to_2023$month_id[i]
  
  # get the fatality count of the current row
  fatality_count <- data_2018_to_2023$outcome[i]
  
  # initialize variables
  outbreak_level <- 0
  fatality_count_previous_month <- NA
  
  ## column "outbreak_level":
  # get the fatality count of the previous month
  if (month_id == 457){
    fatality_count_previous_month <- data_from_1990$ged_sb[which(data_from_1990$month_id == month_id - 1 & 
                                                                   data_from_1990$country_id == country_id)]
  } else{
    fatality_count_previous_month <- data_2018_to_2023$outcome[which(data_2018_to_2023$month_id == month_id - 1 & 
                                                                      data_2018_to_2023$country_id == country_id)]
  }
  # outbreak level
  if(fatality_count_previous_month == 0 & fatality_count > 0){
    outbreak_level <- fatality_count - fatality_count_previous_month
  }
  


  
  if(outbreak_level > 0){
    
    ## column "peace_months_prior":
    fatalities_month_prior <- fatality_count_previous_month
    prior_outbreak_month_id <- month_id
    peace_months_prior_count <- 0
    
    
    while(fatalities_month_prior == 0){
        peace_months_prior_count <- peace_months_prior_count + 1
        
        
        prior_outbreak_month_id <- prior_outbreak_month_id - 1
        
        if (prior_outbreak_month_id <= 457){
          fatalities_month_prior <- data_from_1990$ged_sb[which(data_from_1990$month_id == prior_outbreak_month_id & 
                                                                         data_from_1990$country_id == country_id)]
        } else{
          fatalities_month_prior <- data_2018_to_2023$outcome[which(data_2018_to_2023$month_id == prior_outbreak_month_id & 
                                                                             data_2018_to_2023$country_id == country_id)]
        }
        
        if(length(fatalities_month_prior) == 0){
          break
        }
    }
    
    ## column "conflict_months_after
    past_outbreak_month_id <- month_id + 1
    
    if(month_id == 528){
      fatalities_month_past = 0
    } else{
      fatalities_month_past <- data_2018_to_2023$outcome[which(data_2018_to_2023$month_id == past_outbreak_month_id & 
                                                                 data_2018_to_2023$country_id == country_id)]
    }
    
    
    conflict_months_past_count <- 0
    
    while(fatalities_month_past > 0){
      conflict_months_past_count <- conflict_months_past_count + 1
      
      past_outbreak_month_id <- past_outbreak_month_id + 1
      
      fatalities_month_past <- data_2018_to_2023$outcome[which(data_2018_to_2023$month_id == past_outbreak_month_id & 
                                                                    data_2018_to_2023$country_id == country_id)]
      
      if(length(fatalities_month_past) == 0){
        break
      }
      
    }
    
    ## add the data to the outbreak_data data frame
    outbreak_data <- rbind(outbreak_data, data.frame(outbreak_level = outbreak_level, 
                                                     peace_months_prior = peace_months_prior_count, 
                                                     conflict_months_after = conflict_months_past_count, 
                                                     month_id = month_id,
                                                     country_id = country_id))
    
  }  
}

## histogram conflict outbreaks
ggplot(outbreak_data, aes(x = outbreak_level)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.9) +
  scale_x_log10() +
  labs(
    title = "Intensity-Based Histogram of Outbreak Levels",
    x = "Outbreak Level",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1)
  )
  



## dataset and plot: number doutbreaks per country vs. avrg. outbreak magnitude
# number of outbreaks per country
data_outbreak_byCtry <- outbreak_data %>%
  group_by(country_id) %>%
  summarise(number_of_outbreaks = n(),
            avg_outbreak_level = mean(outbreak_level, na.rm = TRUE))

data_outbreak_byCtry_avg_magnitude <- data_outbreak_byCtry %>%
  group_by(number_of_outbreaks) %>%
  summarise(
    avg_outbreak_level = mean(avg_outbreak_level, na.rm = TRUE) 
  )

# plot
ggplot(data_outbreak_byCtry_avg_magnitude, aes(x = number_of_outbreaks, y = avg_outbreak_level)) + 
  geom_line(color = "steelblue", size = 1) +  
  geom_point(color = "navy", size = 3) +  
  labs(
    title = "Average Magnitude of Conflict-Outbreaks",
    x = "Number of Outbreaks per Country",
    y = "Average Magnitude of Outbreak"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1)
  )



## dataset and plot: time (in months) vs. number of outbreaks
# number of outbreaks per month
data_outbreak_byMonth <- outbreak_data %>%
  group_by(month_id) %>%
  summarise(number_of_outbreak = n(),
            number_of_outbreak_magn_greq_10 = sum(outbreak_level >= 10),
            number_of_outbreak_magn_greq_100 = sum(outbreak_level >= 100))

# row for month 523 (no outbreaks in this month)
new_row <- data.frame(
  month_id = 523,
  number_of_outbreak = 0,
  number_of_outbreak_magn_greq_10 = 0,
  number_of_outbreak_magn_greq_100 = 0
)

# add row
data_outbreak_byMonth <- rbind(
  data_outbreak_byMonth[1:66, ],
  new_row,
  data_outbreak_byMonth[67:nrow(data_outbreak_byMonth), ]
)

# Modify dataset for grouping
data_outbreak_byMonth_long <- data_outbreak_byMonth %>%
  pivot_longer(
    cols = c(
      number_of_outbreak,
      number_of_outbreak_magn_greq_10,
      number_of_outbreak_magn_greq_100
    ),
    names_to = "outbreak_type",
    values_to = "count"
  )

ggplot(data_outbreak_byMonth_long, aes(x = month_id, y = count, group = outbreak_type)) +
  geom_line(aes(color = outbreak_type), size = 0.8) +
  geom_area(aes(fill = outbreak_type), position = "identity", alpha = 0.4, show.legend = FALSE) +
  scale_color_manual(
    values = c(
      "number_of_outbreak" = "black",
      "number_of_outbreak_magn_greq_10" = "blue",
      "number_of_outbreak_magn_greq_100" = "red"
    ),
    labels = c(">= 0", ">= 10", ">= 100")
  ) +
  scale_fill_manual(
    values = c(
      "number_of_outbreak" = "lightgrey",
      "number_of_outbreak_magn_greq_10" = "lightblue",
      "number_of_outbreak_magn_greq_100" = "salmon"
    )
  ) +
  labs(
    title = "Monthly Outbreak Magnitude for all Countries",
    x = "Month ID",
    y = "Number of Outbreaks",
    color = "Outbreak\nMagnitude"
  ) +
  scale_x_continuous(
    breaks = c(457, 469, 481, 493, 505, 517, 528),
    labels = c("01.2018", "01.2019", "01.2020", "01.2021", "01.2022", "01.2023", "12.2023")
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.1)
  )
