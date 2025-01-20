

#### --------------------------------------------------------
# visualisation of onset probabilities of the years 2018-2023 
#### --------------------------------------------------------


# load packages
library(arrow)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(pbapply)


# path to directory on "share-alle"
data_path <- "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS Data/"

# actual data from 2018-2023 in the directory
files_actuals_from18 <- list.files(data_path, pattern = "cm_actuals_\\d{4}\\.parquet", full.names = TRUE)
# read all files and bind them into a single data frame
data_2018_to_2023 <- do.call(rbind, lapply(files_actuals_from18, arrow::read_parquet))

# feature data from 1990
path_features <- paste0(data_path,"cm_features.parquet")
data_from_1990 <- arrow::read_parquet(path_features)
# drop all columns but "country_id", "month_id" and "ged_sb"
data_from_1990 <- data_from_1990 %>% select(country_id, month_id, ged_sb)


## ongoing peace --------------------------------------------------------

# create dataset for all ongoing peace instances
# this dataset is used later on for the predicted onset probability distribution
ongoing_peace_data <- data.frame(fatalities = integer(),
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
  
  if(fatality_count_previous_month == 0 & fatality_count == 0){
  
    ## add the data to data frame
    ongoing_peace_data <- rbind(ongoing_peace_data, data.frame(fatalities = fatality_count,
                                                               month_id = month_id, 
                                                               country_id = country_id))
  }  
}



## conflict onset: peace_months_prior >= 1 --------------------------------------------------------
# create new dataset that tracks each outbreak in the period 2018-2023
# it has the columns "outbreak_level", "peace_months_prior", "conflict_months_after", "month_id", "country_id"

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

## histogram peace months prior
ggplot(outbreak_data, aes(x = peace_months_prior)) +
  geom_histogram(binwidth = 0.05, fill = "#90EE90", color = "black", alpha = 0.9) +
  scale_x_log10() +
  labs(
    title = "Months of Peace Prior to Conflict Onsets (2018-2023)",
    x = "Months in Peace Prior to Conflict Onset",
    y = "Frequency"
  ) +
  theme_bw()


## histogram conflict onset magnitude
ggplot(outbreak_data, aes(x = outbreak_level)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black", alpha = 0.9) +
  scale_x_log10() +
  labs(
    title = "Intensity of Conflict Onsets: Fatalities in Initial Month (2018-2023)",
    subtitle = "Direct Onset: 2018-2023",
    x = "Fatalities",
    y = "Frequency"
  ) +
  theme_bw()




#### --------------------------------------------------------
# distribution of onset probabilities
#### --------------------------------------------------------

## define functions for the different distribution plots
create_density_plot <- function(data, title, subtitle, fill_color = "steelblue") {
  ggplot(data, aes(x = probability_gr_0)) +
    geom_density(color = "black", fill = fill_color, size = 1, alpha = 0.6) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Onset Probability",
      y = "Density"
    ) +
    theme_bw()
}

create_boxplot <- function(data, title, subtitle, fill_color = "steelblue", width = 0.6, limits = NULL) {
  p <- ggplot(data, aes(y = probability_gr_0)) +
    geom_boxplot(fill = fill_color, alpha = 0.6, color = "black", size = 0.6, outlier.color = "red", outlier.size = 2, width = width) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Onset Probability"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  if (!is.null(limits)) {
    p <- p + scale_y_continuous(limits = limits)
  }
  return(p)
}

create_facet_density <- function(data, title, subtitle, fill_color = "steelblue") {
  ggplot(data, aes(x = probability_gr_0)) +
    geom_density(fill = fill_color, alpha = 0.6, size = 0.8) +
    facet_wrap(~ model, scales = "free_y") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Onset Probability",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )
}

#create_facet_boxplot









# path to directory
base_path <- "C:/Users/fn3307/Documents/data views/all_available_cm_predictions"








# subfolders for all models
model_dirs <- list.dirs(base_path, recursive = FALSE)

# returns one dataframe of all years from one model
process_model <- function(model_path) {
  # all directories of the years
  year_dirs <- list.dirs(file.path(model_path, "cm"), recursive = FALSE)
  
  # 
  greater_zero_data <- year_dirs %>%
    map(list.files, full.names = TRUE, pattern = "\\.parquet") %>%  # performs list.files on every dir. of year_dirs; returns the full filepath names
    flatten_chr() %>%                                              # converts list of filepaths into one vector
    map_dfr(arrow::read_parquet)                                   # performs read_parquet on every element in the vector
  
  return(greater_zero_data) # dataframe with all parquet files in year_dirs
}

# create list for all models
list_all_data <- model_dirs %>%
  set_names(basename(.)) %>% # set names of the list entries, basenames last name of the filepath in model_dirs
  map(process_model) # perform process_model function



# list with dataframes containing P(Y > 0)
list_prob_onset <- list()
list_model_names <- names(list_all_data)

# list of months (Jan. 18 until Dec. 23)
month_list <- seq(457,457 + 12*6 - 1, by = 1)
# country_list (equal for each model)
country_list <- unique(unlist(list_all_data$bodentien_rueter_negbin$country_id))


for (j in seq_along(list_all_data)) {
  
  model = list_all_data[[j]]
  
  prob_gr_zero_model_data <- model %>%
    filter(month_id %in% month_list, country_id %in% country_list) %>% # keep all rows that are in month_list and country_list
    group_by(month_id, country_id) %>%  # group by month and country (each combination of month_id and country_id is one group)
    summarise( # summarise creates new dataframe and calculates values for each group
      # returns one row for each combination of grouping variables
      prob_gr_0 = 1 - sum(outcome == 0) / n(),
      .groups = "drop"  #  result dataframe is not grouped (mutate is performed for the whole dataframe)
    ) %>%
    mutate(model = list_model_names[j]) # adds new column
  
  list_prob_onset[[list_model_names[j]]] <- prob_gr_zero_model_data  
  
  cat("\rFinished", j, "of", length(list_all_data))
}






###############
# density plots
###############


## ongoing peace --------------------------------------------------------

# join outbreak data and the prob_gr_0 observations from prob_data
for (model_name in names(list_prob_onset)) {
  
  prob_data_ongoing_peace <- list_prob_onset[[model_name]]
  
  prob_data_ongoing_peace <- prob_data_ongoing_peace %>%
    rename(!!paste0("prob_gr_0_", model_name) := prob_gr_0)
  
  prob_data_ongoing_peace <- prob_data_ongoing_peace %>%
    select(-model)
  
  ongoing_peace_data <- ongoing_peace_data %>%
    left_join(prob_data_ongoing_peace, by = c("month_id", "country_id"))
}

# delete column for zero model
ongoing_peace_data <- ongoing_peace_data %>%
  select(-prob_gr_0_zero)

prob_columns_ongoing_peace <- ongoing_peace_data %>%
  select(starts_with("prob_gr_0"))

# long format for the density plots
prob_long_ongoing_peace <- prob_columns_ongoing_peace %>%
  pivot_longer(
    cols = everything(),
    names_to = "model",
    values_to = "probability_gr_0"
  ) %>%
  mutate(model = sub("^prob_gr_0_", "", model))


# density plot for every observation over all models
ggplot(prob_long_ongoing_peace, aes(x = probability_gr_0)) +
  geom_density(color="black",fill="#90EE90", size = 1, alpha = 0.6) +
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: All Models",
    subtitle = "Ongoing Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_bw()


# density for each model
ggplot(prob_long_ongoing_peace, aes(x = probability_gr_0)) +
  geom_density(fill = "#90EE90", alpha = 0.6, size = 0.8) + 
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: Individual Models",
    subtitle = "Ongoing Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )






## conflict onset: peace_months_prior >= 1 --------------------------------------------------------

# join outbreak data and the prob_gr_0 observations from prob_data











# benötigt man das hier jedes mal???


for (model_name in names(list_prob_onset)) {
  
  prob_data <- list_prob_onset[[model_name]]
  
  prob_data <- prob_data %>%
    rename(!!paste0("prob_gr_0_", model_name) := prob_gr_0)
  
  prob_data <- prob_data %>%
    select(-model)
  
  outbreak_data <- outbreak_data %>%
    left_join(prob_data, by = c("month_id", "country_id"))
}

# delete column for zero model
outbreak_data <- outbreak_data %>%
  select(-prob_gr_0_zero)

prob_columns <- outbreak_data %>%
  select(starts_with("prob_gr_0"))













# long format for the density plots
prob_long <- prob_columns %>%
  pivot_longer(
    cols = everything(),
    names_to = "model",
    values_to = "probability_gr_0"
  ) %>%
  mutate(model = sub("^prob_gr_0_", "", model))


# density plot for every observation over all models
ggplot(prob_long, aes(x = probability_gr_0)) +
  geom_density(color="black",fill="steelblue", size = 1, alpha = 0.6) +
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: All Models",
    subtitle = paste0("Direct Onset: 2018-2023 (", length(outbreak_data$outbreak_level), " Onsets)"),
    x = "onset probability",
    y = "density"
  ) +
  theme_bw()

# boxplot for every observation over all models
ggplot(prob_long, aes(y = probability_gr_0)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, color = "black", size = 0.8, outlier.color = "red", outlier.size = 2) +
  labs(
    title = "Boxplot of Onset Probabilities 2018-2023: All Models",
    subtitle = paste0("Direct Onset: 2018-2023 (", length(outbreak_data$outbreak_level), " Onsets)"),
    x = NULL, 
    y = "Onset Probability"
  ) +
  theme_bw()

# density for each model
ggplot(prob_long, aes(x = probability_gr_0)) +
  geom_density(fill = "steelblue", alpha = 0.6, size = 0.8) + 
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: Individual Models",
    subtitle = paste0("Direct Onset: 2018-2023 (", length(outbreak_data$outbreak_level), " Onsets)"),
    x = "onset probability",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# boxplot for each model
ggplot(prob_long, aes(x = probability_gr_0)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6, color = "black", size = 0.6, outlier.color = "red", outlier.size = 1) +
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Boxplot of Onset Probabilities 2018-2023: Individual Models",
    subtitle = paste0("Direct Onset: 2018-2023 (", length(outbreak_data$outbreak_level), " Onsets)"),
    x = "Onset Probability",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )


## conflict onset: peace_months_prior >= 12 --------------------------------------------------------

# join outbreak data and the prob_gr_0 observations from prob_data
for (model_name in names(list_prob_onset)) {
  
  prob_data_12m_peace <- list_prob_onset[[model_name]]
  
  prob_data_12m_peace <- prob_data_12m_peace %>%
    rename(!!paste0("prob_gr_0_", model_name) := prob_gr_0)
  
  prob_data_12m_peace <- prob_data_12m_peace %>%
    select(-model)
  
  outbreak_data_12m_peace <- outbreak_data_12m_peace %>%
    left_join(prob_data_12m_peace, by = c("month_id", "country_id"))
}

# delete column for zero model
outbreak_data_12m_peace <- outbreak_data_12m_peace %>%
  select(-prob_gr_0_zero)

prob_columns_12m_peace <- outbreak_data_12m_peace %>%
  select(starts_with("prob_gr_0"))

# long format for the density plots
prob_long_12m_peace <- prob_columns_12m_peace %>%
  pivot_longer(
    cols = everything(),
    names_to = "model",
    values_to = "probability_gr_0"
  ) %>%
  mutate(model = sub("^prob_gr_0_", "", model))


# density plot for every observation over all models
ggplot(prob_long_12m_peace, aes(x = probability_gr_0)) +
  geom_density(color="black",fill="#2F4F4F", size = 1, alpha = 0.6) +
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: All Models",
    subtitle = paste0("One Year Prolonged Peace Onset: 2018-2023 (", length(outbreak_data_12m_peace$outbreak_level), " Onsets)"),
    x = "onset probability",
    y = "density"
  ) +
  theme_bw()




# boxplot for every observation over all models
ggplot(prob_long_12m_peace, aes(y = probability_gr_0)) +
  geom_boxplot(fill = "#2F4F4F", alpha = 0.6, color = "black", size = 0.8, outlier.color = "red", outlier.size = 2) +
  labs(
    title = "Boxplot of Onset Probabilities 2018-2023: All Models",
    subtitle = paste0("One Year Prolonged Peace Onset: 2018-2023 (", length(outbreak_data_12m_peace$outbreak_level), " Onsets)"),
    x = NULL, 
    y = "Onset Probability"
  ) +
  theme_bw()

# density for each model
ggplot(prob_long_12m_peace, aes(x = probability_gr_0)) +
  geom_density(fill = "#2F4F4F", alpha = 0.6, size = 0.8) + 
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: Individual Models",
    subtitle = paste0("One Year Prolonged Peace Onset: 2018-2023 (", length(outbreak_data_12m_peace$outbreak_level), " Onsets)"),
    x = "onset probability",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# boxplot for each model
ggplot(prob_long_12m_peace, aes(x = probability_gr_0)) +
  geom_boxplot(fill = "#2F4F4F", alpha = 0.6, color = "black", size = 0.6, outlier.color = "red", outlier.size = 1) +
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Boxplot of Onset Probabilities 2018-2023: Individual Models",
    subtitle = paste0("One Year Prolonged Peace Onset: 2018-2023 (", length(outbreak_data_12m_peace$outbreak_level), " Onsets)"),
    x = "Onset Probability",
    y = NULL
  ) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )



## previous month has peace --------------------------------------------------------------------------

prob_long_month_tMinus1_zero <- rbind(prob_long, prob_long_ongoing_peace)


# density plot for every observation over all models
ggplot(prob_long_month_tMinus1_zero, aes(x = probability_gr_0)) +
  geom_density(color="black",fill="#5F2F4F", size = 1, alpha = 0.6) +
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: All Models",
    subtitle = "Previous Month Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_bw()

# boxplot for every observation over all models
ggplot(prob_long_month_tMinus1_zero, aes(y = probability_gr_0)) +
  geom_boxplot(fill = "#5F2F4F", alpha = 0.6, color = "black", size = 0.6, outlier.color = "red", outlier.size = 2, width = 0.4) +
  labs(
    title = "Boxplot of Onset Probabilities 2018-2023: All Models",
    subtitle = "Previous Month Peace: 2018-2023",
    x = NULL, 
    y = "Onset Probability"
  ) +
  scale_x_continuous(
    limits = c(-0.4, 0.4)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()
  )

# density for each model
ggplot(prob_long_month_tMinus1_zero, aes(x = probability_gr_0)) +
  geom_density(fill = "#5F2F4F", alpha = 0.6, size = 0.8) + 
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Distribution of Onset Probabilities 2018-2023: Individual Models",
    subtitle = "Previous Month Peace: 2018-2023",
    x = "onset probability",
    y = "density"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )

# boxplot for each model
ggplot(prob_long_month_tMinus1_zero, aes(x = probability_gr_0)) +
  geom_boxplot(fill = "#5F2F4F", alpha = 0.6, color = "black", size = 0.6, outlier.color = "red", outlier.size = 1) +
  facet_wrap(~ model, scales = "free_y") +  # facet for each model
  labs(
    title = "Boxplot of Onset Probabilities 2018-2023: Individual Models",
    subtitle = "Previous Month Peace: 2018-2023",
    x = "Onset Probability",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold", size = 14)
  )
