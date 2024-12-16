# distribution of onset probabilities

# load packages
library(arrow)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(pbapply)

# path to directory
base_path <- "C:/Users/fn3307/Documents/data views/all_available_cm_predictions"

# subfolders for all models
model_dirs <- list.dirs(base_path, recursive = FALSE)

list_all_data <- list()

for (model in model_dirs) {
  year_dirs <- list.dirs(file.path(model, "cm"), recursive = FALSE)
  outbreak_data <- data.frame()
  
  print(basename(model))
  
  for (year in year_dirs) {
    
    files <- list.files(year, full.names = TRUE, pattern = "\\.parquet")
    
    # read all files and bind them into a single data frame
    model_data <- do.call(rbind, lapply(files, arrow::read_parquet))
    
    outbreak_data <- rbind(outbreak_data, model_data)
  }
  
  list_all_data[[basename(model)]] <- outbreak_data
}


# list with dataframes containing P(Y > 0)
list_prob_onset <- list()
list_model_names <- names(list_all_data)

i <- 1
for (model in list_all_data) {
  month_list <- unique(unlist(model$month_id))
  month_list <- head(month_list, -12)
  
  prob_gr_zero_model_data <- data.frame(month_id = numeric(), country_id = numeric(), 
                                        prob_gr_0 = numeric(), model = character())
  
  prob_gr_zero_model_data <- model %>%
    filter(month_id %in% month_list, country_id %in% country_list) %>% # keep all rows that are in month_list and country_list
    group_by(month_id, country_id) %>%  # group by month and country (each combination of month_id and country_id is one group)
    summarise( # summarise creates new dataframe and calculates values for each group
      # returns one row for each combination of grouping variables
      prob_gr_0 = 1 - sum(outcome == 0) / n(),
      .groups = "drop"  #  result dataframe is not grouped (mutate is performed for the whole dataframe)
    ) %>%
    mutate(model = list_model_names[i]) # adds new column
  
  list_prob_onset[[list_model_names[i]]] <- prob_gr_zero_model_data      
  
  i <- i + 1
  
  cat("\rFinished", i-1, "of", length(list_all_data))
}



data_bodentien_rueter <- list_prob_onset$bodentien_rueter_negbin


# month country list combination for onset and ongoing peace

# group each list df's by those combinations and plot dist. per model








































# two lists with dataframes containing P(Y > 0) for all
# a) onset-instances and b) ongoing peace instances 

data_bodentien_rueter <- list_all_data$bodentien_rueter_negbin
month_list <- unique(unlist(list_all_data$bodentien_rueter_negbin$month_id))
month_list <- head(month_list, -12)

country_list <- unique(unlist(list_all_data$bodentien_rueter_negbin$country_id))

prob_gr_zero_model_data <- data.frame(month_id = numeric(), country_id = numeric(), prob_gr_0 = numeric(), model = character())

country_list <- c(57, 70, 246)
month_list <- c(457, 458, 459)

prob_gr_zero_model_data_2 <- data_bodentien_rueter %>%
  filter(month_id %in% month_list, country_id %in% country_list) %>% # keep all rows that are in month_list and country_list
  group_by(month_id, country_id) %>%  # group by month and country (each combination of month_id and country_id is one group)
  summarise( # summarise creates new dataframe and calculates values for each group
    # returns one row for each combination of grouping variables
    prob_gr_0 = 1 - sum(outcome == 0) / n(),
    .groups = "drop"  #  result dataframe is not grouped (mutate is performed for the whole dataframe)
  ) %>%
  mutate(model = list_model_names[1]) # adds new column


