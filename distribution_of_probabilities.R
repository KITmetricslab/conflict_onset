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

# returns one dataframe of all years from one model
process_model <- function(model_path) {
  # all directories of the years
  year_dirs <- list.dirs(file.path(model_path, "cm"), recursive = FALSE)
  
  # 
  outbreak_data <- year_dirs %>%
    map(list.files, full.names = TRUE, pattern = "\\.parquet") %>%  # performs list.files on every dir. of year_dirs; returns the full filepath names
    flatten_chr() %>%                                              # converts list of filepaths into one vector
    map_dfr(arrow::read_parquet)                                   # performs read_parquet on every element in the vector
  
  return(outbreak_data) # dataframe with all parquet files in year_dirs
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



# month country list combination for onset and ongoing peace

# group each list df's by those combinations and plot dist. per model








