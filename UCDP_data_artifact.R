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
library(knitr)

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

# Are fatalities coinsiding for every country-month pair between GED and VIEWS?
monthly_comparison %>%
  filter(difference != 0 | is.na(difference))

# Do all GED country-month pairs exist in VIEWS?
GED_monthly %>%
  anti_join(
    observations_18_23,
    by = c("country_id", "month_id_end" = "month_id")
  )

# Check for duplicate country-month pairs
GED_monthly %>%
  count(country_id, month_id_end) %>%
  filter(n > 1)

observations_18_23 %>%
  count(country_id, month_id) %>%
  filter(n > 1)


# affected country-month-pairs
monthly_affected <- monthly_comparison %>%
  filter(n_multimonth > 0) %>%
  arrange(desc(fatalities_multimonth))

# individual multimonth events
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
















## analysis of VIEWS multimonth aggregation ---------------------------

## -----
## Overall impact on VIEWS actuals
## -----
aggregation_summary <- monthly_comparison %>%
  summarise(
    n_country_months = n(),
    n_affected_months = sum(n_multimonth > 0),
    n_affected_countries = n_distinct(country_id[n_multimonth > 0]),
    n_multimonth_events = sum(n_multimonth),
    fatalities_multimonth = sum(fatalities_multimonth),
    fatalities_total = sum(outcome),
    share_multimonth = 100 * fatalities_multimonth / fatalities_total
  )

aggregation_summary

## ------
## Impact by country table
## -----
country_summary <- monthly_comparison %>%
  group_by(country_id) %>%
  summarise(
    n_affected_months = sum(n_multimonth > 0),
    n_multimonth_events = sum(n_multimonth),
    fatalities_multimonth = sum(fatalities_multimonth),
    fatalities_total = sum(outcome),
    .groups = "drop"
  ) %>%
  mutate(
    share_multimonth = if_else(
      fatalities_total > 0,
      100 * fatalities_multimonth / fatalities_total,
      NA_real_
    )
  ) %>%
  left_join(
    country_id_list %>% select(country_id, name),
    by = "country_id"
  ) %>%
  arrange(desc(fatalities_multimonth))

## Add others row
country_top10_others <- country_summary %>%
  mutate(
    name = if_else(row_number() <= 10, name, "Others")
  ) %>%
  group_by(name) %>%
  summarise(
    n_affected_months = sum(n_affected_months),
    n_multimonth_events = sum(n_multimonth_events),
    fatalities_multimonth = sum(fatalities_multimonth),
    fatalities_total = sum(fatalities_total),
    .groups = "drop"
  ) %>%
  mutate(
    share_multimonth = if_else(
      fatalities_total > 0,
      100 * fatalities_multimonth / fatalities_total,
      NA_real_
    )
  ) %>%
  arrange(name == "Others", desc(fatalities_multimonth))

## Add Total Row
top10_ids <- country_summary %>%
  slice_head(n = 10) %>%
  pull(country_id)

event_durations <- GED_multimonth_details %>%
  left_join(
    country_summary %>% select(country_id, name),
    by = "country_id"
  ) %>%
  mutate(
    name = if_else(country_id %in% top10_ids, name, "Others")
  )

# Add Total group
event_durations_total <- bind_rows(
  event_durations,
  event_durations %>% mutate(name = "Total")
)

total_row <- country_top10_others %>%
  summarise(
    name = "Total",
    n_affected_months = sum(n_affected_months),
    n_multimonth_events = sum(n_multimonth_events),
    fatalities_multimonth = sum(fatalities_multimonth),
    fatalities_total = sum(fatalities_total)
  ) %>%
  mutate(
    share_multimonth =
      100 * fatalities_multimonth / fatalities_total
  )

## Add event duration stats
duration_stats <- event_durations_total %>%
  group_by(name) %>%
  summarise(
    Mean = mean(ged_event_months),
    SD = sd(ged_event_months),
    Max = max(ged_event_months),
    .groups = "drop"
  ) %>%
  mutate(
    duration_stats = sprintf(
      "(%.1f, %.1f, %d)",
      Mean, SD, Max
    )
  ) %>%
  select(name, duration_stats)

## final table
country_summary_final <- bind_rows(
  country_top10_others,
  total_row
) %>%
  left_join(duration_stats, by = "name")


## Export to LaTeX
latex_table <- knitr::kable(
  country_top10_others,
  format = "latex",
  booktabs = TRUE,
  caption = "Top 10 countries and others by multimonth fatalities (2018–2023)",
  row.names = FALSE,
  digits = 1
)

cat(latex_table)

## Distribution of multimonth event duration (not included in paper)
duration_distr <- bind_rows(
  event_durations,
  event_durations %>% mutate(name = "Total")
) %>%
  group_by(name) %>%
  summarise(
    Q1 = quantile(ged_event_months, 0.25),
    Median = median(ged_event_months),
    Mean = mean(ged_event_months),
    Q3 = quantile(ged_event_months, 0.75),
    SD = sd(ged_event_months),
    .groups = "drop"
  )

## -----
## country-month pairs with highest multimonth shares table (not included in paper)
## -----
min_outcome <- 1000

top_country_months <- monthly_affected %>%
  filter(
    outcome >= min_outcome,
    fatalities_multimonth > 0
  ) %>%
  left_join(
    country_id_list %>% select(country_id, name),
    by = "country_id"
  ) %>%
  left_join(month_lookup, by = "month_id") %>%
  mutate(
    label = paste0(name, " | ", format(month_date, "%Y-%m"))
  ) %>%
  arrange(desc(share_multimonth)) %>%
  slice_head(n = 12)

top_country_months %>%
  select(name, month_date, outcome, n_multimonth,
         fatalities_multimonth, share_multimonth)







##########################


## monthly development across test years ----------------------------

time_summary <- GED_multimonth_details %>%
  mutate(
    month_date = lubridate::floor_date(ged_date_end, "month")
  ) %>%
  group_by(month_date) %>%
  summarise(
    n_events = n(),
    fatalities_multimonth = sum(best, na.rm = TRUE),
    fatalities_per_event = mean(best, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(
    month_lookup %>% select(month_date),
    by = "month_date"
  ) %>%
  mutate(
    n_events = coalesce(n_events, 0L),
    fatalities_multimonth = coalesce(fatalities_multimonth, 0)
  ) %>%
  arrange(month_date)

time_long <- time_summary %>%
  pivot_longer(
    cols = c(n_events, fatalities_multimonth, fatalities_per_event),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("n_events", "fatalities_multimonth",
                 "fatalities_per_event"),
      labels = c("Number of events",
                 "Multimonth fatalities",
                 "Average fatalities per event")
    )
  )

p4 <- ggplot(
  time_long,
  aes(x = month_date, y = value)
) +
  geom_col(fill = "steelblue", width = 25) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  labs(
    title = "Temporal distribution of multimonth events",
    subtitle = "Events and fatalities assigned to their ending month",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

p4

############################################################################

## Publication plot settings ------------------------------------------

publication_dpi <- 300
publication_plot_width_px <- 3600

publication_base_size_pt <- 12
publication_axis_title_size_pt <- 13
publication_legend_key_spacing_y_pt <- 2

theme_fontsize <- ggplot2::theme(
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  text = element_text(size = publication_base_size_pt),
  axis.title = element_text(size = publication_axis_title_size_pt),
  axis.text = element_text(size = publication_base_size_pt),
  legend.title = element_text(size = publication_base_size_pt),
  legend.text = element_text(size = publication_base_size_pt),
  legend.key.width = grid::unit(12, "pt"),
  legend.key.height = grid::unit(12, "pt"),
  legend.key.spacing.y = grid::unit(
    publication_legend_key_spacing_y_pt, "pt"
  )
)

# Same colors as fatalities_worldwide_plot
fatality_country_colours <- c(
  "Ethiopia" = "#009E73",
  "Ukraine" = "#E69F00",
  "Afghanistan" = "#D55E00",
  "Yemen" = "#0072B2",
  "Syria" = "#51103C"
)

top5_names <- names(fatality_country_colours)

# Month labels
month_lookup_vec <- setNames(
  format(month_lookup$month_date, "%m-%Y"),
  month_lookup$month_id
)

month_breaks <- seq(
  min(month_lookup$month_id),
  max(month_lookup$month_id),
  by = 6
)

# Year boundaries
month_boundaries <- month_lookup %>%
  filter(
    lubridate::month(month_date) == 12,
    month_id < max(month_id)
  ) %>%
  pull(month_id) + 0.5

# Same numeric formatting as reference plot
number_labels <- function(x) {
  format(
    x,
    big.mark = ".",
    decimal.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


## Prepare plot data --------------------------------------------------

# Identify countries not among the highlighted five
other_ids <- sort(
  setdiff(
    unique(monthly_comparison$country_id),
    country_id_list %>%
      filter(name %in% top5_names) %>%
      pull(country_id)
  )
)

other_keys <- paste0("Others__", other_ids)

# Alternate grey shades for individual other countries
other_colours <- setNames(
  rep(c("grey30", "grey65"), length.out = length(other_keys)),
  other_keys
)

# Colors used in the actual bars
plot_colours <- c(
  "Single-month fatalities" = "white", ##D8C5E5
  other_colours,
  fatality_country_colours,
  "Others" = "grey30"
)

# Stacking order: bottom to top
stack_levels <- c(
  "Single-month fatalities",
  other_keys,
  rev(top5_names),
  "Others"
)

# Prepare country-month data
plot_data <- monthly_comparison %>%
  left_join(
    country_id_list %>%
      select(country_id, name),
    by = "country_id"
  ) %>%
  mutate(
    country_category = if_else(
      name %in% top5_names,
      name,
      paste0("Others__", country_id)
    ),
    fatalities_singlemonth = outcome - fatalities_multimonth
  )

# Validate underlying data
stopifnot(
  !anyNA(plot_data$name),
  !anyNA(plot_data$difference),
  all(plot_data$difference == 0),
  all(plot_data$fatalities_singlemonth >= 0)
)


## Create monthly stacked-bar data ------------------------------------

# Upper panel: number of multimonth events
events_data <- plot_data %>%
  transmute(
    month_id,
    country_category,
    value = n_multimonth
  )

# Lower panel: multimonth fatalities by country
multimonth_fatalities_data <- plot_data %>%
  transmute(
    month_id,
    country_category,
    value = fatalities_multimonth
  )

# Lower panel: all single-month fatalities in one white segment
singlemonth_fatalities_data <- plot_data %>%
  group_by(month_id) %>%
  summarise(
    value = sum(fatalities_singlemonth),
    .groups = "drop"
  ) %>%
  mutate(
    country_category = "Single-month fatalities"
  )

# Combine both components into the SAME monthly bars
fatalities_data <- bind_rows(
  singlemonth_fatalities_data,
  multimonth_fatalities_data
)


## Calculate stacked segment positions --------------------------------

# This ensures every component has exactly the same x position
# and each segment starts where the previous one ends.

make_segments <- function(data) {
  
  data %>%
    mutate(
      country_category = factor(
        country_category,
        levels = stack_levels
      )
    ) %>%
    arrange(month_id, country_category) %>%
    group_by(month_id) %>%
    mutate(
      ymin = cumsum(value) - value,
      ymax = cumsum(value),
      xmin = month_id - 0.45,
      xmax = month_id + 0.45
    ) %>%
    ungroup() %>%
    filter(value > 0)
}

events_segments <- make_segments(events_data)
fatalities_segments <- make_segments(fatalities_data)


## Verify monthly bar totals ------------------------------------------

# The sum of all stacked segments must equal VIEWS actuals
height_check <- fatalities_data %>%
  group_by(month_id) %>%
  summarise(
    reconstructed = sum(value),
    .groups = "drop"
  ) %>%
  left_join(
    observations_18_23 %>%
      group_by(month_id) %>%
      summarise(
        outcome = sum(outcome),
        .groups = "drop"
      ),
    by = "month_id"
  )

stopifnot(
  !anyNA(height_check$outcome),
  all(height_check$reconstructed == height_check$outcome)
)


## Shared plot theme --------------------------------------------------

theme_multimonth <- theme_minimal() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks = element_line(colour = "black"),
    legend.title = element_blank()
  ) +
  theme_fontsize


# Identical x-axis for both panels
month_axis <- scale_x_continuous(
  breaks = month_breaks,
  labels = month_lookup_vec[as.character(month_breaks)]
)


## Upper panel: multimonth events -------------------------------------

p_events <- ggplot(
  events_segments,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = country_category
  )
) +
  geom_rect(colour = NA) +
  geom_vline(
    xintercept = month_boundaries,
    linetype = "dashed",
    colour = "black"
  ) +
  scale_fill_manual(
    values = plot_colours,
    guide = "none"
  ) +
  month_axis +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    labels = number_labels,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Ending multimonth events"
  ) +
  theme_multimonth

# Both panels retain their month tick labels
# but only the lower panel receives the x-axis title.


## Lower panel: all fatalities ---------------------------------------

p_fatalities <- ggplot(
  fatalities_segments,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = country_category
  )
) +
  ggpattern::geom_rect_pattern(
    aes(
      pattern = country_category == "Single-month fatalities"
    ),
    colour = NA,
    pattern_fill = "black",
    pattern_colour = NA,
    pattern_angle = 45,
    pattern_density = 0.21, #0.19
    pattern_spacing = 0.007
  ) +
  ggpattern::scale_pattern_manual(
    values = c("TRUE" = "stripe", "FALSE" = "none"),
    guide = "none"
  ) +
  geom_vline(
    xintercept = month_boundaries,
    linetype = "dashed",
    colour = "black"
  ) +
  scale_fill_manual(
    values = plot_colours,
    guide = "none"
  ) +
  month_axis +
  scale_y_continuous(
    breaks = seq(0, 150000, length.out = 6),
    labels = number_labels,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Month",
    y = "Fatalities"
  ) +
  theme_multimonth

## Create correctly colored shared legend ----------------------------

# Custom legend key:
# - highlighted countries: their original colors
# - Others: dark/light grey split
# - Single-month: white with black border

draw_key_multimonth <- function(data, params, size) {
  
  fill_value <- tolower(as.character(data$fill[1]))
  
  # Others: two grey shades
  if (fill_value %in% c("grey30", "#4d4d4d")) {
    return(grid::grobTree(
      grid::rectGrob(
        y = 0.75,
        height = 0.5,
        gp = grid::gpar(fill = "grey30", col = NA)
      ),
      grid::rectGrob(
        y = 0.25,
        height = 0.5,
        gp = grid::gpar(fill = "grey65", col = NA)
      )
    ))
  }
  
  # Single-month fatalities: white with black diagonal stripes
  if (fill_value %in% c("white", "#ffffff")) {
    
    stripe_x <- seq(-0.75, 0.75, by = 0.25)
    
    return(grid::grobTree(
      grid::rectGrob(
        gp = grid::gpar(fill = "white", col = NA)
      ),
      grid::segmentsGrob(
        x0 = grid::unit(stripe_x, "npc"),
        y0 = grid::unit(0, "npc"),
        x1 = grid::unit(stripe_x + 1, "npc"),
        y1 = grid::unit(1, "npc"),
        gp = grid::gpar(
          col = "black",
          lwd = 0.9
        )
      ),
      vp = grid::viewport(clip = "on")
    ))
  }
  
  # Highlighted countries: original colors
  grid::rectGrob(
    gp = grid::gpar(
      fill = data$fill[1],
      col = NA
    )
  )
}


# Create legend separately to prevent grey legend-key problems
legend_levels <- c(
  top5_names,
  "Others",
  "Single-month fatalities"
)

legend_colours <- c(
  fatality_country_colours,
  "Others" = "grey30",
  "Single-month fatalities" = "white"
)

legend_data <- tibble(
  country_category = factor(
    legend_levels,
    levels = legend_levels
  ),
  value = rep(1, length(legend_levels))
)

legend_source <- ggplot(
  legend_data,
  aes(
    x = 1,
    y = value,
    fill = country_category
  )
) +
  geom_col(
    key_glyph = draw_key_multimonth,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = legend_colours,
    breaks = legend_levels,
    limits = legend_levels,
    drop = FALSE,
    name = NULL
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      override.aes = list(alpha = 1)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(
      size = publication_base_size_pt
    ),
    legend.key.width = grid::unit(12, "pt"),
    legend.key.height = grid::unit(12, "pt"),
    legend.key.spacing.y = grid::unit(
      publication_legend_key_spacing_y_pt, "pt"
    )
  )

shared_legend <- cowplot::get_legend(legend_source)


## Combine aligned panels ---------------------------------------------

# Use patchwork for consistent alignment of y-axis titles
# and plotting areas, as in the VIEWS figures

combined_panels <- (
  p_events /
    p_fatalities
) +
  patchwork::plot_layout(
    heights = c(1, 1.4)
  )

# Add the separately generated legend
p4 <- cowplot::plot_grid(
  combined_panels,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.19)
)


## Export for LaTeX ---------------------------------------------------

ggplot2::ggsave(
  filename = "Plot_VIEWS_multimonth.png",
  plot = p4,
  width = 3600,
  height = 2400,
  units = "px",
  dpi = 300,
  device = ragg::agg_png,
  bg = "white"
)







##############################################################################









##################################################################################

## Scatter plot: fatalities per multimonth event -----------------------

# Top 5 countries by total multimonth fatalities
top5_countries <- country_summary %>%
  slice_max(
    order_by = fatalities_multimonth,
    n = 5,
    with_ties = FALSE
  ) %>%
  select(country_id, name)


## Colors -------------------------------------------------------------

# Retain the country colors from fatalities_worldwide_plot.
# Azerbaijan receives an additional color.
scatter_colours <- c(
  fatality_country_colours,
  "Azerbaijan" = "#CC79A7",
  "Other countries" = "grey70"
)

scatter_legend_levels <- c(
  top5_countries$name,
  "Other countries"
)

stopifnot(
  all(top5_countries$name %in% names(scatter_colours))
)


## Prepare scatter data -----------------------------------------------

scatter_data <- GED_multimonth_details %>%
  filter(best > 0) %>%
  left_join(
    top5_countries,
    by = "country_id"
  ) %>%
  mutate(
    highlight_country = if_else(
      is.na(name),
      "Other countries",
      name
    ),
    highlight_country = factor(
      highlight_country,
      levels = scatter_legend_levels
    )
  )

# Use identical y-axis and point-size scales in both panels.
scatter_y_max <- max(
  200000,
  max(scatter_data$best, na.rm = TRUE)
)

scatter_size_max <- max(
  12,
  max(scatter_data$ged_event_months, na.rm = TRUE)
)


## Function to create one 3-year panel -------------------------------

make_scatter_panel <- function(
    first_month,
    last_month,
    show_x_title = FALSE
) {
  
  # Same six-month labeling as the VIEWS figures
  x_breaks <- month_breaks[
    month_breaks >= first_month &
      month_breaks <= last_month
  ]
  
  # Only show year boundaries within the panel
  year_lines <- month_boundaries[
    month_boundaries > first_month &
      month_boundaries < last_month
  ]
  
  panel_data <- scatter_data %>%
    filter(
      dplyr::between(
        month_id_end,
        first_month,
        last_month
      )
    )
  
  ggplot(
    panel_data,
    aes(
      x = month_id_end,
      y = best,
      size = ged_event_months,
      colour = highlight_country
    )
  ) +
    
    # Events assigned to their ending month
    geom_point(
      position = position_jitter(
        width = 0.20,
        height = 0,
        seed = 42
      ),
      alpha = 0.90
    ) +
    
    # Dashed year boundaries
    geom_vline(
      xintercept = year_lines,
      linetype = "dashed",
      colour = "black"
    ) +
    
    # Country colors
    scale_colour_manual(
      values = scatter_colours,
      breaks = scatter_legend_levels,
      limits = scatter_legend_levels,
      drop = FALSE,
      name = NULL
    ) +
    
    # Same point-size encoding in both panels
    scale_size_continuous(
      name = "Event duration\n(months)",
      limits = c(2, scatter_size_max),
      range = c(1.4, 5.8),
      breaks = c(2, 4, 6, 8, 10, 12)
    ) +
    
    # Same month format as p_fatalities
    scale_x_continuous(
      breaks = x_breaks,
      labels = month_lookup_vec[as.character(x_breaks)],
      limits = c(
        first_month - 0.5,
        last_month + 0.5
      ),
      expand = expansion(mult = 0.05)
    ) +
    
    # Identical pseudo-log scale for both panels
    scale_y_log10(
      limits = c(0.3, scatter_y_max),
      breaks = c(
        1, 10, 100, 1000,
        10000, 100000
      ),
      labels = number_labels,
      expand = expansion(mult = c(0, 0.05))
    ) +
    
    # Legend design
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          colour = "grey60",
          alpha = 1
        )
      ),
      colour = guide_legend(
        order = 2,
        override.aes = list(
          size = 4.8,
          alpha = 1
        )
      )
    ) +
    
    labs(
      x = if (show_x_title) "Ending month" else NULL,
      y = "Fatalities per event (log10 scale)"
    ) +
    
    theme_multimonth +
    
    theme(
      # Restore the title of the event-duration legend
      legend.title = element_text(
        size = publication_base_size_pt
      ),
      
      legend.text = element_text(
        size = publication_base_size_pt
      ),
      
      legend.key.spacing.y = grid::unit(
        publication_legend_key_spacing_y_pt,
        "pt"
      ),
      
      # No year headings
      strip.text = element_blank(),
      
      # Keep grid lines horizontal only
      panel.grid.minor = element_blank(),
      
      legend.position = "right"
    )
}


## Upper panel: January 2018 – December 2020 -------------------------

p_scatter_18_20 <- make_scatter_panel(
  first_month = 457,
  last_month = 492,
  show_x_title = FALSE
)


## Lower panel: January 2021 – December 2023 -------------------------

p_scatter_21_23 <- make_scatter_panel(
  first_month = 493,
  last_month = 528,
  show_x_title = TRUE
)


## Combine panels -----------------------------------------------------

p3 <- (
  p_scatter_18_20 /
    p_scatter_21_23
) +
  patchwork::plot_layout(
    guides = "collect",
    axis_titles = "collect",
    heights = c(1, 1)
  ) &
  theme(
    legend.position = "right"
  )

# Display
p3


## Export for LaTeX ---------------------------------------------------

ggplot2::ggsave(
  filename = "Plot_VIEWS_multimonth_scatter.png",
  plot = p3,
  width = publication_plot_width_px,
  height = 2400,
  units = "px",
  dpi = publication_dpi,
  device = ragg::agg_png,
  bg = "white"
)


# Only multimonth events with positive fatality estimates (best > 0) are displayed to allow for a logarithmic y-axis.

GED_multimonth_details %>%
  summarise(
    n_total = n(),
    n_positive = sum(best > 0, na.rm = TRUE),
    n_zero = sum(best == 0, na.rm = TRUE),
    n_missing = sum(is.na(best))
  )

# 138 sind 0
