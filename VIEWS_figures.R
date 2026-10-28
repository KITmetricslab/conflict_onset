# ==============================================================================
# VIEWS competition: publication figures
# ==============================================================================
#
# Purpose
# -------
# Reproduce only the eight figures used in the manuscript:
#
#   1. Plot_VIEWS_fatalities.png
#   2. Plot_VIEWS_triple_country_greq25.png
#   3. Plot_VIEWS_triple_country_top5_per_measure_greq25.png
#   4. Plot_VIEWS_triple_situation_greq25.png
#   5. Plot_VIEWS_double_thrBS_greq25.png
#   6. Plot_VIEWS_ROC_greq25.png
#   7. Plot_VIEWS_reliability_greq25.png
#   8. Plot_VIEWS_triple_model_ranking_greq25.png
#
# The script deliberately excludes exploratory calculations, log-score analyses,
# score-integral checks, model diagnostics, and figures not listed above.
#
# Required external data
# ----------------------
# The script automatically uses the KIT Windows/macOS network-share path
# containing:
#   cm_actuals_2018.parquet, ..., cm_actuals_2023.parquet
#   cm_features.parquet
#   all_available_cm_predictions/<model>/cm/window=Y20YY/<model>_cm_Y20YY.parquet
#
# Required repository files
# -------------------------
# country_list.csv
#
# Required cached score objects
# -----------------------------
# By default these are read from ./output. Override with VIEWS_CACHE_DIR.
#
#   models_crps.RData
#     object: models_crps
#
#   models_twcrps.RData
#     object: models_twcrps
#
#   all_models_brier_threshold_extended.RData
#     object: all_models_brier_full
#
# The first two caches contain country-month CRPS values. The third contains
# the mean Brier score as a function of threshold. Keeping these published
# intermediates avoids recomputing expensive score objects in the figure script.
#
# Optional environment variables
# ------------------------------
# Repository-relative paths:
#   cache directory:      output/
#   figure directory:     plots_final/
#   country lookup:       country_list.csv
#
# ==============================================================================


# ---- Packages ----------------------------------------------------------------

required_packages <- c(
  "arrow",
  "dplyr",
  "ggplot2",
  "readr",
  "purrr",
  "precrec",
  "reliabilitydiag",
  "triptych",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(precrec)
  library(reliabilitydiag)
  library(triptych)
  library(patchwork)
})


# ---- Configuration -----------------------------------------------------------

## VIEWS data location.
## This reproduces the platform-specific paths used in the original analysis.
os <- Sys.info()[["sysname"]]

data_dir <- if (os == "Windows") {
  "//stat-meth-file1.stat.kit.edu/share-alle/Data/VIEWS/"
} else if (os == "Darwin") {
  "/Volumes/share-alle/Data/VIEWS/"
} else {
  stop(
    "No default VIEWS data path is configured for operating system: ",
    os
  )
}

## Resolve repository-relative paths from the location of this script rather
## than from the process working directory. This makes the script robust when
## run via RStudio, Rscript, or from another directory.
get_script_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  
  if (length(file_arg) > 0L) {
    return(
      dirname(
        normalizePath(
          sub("^--file=", "", file_arg[[1L]]),
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }
  
  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()
  ) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    
    if (nzchar(active_path)) {
      return(
        dirname(
          normalizePath(
            active_path,
            winslash = "/",
            mustWork = TRUE
          )
        )
      )
    }
  }
  
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  )
}

script_dir <- get_script_dir()

cache_dir <- file.path(script_dir, "output")
publication_plot_dir <- file.path(script_dir, "plots_final")
country_list_path <- file.path(script_dir, "country_list.csv")

if (!dir.exists(data_dir)) {
  stop(
    "VIEWS data directory not found: ",
    data_dir,
    "\nPlease mount the KIT share-alle network drive first."
  )
}

dir.create(publication_plot_dir, recursive = TRUE, showWarnings = FALSE)

round_samples <- TRUE

# Manuscript definition: conflict at threshold k = 24 means fatalities > 24.
fatality_threshold <- 24L
greq_id <- "greq25"

selected_models <- c(
  "bodentien_rueter_negbin",
  "conflictforecast_v2",
  "submission_muchlinski_thornhill",
  "ShapeFinder",
  "submission_final_omm",
  "unito_transformer",
  "conflictology",
  "zero"
)

model_labels <- c(
  "bodentien_rueter_negbin" = "BR NB",
  "conflictforecast_v2" = "CFLT RF",
  "submission_muchlinski_thornhill" = "MT ZeroInfl GAM",
  "ShapeFinder" = "PACE ShapeFinder",
  "submission_final_omm" = "RV O MM",
  "unito_transformer" = "UNITO NB Transformer",
  "conflictology" = "VIEWS Conflictology",
  "zero" = "VIEWS Zero"
)

selected_colors <- c(
  "#009682",
  "#4664aa",
  "#23a1e0",
  "black",
  "#df9b1b",
  "#8cb63c",
  "#a22223",
  "#a3107c"
)

selected_model_colors <- setNames(selected_colors, selected_models)
selected_model_labels <- model_labels[selected_models]

mean_brier_score_label <- bquote(
  "Mean Brier Score for" ~ italic(k) == .(fatality_threshold)
)


# ---- Publication export settings --------------------------------------------

publication_dpi <- 300
publication_plot_width_px <- 3600
publication_height_scale <- 2 / 3

single_panel_height_px <- round((4 / 5) * 1.2 * 1322)
triple_country_height_px <- 3 * single_panel_height_px
triple_situation_height_px <- round(
  3900 * publication_plot_width_px / 4222
)
double_threshold_height_px <- round(
  2600 * publication_plot_width_px / 4222
)
roc_plot_height_px <- round(
  1470 * publication_plot_width_px / 3400
)
reliability_plot_height_px <- round(
  1800 * publication_plot_width_px / 3400
)
model_ranking_height_px <- round(
  (3 / 5) * 2000 * publication_plot_width_px / 4222
)
triple_model_ranking_height_px <- 3 * model_ranking_height_px

publication_base_size_pt <- 12
publication_axis_title_size_pt <- 13
publication_panel_title_size_pt <- 12

## Use exactly the same legend-key dimensions in every publication figure.
publication_legend_key_width_pt <- 12
publication_legend_key_height_pt <- 12
publication_legend_key_spacing_y_pt <- 2

pt_to_mm <- 1 / 2.845276
publication_geom_text_size <- publication_base_size_pt * pt_to_mm

theme_no_plot_title <- ggplot2::theme(
  plot.title = ggplot2::element_blank(),
  plot.subtitle = ggplot2::element_blank()
)

theme_publication_text <- ggplot2::theme(
  text = ggplot2::element_text(size = publication_base_size_pt),
  axis.title = ggplot2::element_text(size = publication_axis_title_size_pt),
  axis.text = ggplot2::element_text(size = publication_base_size_pt),
  legend.title = ggplot2::element_text(size = publication_base_size_pt),
  legend.text = ggplot2::element_text(size = publication_base_size_pt),
  legend.key.width = grid::unit(publication_legend_key_width_pt, "pt"),
  legend.key.height = grid::unit(publication_legend_key_height_pt, "pt"),
  legend.key.spacing.y = grid::unit(publication_legend_key_spacing_y_pt, "pt"),
  strip.text = ggplot2::element_text(size = publication_base_size_pt)
)

theme_fontsize <- ggplot2::theme(
  plot.title = ggplot2::element_blank(),
  plot.subtitle = ggplot2::element_blank(),
  text = ggplot2::element_text(size = publication_base_size_pt),
  axis.title = ggplot2::element_text(size = publication_axis_title_size_pt),
  axis.text = ggplot2::element_text(size = publication_base_size_pt),
  legend.title = ggplot2::element_text(size = publication_base_size_pt),
  legend.text = ggplot2::element_text(size = publication_base_size_pt),
  legend.key.width = grid::unit(publication_legend_key_width_pt, "pt"),
  legend.key.height = grid::unit(publication_legend_key_height_pt, "pt"),
  legend.key.spacing.y = grid::unit(publication_legend_key_spacing_y_pt, "pt"),
  strip.text = ggplot2::element_text(size = publication_base_size_pt)
)

extract_grob_labels <- function(grob) {
  labels <- character()
  
  if (!is.null(grob$label)) {
    labels <- c(labels, as.character(grob$label))
  }
  
  if (!is.null(grob$children)) {
    for (child_name in names(grob$children)) {
      labels <- c(labels, extract_grob_labels(grob$children[[child_name]]))
    }
  }
  
  if (!is.null(grob$grobs)) {
    for (child_grob in grob$grobs) {
      labels <- c(labels, extract_grob_labels(child_grob))
    }
  }
  
  labels
}

apply_split_others_legend <- function(
    plot_object,
    others_label = "Others",
    top_fill = "grey30",
    bottom_fill = "grey65"
) {
  # Convert the finished plot to a grob and replace only the legend key that
  # belongs to "Others". This leaves all ordinary country legend keys intact.
  plot_grob <- if (inherits(plot_object, "patchwork")) {
    patchwork::patchworkGrob(plot_object)
  } else {
    ggplot2::ggplotGrob(plot_object)
  }
  
  guide_box_indices <- which(grepl("^guide-box", plot_grob$layout$name))
  modified_any <- FALSE
  
  for (guide_box_idx in guide_box_indices) {
    guide_box <- plot_grob$grobs[[guide_box_idx]]
    
    # Empty guide positions are zeroGrobs; only actual guide gtables matter.
    if (!inherits(guide_box, "gtable")) {
      next
    }
    
    for (guide_i in seq_along(guide_box$grobs)) {
      guide_table <- guide_box$grobs[[guide_i]]
      
      if (!inherits(guide_table, "gtable")) {
        next
      }
      
      label_indices <- which(grepl("label", guide_table$layout$name))
      
      if (length(label_indices) == 0L) {
        next
      }
      
      others_label_indices <- label_indices[
        vapply(
          label_indices,
          function(idx) {
            others_label %in% extract_grob_labels(guide_table$grobs[[idx]])
          },
          logical(1)
        )
      ]
      
      if (length(others_label_indices) == 0L) {
        next
      }
      
      target_rows <- unique(guide_table$layout$t[others_label_indices])
      
      for (target_row in target_rows) {
        key_indices <- which(
          grepl("^key", guide_table$layout$name) &
            guide_table$layout$t == target_row
        )
        
        for (key_idx in key_indices) {
          # Split the normal legend key horizontally: dark grey on top,
          # light grey on the bottom.
          guide_table$grobs[[key_idx]] <- grid::grobTree(
            grid::rectGrob(
              x = 0.5,
              y = 0.75,
              width = 1,
              height = 0.5,
              gp = grid::gpar(fill = top_fill, col = NA)
            ),
            grid::rectGrob(
              x = 0.5,
              y = 0.25,
              width = 1,
              height = 0.5,
              gp = grid::gpar(fill = bottom_fill, col = NA)
            )
          )
        }
        
        modified_any <- TRUE
      }
      
      guide_box$grobs[[guide_i]] <- guide_table
    }
    
    plot_grob$grobs[[guide_box_idx]] <- guide_box
  }
  
  if (modified_any) {
    plot_grob
  } else {
    plot_object
  }
}

# Legend key used for country stacks. Ordinary categories get the standard
# rectangular key. The visible "Others" category uses grey30 in the scale and
# is shown as one normal-sized key split horizontally into dark and light grey.
draw_key_country_fill <- function(data, params, size) {
  fill_value <- if (is.null(data$fill)) NA_character_ else as.character(data$fill)
  is_others_key <- !is.na(fill_value) &&
    toupper(fill_value) %in% c("GREY30", "#4D4D4D")
  
  if (!is_others_key) {
    return(ggplot2::draw_key_rect(data, params, size))
  }
  
  grid::grobTree(
    grid::rectGrob(
      x = 0.5, y = 0.75,
      width = grid::unit(1, "npc"),
      height = grid::unit(0.5, "npc"),
      gp = grid::gpar(fill = "grey30", col = NA)
    ),
    grid::rectGrob(
      x = 0.5, y = 0.25,
      width = grid::unit(1, "npc"),
      height = grid::unit(0.5, "npc"),
      gp = grid::gpar(fill = "grey65", col = NA)
    )
  )
}

save_plot_for_latex <- function(
    filename,
    plot_object,
    height_px,
    height_scale = publication_height_scale
) {
  plot_for_export <- if (inherits(plot_object, "patchwork")) {
    plot_object & theme_publication_text
  } else {
    plot_object + theme_publication_text
  }
  
  ggplot2::ggsave(
    filename = filename,
    plot = plot_for_export,
    width = publication_plot_width_px,
    height = round(height_scale * height_px),
    dpi = publication_dpi,
    units = "px",
    bg = "white"
  )
}


# ---- Validation helpers ------------------------------------------------------

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  invisible(path)
}

require_columns <- function(data, columns, object_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(
      object_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(data)
}

assert_unique_key <- function(data, key, object_name) {
  duplicates <- data %>%
    group_by(across(all_of(key))) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n > 1L)
  
  if (nrow(duplicates) > 0L) {
    stop(
      object_name,
      " contains duplicated ",
      paste(key, collapse = " + "),
      " keys."
    )
  }
  invisible(data)
}

load_rdata_object <- function(path, object_name) {
  require_file(path)
  env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = env)
  
  if (!object_name %in% loaded_names) {
    stop(
      "Expected object '",
      object_name,
      "' in ",
      path,
      ". Found: ",
      paste(loaded_names, collapse = ", ")
    )
  }
  
  env[[object_name]]
}

validate_named_model_list <- function(x, object_name) {
  if (is.null(names(x))) {
    stop(object_name, " must be a named list.")
  }
  
  missing_models <- setdiff(selected_models, names(x))
  if (length(missing_models) > 0L) {
    stop(
      object_name,
      " is missing selected model(s): ",
      paste(missing_models, collapse = ", ")
    )
  }
  
  invisible(x)
}


# ---- Observations and country lookup -----------------------------------------

actual_files <- file.path(
  data_dir,
  paste0("cm_actuals_", 2018:2023, ".parquet")
)

invisible(lapply(actual_files, require_file))

observations_18_23 <- purrr::map_dfr(
  actual_files,
  arrow::read_parquet
) %>%
  dplyr::select(country_id, month_id, outcome) %>%
  dplyr::arrange(country_id, month_id)

require_columns(
  observations_18_23,
  c("country_id", "month_id", "outcome"),
  "observations_18_23"
)
assert_unique_key(
  observations_18_23,
  c("country_id", "month_id"),
  "observations_18_23"
)

if (anyNA(observations_18_23$outcome)) {
  stop("observations_18_23 contains missing outcomes.")
}

if (!file.exists(country_list_path)) {
  stop(
    "country_list.csv not found next to the analysis script. Expected: ",
    country_list_path,
    "\nPlace country_list.csv in the same repository directory as this R script."
  )
}

country_id_list <- readr::read_csv(
  country_list_path,
  show_col_types = FALSE
) %>%
  dplyr::select(country_id, name)

require_columns(country_id_list, c("country_id", "name"), "country_id_list")
assert_unique_key(country_id_list, "country_id", "country_id_list")

if (anyNA(country_id_list$name)) {
  stop("country_list.csv contains missing country names.")
}

month_ids <- sort(unique(observations_18_23$month_id))

if (
  length(month_ids) != 72L ||
  !all(diff(month_ids) == 1L)
) {
  stop(
    "Expected 72 consecutive monthly IDs for January 2018 to December 2023."
  )
}

month_calendar <- data.frame(
  month_id = month_ids,
  date = seq(
    from = as.Date("2018-01-01"),
    to = as.Date("2023-12-01"),
    by = "month"
  )
) %>%
  mutate(
    month_label = format(date, "%m-%Y"),
    month_number = as.integer(format(date, "%m"))
  )

month_lookup_vec <- setNames(
  month_calendar$month_label,
  month_calendar$month_id
)

observation_keys <- observations_18_23 %>%
  distinct(country_id, month_id)


# ---- Predictive event probabilities for the eight selected models ------------

prediction_year_suffixes <- sprintf("%02d", 18:23)

prediction_paths_for_model <- function(model_name) {
  file.path(
    data_dir,
    "all_available_cm_predictions",
    model_name,
    "cm",
    paste0("window=Y20", prediction_year_suffixes),
    paste0(model_name, "_cm_Y20", prediction_year_suffixes, ".parquet")
  )
}

compute_onset_probabilities <- function(model_name) {
  paths <- prediction_paths_for_model(model_name)
  invisible(lapply(paths, require_file))
  
  predictions <- purrr::map_dfr(
    paths,
    arrow::read_parquet
  ) %>%
    dplyr::select(country_id, month_id, outcome) %>%
    semi_join(observation_keys, by = c("country_id", "month_id"))
  
  if (round_samples) {
    predictions <- predictions %>%
      mutate(outcome = round(outcome, 0))
  }
  
  if (anyNA(predictions$outcome)) {
    stop("Prediction outcomes contain NA for model: ", model_name)
  }
  
  draw_counts <- predictions %>%
    count(country_id, month_id, name = "n_draws")
  
  if (any(draw_counts$n_draws != 1000L)) {
    stop(
      "Expected exactly 1000 predictive draws per country-month for model: ",
      model_name
    )
  }
  
  missing_keys <- anti_join(
    observation_keys,
    draw_counts %>% select(country_id, month_id),
    by = c("country_id", "month_id")
  )
  
  if (nrow(missing_keys) > 0L) {
    stop("Predictions are incomplete for model: ", model_name)
  }
  
  predictions %>%
    group_by(country_id, month_id) %>%
    summarise(
      onset_prob_pred = mean(outcome > fatality_threshold),
      .groups = "drop"
    )
}

message("Computing event probabilities for the eight selected models...")
onset_probabilities <- setNames(
  lapply(selected_models, compute_onset_probabilities),
  selected_models
)


# ---- Load cached CRPS values and construct one analysis table ----------------

models_crps <- load_rdata_object(
  file.path(cache_dir, "models_crps.RData"),
  "models_crps"
)
models_twcrps <- load_rdata_object(
  file.path(cache_dir, "models_twcrps.RData"),
  "models_twcrps"
)

validate_named_model_list(models_crps, "models_crps")
validate_named_model_list(models_twcrps, "models_twcrps")

models_crps <- models_crps[selected_models]
models_twcrps <- models_twcrps[selected_models]

score_long <- purrr::map_dfr(
  selected_models,
  function(model_name) {
    crps_df <- models_crps[[model_name]] %>%
      dplyr::select(country_id, month_id, crps)
    
    twcrps_df <- models_twcrps[[model_name]] %>%
      dplyr::select(country_id, month_id, twcrps)
    
    assert_unique_key(
      crps_df,
      c("country_id", "month_id"),
      paste0("CRPS data for ", model_name)
    )
    assert_unique_key(
      twcrps_df,
      c("country_id", "month_id"),
      paste0("twCRPS data for ", model_name)
    )
    
    out <- observation_keys %>%
      left_join(crps_df, by = c("country_id", "month_id")) %>%
      left_join(twcrps_df, by = c("country_id", "month_id")) %>%
      left_join(
        onset_probabilities[[model_name]],
        by = c("country_id", "month_id")
      ) %>%
      left_join(
        observations_18_23 %>%
          select(country_id, month_id, actual = outcome),
        by = c("country_id", "month_id")
      ) %>%
      mutate(
        model = model_name,
        actual_conflict = actual > fatality_threshold,
        brier_onset = (as.numeric(actual_conflict) - onset_prob_pred)^2
      )
    
    if (anyNA(out[, c("crps", "twcrps", "onset_prob_pred", "actual", "brier_onset")])) {
      stop("Incomplete score data for model: ", model_name)
    }
    
    out
  }
)

# One global country ordering for every grey "Others" block.
# Overall CRPS is the mean CRPS across all months and selected models for each
# country. This ranking is computed once and reused in the fatalities figure
# and in every country-score panel.
global_country_crps_order <- score_long %>%
  group_by(country_id) %>%
  summarise(
    overall_crps = mean(crps),
    .groups = "drop"
  ) %>%
  arrange(desc(overall_crps), country_id) %>%
  left_join(country_id_list, by = "country_id")

if (anyNA(global_country_crps_order$name)) {
  stop("Missing country names in the global CRPS ordering.")
}

global_crps_country_ids <- global_country_crps_order$country_id

rm(models_crps, models_twcrps, onset_probabilities)


# ---- Conflict situation at the annual forecast origin ------------------------
#
# For each calendar year, the reference state is conflict status in the
# preceding October. Thus:
#   peace         = no conflict now, no conflict at preceding October
#   onset         = conflict now, no conflict at preceding October
#   conflict      = conflict now, conflict at preceding October
#   deescalation  = no conflict now, conflict at preceding October
#
# "Previous conflict" and "No previous conflict" in the ranking therefore refer
# to this preceding-October state.

features_path <- file.path(data_dir, "cm_features.parquet")
require_file(features_path)

observations_2017 <- arrow::read_parquet(features_path) %>%
  select(month_id, country_id, ged_sb) %>%
  filter(month_id %in% 445:456) %>%
  transmute(
    country_id,
    month_id,
    outcome = ged_sb
  )

historical_outcomes <- bind_rows(
  observations_2017,
  observations_18_23 %>% select(country_id, month_id, outcome)
)

first_forecast_month_id <- min(month_ids)
first_reference_october_id <- first_forecast_month_id - 3L

conflict_situations <- observations_18_23 %>%
  transmute(
    country_id,
    month_id,
    current_conflict = outcome > fatality_threshold,
    reference_october_id =
      first_reference_october_id +
      12L * ((month_id - first_forecast_month_id) %/% 12L)
  ) %>%
  left_join(
    historical_outcomes %>%
      transmute(
        country_id,
        reference_october_id = month_id,
        previous_conflict = outcome > fatality_threshold
      ),
    by = c("country_id", "reference_october_id")
  ) %>%
  mutate(
    situation_month = case_when(
      !current_conflict & !previous_conflict ~ "peace",
      current_conflict & !previous_conflict ~ "onset",
      current_conflict & previous_conflict ~ "conflict",
      !current_conflict & previous_conflict ~ "deescalation",
      TRUE ~ NA_character_
    )
  ) %>%
  select(country_id, month_id, situation_month)

if (anyNA(conflict_situations$situation_month)) {
  stop("Could not classify every country-month into a conflict situation.")
}

score_long <- score_long %>%
  left_join(
    conflict_situations,
    by = c("country_id", "month_id")
  )

if (anyNA(score_long$situation_month)) {
  stop("Conflict-situation join failed for at least one score row.")
}

previous_no_conflict_data <- score_long %>%
  filter(situation_month %in% c("peace", "onset")) %>%
  mutate(actual_binary = as.integer(actual_conflict))


# ---- Shared calendar annotations ---------------------------------------------

month_boundaries <- month_calendar %>%
  filter(month_number == 12L, month_id < max(month_id)) %>%
  pull(month_id) + 0.5

country_breaks <- seq(
  min(month_ids),
  max(month_ids),
  by = 6
)


# ==============================================================================
# Figure 1: worldwide fatalities
# ==============================================================================

top_country_fatalities <- observations_18_23 %>%
  group_by(country_id) %>%
  summarise(
    total_fatalities = sum(outcome),
    .groups = "drop"
  ) %>%
  arrange(desc(total_fatalities), country_id) %>%
  slice_head(n = 5L) %>%
  left_join(country_id_list, by = "country_id")

if (anyNA(top_country_fatalities$name)) {
  stop("Missing country names for the top-five fatality countries.")
}

top_fatality_country_ids <- top_country_fatalities$country_id
top_fatality_country_names <- top_country_fatalities$name

# Fixed publication palette. The assertion guards against silently assigning
# these colours to a different set of countries if the source data change.
fatality_country_colours <- c(
  "Ethiopia" = "#009E73",
  "Ukraine" = "#E69F00",
  "Afghanistan" = "#D55E00",
  "Yemen" = "#0072B2",
  "Syria" = "#51103C",
  "Others" = "grey30"
)

if (!setequal(
  top_fatality_country_names,
  setdiff(names(fatality_country_colours), "Others")
)) {
  stop(
    "The five countries with the largest total fatalities no longer match ",
    "the fixed publication palette. Found: ",
    paste(top_fatality_country_names, collapse = ", ")
  )
}

# Keep "Others" as one legend entry, but show its internal country
# composition in the bars. The grey countries follow the single global CRPS
# ordering used throughout the manuscript figures: largest overall CRPS first,
# then descending overall CRPS. Highlighted countries are removed while the
# relative order of all remaining countries is preserved.
fatality_other_country_order <- global_country_crps_order %>%
  filter(!country_id %in% top_fatality_country_ids)

if (anyNA(fatality_other_country_order$name)) {
  stop("Missing country names among countries grouped as Others.")
}

fatality_other_keys <- paste0(
  "Others__",
  fatality_other_country_order$country_id
)

fatality_other_colours <- setNames(
  rep(
    c("grey30", "grey65"),
    length.out = length(fatality_other_keys)
  ),
  fatality_other_keys
)

fatality_highlight_colours <- fatality_country_colours[
  names(fatality_country_colours) != "Others"
]

# "Others" itself is an unused dummy level used only for one legend key.
# The plotted grey block is drawn country by country.
fatality_stack_colours <- c(
  fatality_highlight_colours,
  "Others" = "grey30",
  fatality_other_colours
)

fatality_stack_levels <- c(
  names(fatality_highlight_colours),
  fatality_other_keys,
  "Others"
)

fatalities_worldwide_plot_data <- observations_18_23 %>%
  left_join(country_id_list, by = "country_id") %>%
  mutate(
    country_category = if_else(
      country_id %in% top_fatality_country_ids,
      name,
      paste0("Others__", country_id)
    ),
    country_category = factor(
      country_category,
      levels = fatality_stack_levels
    )
  ) %>%
  group_by(month_id, country_category) %>%
  summarise(
    n_fatalities = sum(outcome),
    .groups = "drop"
  ) %>%
  bind_rows(
    tibble(
      month_id = min(month_ids),
      country_category = factor("Others", levels = fatality_stack_levels),
      n_fatalities = 0
    )
  )

october_tops <- fatalities_worldwide_plot_data %>%
  left_join(
    month_calendar %>% select(month_id, month_number),
    by = "month_id"
  ) %>%
  filter(month_number == 10L) %>%
  group_by(month_id) %>%
  summarise(
    y_pos = sum(n_fatalities),
    .groups = "drop"
  )

fatalities_worldwide_plot <- ggplot(
  fatalities_worldwide_plot_data,
  aes(
    x = month_id,
    y = n_fatalities,
    fill = country_category
  )
) +
  geom_col(key_glyph = draw_key_country_fill) +
  geom_vline(
    xintercept = month_boundaries,
    linetype = "dashed",
    colour = "black"
  ) +
  geom_point(
    data = october_tops,
    aes(x = month_id, y = y_pos + 3000),
    inherit.aes = FALSE,
    shape = 6,
    size = 2.5,
    fill = "black"
  ) +
  scale_fill_manual(
    values = fatality_stack_colours,
    breaks = c(names(fatality_highlight_colours), "Others"),
    labels = c(names(fatality_highlight_colours), "Others"),
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    name = "Fatalities",
    breaks = seq(0, 150000, length.out = 6),
    labels = function(x) {
      format(
        x,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE,
        trim = TRUE
      )
    }
  ) +
  scale_x_continuous(
    name = "Month",
    breaks = country_breaks,
    labels = month_lookup_vec[as.character(country_breaks)]
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank(),
    legend.key.spacing.y = grid::unit(publication_legend_key_spacing_y_pt, "pt"),
    axis.ticks = element_line(colour = "black")
  ) +
  theme_fontsize

save_plot_for_latex(
  filename = file.path(publication_plot_dir, "Plot_VIEWS_fatalities.png"),
  plot_object = fatalities_worldwide_plot + theme_no_plot_title,
  height_px = single_panel_height_px,
  height_scale = publication_height_scale * (5 / 4)
)


# ==============================================================================
# Figures 2 and 3: country contributions
# ==============================================================================

# First average scores over the eight selected models for each country-month.
# This reproduces the aggregation used in the original country figures.
country_scores_by_month <- score_long %>%
  group_by(month_id, country_id) %>%
  summarise(
    crps = mean(crps),
    twcrps = mean(twcrps),
    brier = mean(brier_onset),
    .groups = "drop"
  )

prepare_country_stack <- function(
    score_data,
    score_column,
    highlighted_country_ids,
    highlighted_country_colours
) {
  highlighted_country_names <- names(highlighted_country_colours)
  
  highlighted_lookup <- country_id_list %>%
    filter(country_id %in% highlighted_country_ids)
  
  if (!setequal(highlighted_lookup$name, highlighted_country_names)) {
    stop(
      "Highlighted country IDs and highlighted colour names do not match for ",
      score_column,
      "."
    )
  }
  
  # Use exactly the same global CRPS ordering in every panel. Highlighted
  # countries are removed, while the relative order of all remaining countries
  # is preserved: largest overall CRPS first, then descending overall CRPS.
  other_country_order <- global_country_crps_order %>%
    filter(!country_id %in% highlighted_country_ids)
  
  if (anyNA(other_country_order$name)) {
    stop(
      "Missing country names among countries grouped as Others for ",
      score_column,
      "."
    )
  }
  
  other_keys <- paste0(
    "Others__",
    other_country_order$country_id
  )
  
  other_colours <- setNames(
    rep(
      c("grey30", "grey65"),
      length.out = length(other_keys)
    ),
    other_keys
  )
  
  # Keep one dummy "Others" level for a single legend entry. The actual bars
  # use one level per country in Others, alternating dark and light grey.
  stack_colours <- c(
    highlighted_country_colours,
    "Others" = "grey30",
    other_colours
  )
  
  stack_levels <- c(
    highlighted_country_names,
    other_keys,
    "Others"
  )
  
  plot_data <- score_data %>%
    left_join(country_id_list, by = "country_id") %>%
    mutate(
      country_category = if_else(
        country_id %in% highlighted_country_ids,
        name,
        paste0("Others__", country_id)
      ),
      country_category = factor(
        country_category,
        levels = stack_levels
      ),
      score_value = .data[[score_column]]
    ) %>%
    group_by(month_id, country_category) %>%
    summarise(
      value = sum(score_value),
      .groups = "drop"
    ) %>%
    bind_rows(
      tibble(
        month_id = min(score_data$month_id),
        country_category = factor("Others", levels = stack_levels),
        value = 0
      )
    )
  
  list(
    data = plot_data,
    colours = stack_colours,
    legend_breaks = c(highlighted_country_names, "Others"),
    legend_labels = c(highlighted_country_names, "Others")
  )
}

make_country_panel <- function(
    country_stack,
    y_label,
    crps_axis = FALSE
) {
  p <- ggplot(
    country_stack$data,
    aes(
      x = month_id,
      y = value,
      fill = country_category
    )
  ) +
    geom_col(key_glyph = draw_key_country_fill) +
    geom_vline(
      xintercept = month_boundaries,
      linetype = "dashed",
      colour = "black"
    ) +
    scale_fill_manual(
      values = country_stack$colours,
      breaks = country_stack$legend_breaks,
      labels = country_stack$legend_labels,
      drop = FALSE,
      name = NULL
    ) +
    scale_x_continuous(
      breaks = country_breaks,
      labels = month_lookup_vec[as.character(country_breaks)]
    ) +
    labs(
      x = "Month",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.key.spacing.y = grid::unit(publication_legend_key_spacing_y_pt, "pt"),
      axis.ticks = element_line(colour = "black")
    ) +
    theme_fontsize +
    theme_no_plot_title
  
  if (crps_axis) {
    p <- p +
      scale_y_continuous(
        breaks = seq(0, 150000, length.out = 6),
        labels = function(x) {
          format(
            x,
            big.mark = ".",
            decimal.mark = ",",
            scientific = FALSE,
            trim = TRUE
          )
        }
      )
  }
  
  p
}
# ---- Figure 2: same five countries as the fatalities figure ------------------

country_crps_fixed <- prepare_country_stack(
  country_scores_by_month,
  "crps",
  top_fatality_country_ids,
  fatality_highlight_colours
)

country_twcrps_fixed <- prepare_country_stack(
  country_scores_by_month,
  "twcrps",
  top_fatality_country_ids,
  fatality_highlight_colours
)

country_brier_fixed <- prepare_country_stack(
  country_scores_by_month,
  "brier",
  top_fatality_country_ids,
  fatality_highlight_colours
)

country_crps_plot <- make_country_panel(
  country_crps_fixed,
  "Mean CRPS for magnitude",
  crps_axis = TRUE
)

country_twcrps_plot <- make_country_panel(
  country_twcrps_fixed,
  "Mean CRPS for relative changes"
)

country_brier_plot <- make_country_panel(
  country_brier_fixed,
  mean_brier_score_label
)
views_triple_country_plot <-
  (country_crps_plot /
     country_twcrps_plot /
     country_brier_plot) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  ) &
  theme(legend.position = "right") &
  theme_no_plot_title

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0("Plot_VIEWS_triple_country_", greq_id, ".png")
  ),
  plot_object = views_triple_country_plot,
  height_px = triple_country_height_px,
  height_scale = publication_height_scale * (5 / 4)
)


# ---- Figure 3: same five coloured countries and ordering as fatalities --------
#
# Use exactly the same highlighted countries, colours and legend order in all
# three panels:
#   Ethiopia, Ukraine, Afghanistan, Yemen, Syria, Others.
#
# This keeps the country encoding identical across fatalities, CRPS, twCRPS
# and Brier score.

prepare_consistent_country_panel <- function(
    score_data,
    score_column,
    y_label,
    crps_axis = FALSE
) {
  country_stack <- prepare_country_stack(
    score_data,
    score_column,
    top_fatality_country_ids,
    fatality_highlight_colours
  )
  
  make_country_panel(
    country_stack,
    y_label,
    crps_axis = crps_axis
  )
}

country_crps_consistent_plot <- prepare_consistent_country_panel(
  country_scores_by_month,
  "crps",
  "Mean CRPS for magnitude",
  crps_axis = TRUE
)

country_twcrps_consistent_plot <- prepare_consistent_country_panel(
  country_scores_by_month,
  "twcrps",
  "Mean CRPS for relative changes"
)

country_brier_consistent_plot <- prepare_consistent_country_panel(
  country_scores_by_month,
  "brier",
  mean_brier_score_label
)

views_triple_country_metric_specific_plot <-
  (country_crps_consistent_plot /
     country_twcrps_consistent_plot /
     country_brier_consistent_plot) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  ) &
  theme(
    legend.position = "right",
    legend.title = element_blank()
  ) &
  theme_no_plot_title

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0(
      "Plot_VIEWS_triple_country_top5_per_measure_",
      greq_id,
      ".png"
    )
  ),
  plot_object = views_triple_country_metric_specific_plot,
  height_px = triple_country_height_px,
  height_scale = publication_height_scale * (5 / 4)
)


# ==============================================================================
# Figure 4: score contributions by conflict situation
# ==============================================================================

n_country_months <- nrow(observations_18_23)

situation_contributions <- score_long %>%
  group_by(model, situation_month) %>%
  summarise(
    CRPS = sum(crps) / n_country_months,
    twCRPS = sum(twcrps) / n_country_months,
    Brier = sum(brier_onset) / n_country_months,
    .groups = "drop"
  )

model_order_for_situation_plot <- situation_contributions %>%
  group_by(model) %>%
  summarise(
    total_CRPS = sum(CRPS),
    .groups = "drop"
  ) %>%
  arrange(desc(total_CRPS)) %>%
  pull(model)

situation_contributions <- situation_contributions %>%
  mutate(
    Model = factor(
      unname(model_labels[model]),
      levels = unname(model_labels[model_order_for_situation_plot])
    )
  )

situation_colours <- c(
  "conflict" = "#a22223",
  "deescalation" = "#d09191",
  "onset" = "#4664aa",
  "peace" = "#a2b2d4"
)

situation_labels <- c(
  "conflict" = "Sustained conflict",
  "deescalation" = "End/interruption of conflict",
  "onset" = "Conflict onset",
  "peace" = "Continued peace"
)

make_situation_panel <- function(score_column, x_label) {
  ggplot(
    situation_contributions,
    aes(
      y = Model,
      x = .data[[score_column]],
      fill = situation_month
    )
  ) +
    geom_col() +
    labs(
      x = x_label,
      y = NULL
    ) +
    scale_fill_manual(
      name = NULL,
      values = situation_colours,
      breaks = names(situation_colours),
      labels = situation_labels
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_blank()
    ) +
    theme_fontsize
}

crps_conflict_situation_plot <- make_situation_panel(
  "CRPS",
  "Mean CRPS for magnitude"
)

twcrps_conflict_situation_plot <- make_situation_panel(
  "twCRPS",
  "Mean CRPS for relative changes"
)

brier_conflict_situation_plot <- make_situation_panel(
  "Brier",
  mean_brier_score_label
)

views_triple_situation_plot <-
  (crps_conflict_situation_plot /
     twcrps_conflict_situation_plot /
     brier_conflict_situation_plot) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  theme_no_plot_title

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0("Plot_VIEWS_triple_situation_", greq_id, ".png")
  ),
  plot_object = views_triple_situation_plot,
  height_px = triple_situation_height_px
)


# ==============================================================================
# Figure 5: mean Brier score by threshold
# ==============================================================================

all_models_brier_full <- load_rdata_object(
  file.path(cache_dir, "all_models_brier_threshold_extended.RData"),
  "all_models_brier_full"
)

require_columns(
  all_models_brier_full,
  c("model", "a", "brier_avg"),
  "all_models_brier_full"
)

threshold_plot_data <- all_models_brier_full %>%
  filter(
    model %in% selected_models,
    a >= 0,
    a <= 200
  ) %>%
  mutate(log_a_plus1 = log1p(a)) %>%
  arrange(model, a)

for (model_name in selected_models) {
  model_thresholds <- threshold_plot_data %>%
    filter(model == model_name) %>%
    pull(a)
  
  if (!all(0:200 %in% model_thresholds)) {
    stop(
      "Threshold cache is incomplete on 0:200 for model: ",
      model_name
    )
  }
}

threshold_ticks <- c(0, 1, 5, 10, 25, 50, 100, 200)
threshold_y_limits <- c(
  0,
  1.02 * max(threshold_plot_data$brier_avg)
)

crps_brier_integrands_plot <- ggplot(
  threshold_plot_data,
  aes(
    x = a,
    y = brier_avg,
    colour = model,
    group = model
  )
) +
  geom_step(
    direction = "hv",
    linewidth = 0.95,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = fatality_threshold,
    linewidth = 0.6,
    linetype = "dashed",
    colour = "grey20"
  ) +
  labs(
    x = "Threshold a",
    y = "Mean Brier Score",
    colour = NULL
  ) +
  scale_x_continuous(
    limits = c(0, 200),
    breaks = threshold_ticks,
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = threshold_y_limits,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = selected_model_colors,
    breaks = selected_models,
    labels = selected_model_labels
  ) +
  theme_bw() +
  theme_fontsize +
  theme_no_plot_title

twcrps_brier_integrands_plot <- ggplot(
  threshold_plot_data,
  aes(
    x = log_a_plus1,
    y = brier_avg,
    colour = model,
    group = model
  )
) +
  geom_step(
    direction = "hv",
    linewidth = 0.95,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = log1p(fatality_threshold),
    linewidth = 0.6,
    linetype = "dashed",
    colour = "grey20"
  ) +
  labs(
    x = "Threshold a (log(a+1) spacing)",
    y = "Mean Brier Score",
    colour = NULL
  ) +
  scale_x_continuous(
    limits = log1p(c(0, 200)),
    breaks = log1p(threshold_ticks),
    labels = threshold_ticks,
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = threshold_y_limits,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = selected_model_colors,
    breaks = selected_models,
    labels = selected_model_labels
  ) +
  theme_bw() +
  theme_fontsize +
  theme_no_plot_title

views_double_threshold_plot <-
  (crps_brier_integrands_plot /
     twcrps_brier_integrands_plot) +
  plot_layout(
    guides = "collect",
    axes = "collect_y",
    axis_titles = "collect_y"
  ) &
  theme(legend.position = "right") &
  theme_no_plot_title

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0("Plot_VIEWS_double_thrBS_", greq_id, ".png")
  ),
  plot_object = views_double_threshold_plot,
  height_px = double_threshold_height_px
)

rm(all_models_brier_full, threshold_plot_data)


# ==============================================================================
# Figure 5: mean Brier score by threshold, a = 200 to 125000
# ==============================================================================

all_models_brier_full_ge200 <- load_rdata_object(
  file.path(cache_dir, "all_models_brier_threshold_extended.RData"),
  "all_models_brier_full"
)

require_columns(
  all_models_brier_full_ge200,
  c("model", "a", "brier_avg"),
  "all_models_brier_full_ge200"
)

threshold_plot_data_ge200 <- all_models_brier_full_ge200 %>%
  filter(
    model %in% selected_models,
    a >= 200,
    a <= 125000
  ) %>%
  mutate(log_a_plus1 = log1p(a)) %>%
  arrange(model, a)

for (model_name in selected_models) {
  model_thresholds <- threshold_plot_data_ge200 %>%
    filter(model == model_name) %>%
    pull(a)
  
  # if (!all(seq.int(200L, 125000L) %in% model_thresholds)) {
  #   stop(
  #     "Threshold cache is incomplete on 200:130000 for model: ",
  #     model_name
  #   )
  # }
}

threshold_ticks_ge200 <- c(
  200, 25000, 50000, 75000, 100000, 125000
)

threshold_log_ticks_ge200 <- c(
  200, 500, 1000, 5000, 10000, 50000, 125000
)

threshold_y_limits_ge200 <- c(
  0,
  1.02 * max(threshold_plot_data_ge200$brier_avg, na.rm = TRUE)
)

crps_brier_integrands_plot_ge200 <- ggplot(
  threshold_plot_data_ge200,
  aes(
    x = a,
    y = brier_avg,
    colour = model,
    group = model
  )
) +
  geom_step(
    direction = "hv",
    linewidth = 0.95,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = fatality_threshold,
    linewidth = 0.6,
    linetype = "dashed",
    colour = "grey20"
  ) +
  labs(
    x = "Threshold a",
    y = "Mean Brier Score",
    colour = NULL
  ) +
  scale_x_continuous(
    limits = c(200, 127000),
    breaks = threshold_ticks_ge200,
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = threshold_y_limits_ge200,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = selected_model_colors,
    breaks = selected_models,
    labels = selected_model_labels
  ) +
  theme_bw() +
  theme_fontsize +
  theme_no_plot_title

twcrps_brier_integrands_plot_ge200 <- ggplot(
  threshold_plot_data_ge200,
  aes(
    x = log_a_plus1,
    y = brier_avg,
    colour = model,
    group = model
  )
) +
  geom_step(
    direction = "hv",
    linewidth = 0.95,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = log1p(fatality_threshold),
    linewidth = 0.6,
    linetype = "dashed",
    colour = "grey20"
  ) +
  labs(
    x = "Threshold a (log(a+1) spacing)",
    y = "Mean Brier Score",
    colour = NULL
  ) +
  scale_x_continuous(
    limits = log1p(c(200, 138000)),
    breaks = log1p(threshold_log_ticks_ge200),
    labels = scales::label_comma()(threshold_log_ticks_ge200),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = threshold_y_limits_ge200,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(
    values = selected_model_colors,
    breaks = selected_models,
    labels = selected_model_labels
  ) +
  theme_bw() +
  theme_fontsize +
  theme_no_plot_title

views_double_threshold_plot_ge200 <-
  (crps_brier_integrands_plot_ge200 /
     twcrps_brier_integrands_plot_ge200) +
  plot_layout(
    guides = "collect",
    axes = "collect_y",
    axis_titles = "collect_y"
  ) &
  theme(legend.position = "right") &
  theme_no_plot_title

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0("Plot_VIEWS_double_thrBS_ge200_", greq_id, ".png")
  ),
  plot_object = views_double_threshold_plot_ge200,
  height_px = double_threshold_height_px
)

rm(all_models_brier_full_ge200, threshold_plot_data_ge200)


# ==============================================================================
# Figure 6: ROC curves after previous non-conflict
# ==============================================================================

roc_data <- previous_no_conflict_data %>%
  select(model, onset_prob_pred, actual_binary)

score_list_selected <- lapply(
  selected_models,
  function(model_name) {
    roc_data %>%
      filter(model == model_name) %>%
      pull(onset_prob_pred)
  }
)

label_list_selected <- lapply(
  selected_models,
  function(model_name) {
    roc_data %>%
      filter(model == model_name) %>%
      pull(actual_binary)
  }
)

scores_joined_selected <- precrec::join_scores(score_list_selected)
labels_joined_selected <- precrec::join_labels(label_list_selected)

mm_selected_models <- precrec::mmdata(
  scores = scores_joined_selected,
  labels = labels_joined_selected,
  modnames = selected_models
)

mm_eval <- precrec::evalmod(mm_selected_models)

auc_selected_models <- precrec::auc(mm_eval) %>%
  filter(curvetypes == "ROC") %>%
  mutate(modnames = as.character(modnames))

if (!setequal(auc_selected_models$modnames, selected_models)) {
  stop("ROC AUC output does not contain exactly the selected models.")
}

auc_labels <- setNames(
  paste0(
    unname(model_labels[auc_selected_models$modnames]),
    ", AUC = ",
    round(auc_selected_models$aucs, 2)
  ),
  auc_selected_models$modnames
)

roc_curve_selected_models <- autoplot(
  mm_eval,
  curvetype = "ROC"
) +
  geom_line(
    linewidth = 1,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = selected_model_colors,
    breaks = selected_models,
    labels = auc_labels[selected_models]
  ) +
  labs(
    x = "FAR",
    y = "HR",
    color = ""
  ) +
  theme_fontsize +
  theme(
    plot.margin = margin(4, 4, 4, 4, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt")
  )

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0("Plot_VIEWS_ROC_", greq_id, ".png")
  ),
  plot_object = roc_curve_selected_models + theme_no_plot_title,
  height_px = roc_plot_height_px
)


# ==============================================================================
# Figure 7: reliability diagrams after previous non-conflict
# ==============================================================================

reliabilitydiag_custom <- function(
    fcst,
    obs,
    path_colour,
    region_level = 0.9,
    region_position = "diagonal",
    region_method = "resampling"
) {
  r_selected <- reliabilitydiag(
    x = fcst,
    y = obs,
    region.level = region_level,
    region.method = region_method,
    region.position = region_position
  )
  
  reliability_plot <- autoplot(
    r_selected,
    params_histogram = list(
      yscale = 0.2,
      drop = "all",
      colour = "black",
      fill = NA
    ))
  
  # Remove package-supplied textual score annotations (e.g. MCB/DSC/UNC).
  is_text_annotation <- vapply(
    reliability_plot$layers,
    function(layer) {
      inherits(layer$geom, "GeomText") ||
        inherits(layer$geom, "GeomLabel")
    },
    logical(1)
  )
  reliability_plot$layers <- reliability_plot$layers[!is_text_annotation]
  
  # The version used for the paper adds an auxiliary fifth segment layer that
  # is not part of the desired publication display. Remove it only if it is
  # actually present and is a segment, avoiding the brittle recycled logical
  # index used in the exploratory script.
  if (
    length(reliability_plot$layers) >= 5L &&
    inherits(reliability_plot$layers[[5L]]$geom, "GeomSegment")
  ) {
    reliability_plot$layers <- reliability_plot$layers[-5L]
  }
  
  data_estim <- triptych::estimates(
    triptych::reliability(
      data.frame(forecast = fcst),
      y = obs
    )
  ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(x)
  
  # Add horizontal pieces explicitly when adjacent x-ordered estimates have
  # the same CEP. Ordering by x is essential: geom_path must follow forecast
  # value, not the ordinate.
  horizontal_segments <- data.frame()
  
  if (nrow(data_estim) > 1L) {
    for (i in seq_len(nrow(data_estim) - 1L)) {
      if (data_estim$CEP[i] == data_estim$CEP[i + 1L]) {
        horizontal_segments <- rbind(
          horizontal_segments,
          data.frame(
            x = data_estim$x[i],
            x_end = data_estim$x[i + 1L],
            CEP = data_estim$CEP[i],
            CEP_end = data_estim$CEP[i + 1L]
          )
        )
      }
    }
  }
  
  p <- reliability_plot +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "none"
    ) +
    geom_path(
      data = data_estim,
      aes(x = .data$x, y = .data$CEP),
      linewidth = 0.9,
      colour = path_colour
    ) +
    geom_segment(
      aes(
        x = 0,
        y = 0,
        xend = 1,
        yend = 1
      ),
      inherit.aes = FALSE,
      colour = "grey50",
      linewidth = 0.4
    )
  
  if (nrow(horizontal_segments) > 0L) {
    p <- p +
      geom_segment(
        data = horizontal_segments,
        aes(
          x = .data$x,
          y = .data$CEP,
          xend = .data$x_end,
          yend = .data$CEP_end
        ),
        linewidth = 1.3,
        colour = path_colour
      )
  } else {
    p <- p +
      geom_point(
        data = data_estim,
        aes(x = .data$x, y = .data$CEP),
        colour = path_colour,
        shape = 19,
        size = 1.7
      ) +
      coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0, 1)
      )
  }
  
  p
}

# Confidence regions use resampling, so fix the RNG state for reproducibility.
set.seed(20260915)

reliability_plots <- setNames(
  lapply(
    selected_models,
    function(model_name) {
      reliability_df <- previous_no_conflict_data %>%
        filter(model == model_name)
      
      reliabilitydiag_custom(
        reliability_df$onset_prob_pred,
        reliability_df$actual_binary,
        selected_model_colors[[model_name]]
      ) +
        labs(
          title = unname(model_labels[[model_name]]),
          x = "Forecast value",
          y = "CEP"
        ) +
        theme_fontsize +
        theme(
          plot.title = element_text(
            size = publication_panel_title_size_pt,
            hjust = 0.5,
            margin = margin(b = 4, unit = "pt")
          ),
          plot.title.position = "panel",
          plot.subtitle = element_blank(),
          axis.title = element_text(
            size = publication_axis_title_size_pt
          ),
          axis.text = element_text(
            size = publication_base_size_pt
          ),
          plot.margin = margin(
            t = 3,
            r = 4,
            b = 1,
            l = 4,
            unit = "pt"
          )
        )
    }
  ),
  selected_models
)

reliability_grid <- patchwork::wrap_plots(
  unname(reliability_plots),
  ncol = 4,
  nrow = 2
) +
  plot_layout(heights = c(1, 1))

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0("Plot_VIEWS_reliability_", greq_id, ".png")
  ),
  plot_object = reliability_grid,
  height_px = reliability_plot_height_px,
  height_scale = 1
)


# ==============================================================================
# Figure 8: model ranking across three scores
# ==============================================================================

make_model_rankings <- function(score_data, situations = NULL) {
  ranking_data <- score_data
  
  if (!is.null(situations)) {
    ranking_data <- ranking_data %>%
      filter(situation_month %in% situations)
  }
  
  score_means <- ranking_data %>%
    group_by(model) %>%
    summarise(
      mean_CRPS = mean(crps),
      mean_twCRPS = mean(twcrps),
      mean_Brier = mean(brier_onset),
      .groups = "drop"
    )
  
  missing_models <- setdiff(selected_models, score_means$model)
  
  if (length(missing_models) > 0L) {
    stop(
      "Ranking data are missing selected model(s): ",
      paste(missing_models, collapse = ", ")
    )
  }
  
  if (any(!is.finite(
    as.matrix(
      score_means[, c("mean_CRPS", "mean_twCRPS", "mean_Brier")]
    )
  ))) {
    stop("At least one conditional mean score is not finite.")
  }
  
  crps_ranks <- score_means %>%
    arrange(mean_CRPS, match(model, selected_models)) %>%
    transmute(model, CRPS = row_number())
  
  twcrps_ranks <- score_means %>%
    arrange(mean_twCRPS, match(model, selected_models)) %>%
    transmute(model, twCRPS = row_number())
  
  brier_ranks <- score_means %>%
    arrange(mean_Brier, match(model, selected_models)) %>%
    transmute(model, BRIER = row_number())
  
  list(
    crps_ranks,
    twcrps_ranks,
    brier_ranks
  ) %>%
    purrr::reduce(left_join, by = "model") %>%
    mutate(Model = unname(model_labels[model])) %>%
    select(Model, CRPS, twCRPS, BRIER)
}

model_rankings_all <- make_model_rankings(score_long)

model_rankings_previous_no_conflict <- make_model_rankings(
  score_long,
  situations = c("peace", "onset")
)

model_rankings_previous_conflict <- make_model_rankings(
  score_long,
  situations = c("conflict", "deescalation")
)

smoothstep <- function(u) {
  3 * u^2 - 2 * u^3
}

plateau_step <- function(t, flat = 0.42) {
  flat <- pmin(pmax(flat, 0), 0.49)
  
  ifelse(
    t <= flat,
    0,
    ifelse(
      t >= 1 - flat,
      1,
      smoothstep((t - flat) / (1 - 2 * flat))
    )
  )
}

make_segment <- function(y0, y1, n = 160L, flat = 0.42) {
  x <- seq(0, 1, length.out = n)
  s <- plateau_step(x, flat = flat)
  
  data.frame(
    x = x,
    y = y0 + (y1 - y0) * s
  )
}

ranking_flat <- 0.33
ranking_tick_length <- 0.03
ranking_number_offset <- 0.055
ranking_axis_colour <- "grey35"
ranking_axis_linewidth <- 0.5
ranking_header_y <- -0.30

ranking_model_colours <- setNames(
  unname(selected_model_colors[selected_models]),
  unname(model_labels[selected_models])
)

brier_ranking_header <- paste0(
  "Brier Score for k=",
  fatality_threshold
)

ranking_base_theme <- theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.background = element_blank(),
    axis.line = element_blank()
  )

make_model_ranking_panel <- function(
    model_rankings,
    ranking_y_title
) {
  left_curves <- purrr::pmap_dfr(
    model_rankings[, c("CRPS", "twCRPS", "Model")],
    function(CRPS, twCRPS, Model) {
      segment <- make_segment(
        CRPS,
        twCRPS,
        flat = ranking_flat
      )
      segment$Model <- Model
      segment
    }
  )
  
  right_curves <- purrr::pmap_dfr(
    model_rankings[, c("twCRPS", "BRIER", "Model")],
    function(twCRPS, BRIER, Model) {
      segment <- make_segment(
        twCRPS,
        BRIER,
        flat = ranking_flat
      )
      segment$Model <- Model
      segment
    }
  )
  
  left_points <- model_rankings %>%
    transmute(
      Model,
      y = CRPS,
      y2 = twCRPS
    )
  
  right_points <- model_rankings %>%
    transmute(
      Model,
      y = twCRPS,
      y2 = BRIER
    )
  
  number_of_models <- nrow(model_rankings)
  ranking_ticks <- data.frame(y = seq_len(number_of_models))
  ranking_limits <- c(number_of_models + 0.4, -0.30)
  
  ranking_left <- ggplot() +
    geom_vline(
      xintercept = c(0, 1),
      colour = ranking_axis_colour,
      linewidth = ranking_axis_linewidth
    ) +
    geom_segment(
      data = ranking_ticks,
      aes(
        x = 0,
        xend = -ranking_tick_length,
        y = y,
        yend = y
      ),
      colour = ranking_axis_colour,
      linewidth = 0.4
    ) +
    geom_text(
      data = ranking_ticks,
      aes(
        x = -ranking_number_offset,
        y = y,
        label = y
      ),
      hjust = 1,
      colour = ranking_axis_colour,
      size = publication_geom_text_size
    ) +
    geom_segment(
      data = ranking_ticks,
      aes(
        x = 1,
        xend = 1 - ranking_tick_length,
        y = y,
        yend = y
      ),
      colour = ranking_axis_colour,
      linewidth = 0.35
    ) +
    geom_line(
      data = left_curves,
      aes(
        x = x,
        y = y,
        colour = Model
      ),
      linewidth = 1.1
    ) +
    geom_point(
      data = left_points,
      aes(
        x = 0,
        y = y,
        colour = Model
      ),
      size = 2
    ) +
    geom_point(
      data = left_points,
      aes(
        x = 1,
        y = y2,
        colour = Model
      ),
      size = 2
    ) +
    geom_text(
      data = model_rankings,
      aes(
        x = -0.15,
        y = CRPS,
        label = Model,
        colour = Model
      ),
      hjust = 1,
      size = publication_geom_text_size
    ) +
    annotate(
      "text",
      x = -0.012,
      y = ranking_header_y,
      label = "CRPS for magnitude",
      fontface = "bold",
      hjust = 1,
      colour = ranking_axis_colour,
      size = publication_geom_text_size
    ) +
    annotate(
      "text",
      x = 1 - 0.012,
      y = ranking_header_y,
      label = "CRPS for relative changes",
      fontface = "bold",
      hjust = 1,
      colour = ranking_axis_colour,
      size = publication_geom_text_size
    ) +
    scale_colour_manual(
      values = ranking_model_colours,
      guide = "none"
    ) +
    scale_x_continuous(
      name = NULL,
      limits = c(-0.15, 1),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = ranking_limits,
      breaks = NULL
    ) +
    coord_cartesian(clip = "off") +
    ranking_base_theme +
    theme(
      axis.title.y = element_text(
        size = publication_axis_title_size_pt,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        margin = margin(r = 145, unit = "pt")
      ),
      plot.margin = margin(6, 0, 6, 8, unit = "pt")
    ) +
    labs(y = ranking_y_title)
  
  ranking_right <- ggplot() +
    labs(y = NULL) +
    geom_vline(
      xintercept = c(0, 1),
      colour = ranking_axis_colour,
      linewidth = ranking_axis_linewidth
    ) +
    geom_segment(
      data = ranking_ticks,
      aes(
        x = 1,
        xend = 1 + ranking_tick_length,
        y = y,
        yend = y
      ),
      colour = ranking_axis_colour,
      linewidth = 0.4
    ) +
    geom_text(
      data = ranking_ticks,
      aes(
        x = 1 + ranking_number_offset,
        y = y,
        label = y
      ),
      hjust = 0,
      colour = ranking_axis_colour,
      size = publication_geom_text_size
    ) +
    geom_line(
      data = right_curves,
      aes(
        x = x,
        y = y,
        colour = Model
      ),
      linewidth = 1.1
    ) +
    geom_point(
      data = right_points,
      aes(
        x = 0,
        y = y,
        colour = Model
      ),
      size = 2
    ) +
    geom_point(
      data = right_points,
      aes(
        x = 1,
        y = y2,
        colour = Model
      ),
      size = 2
    ) +
    annotate(
      "text",
      x = 1 - 0.012,
      y = ranking_header_y,
      label = brier_ranking_header,
      fontface = "bold",
      hjust = 1,
      colour = ranking_axis_colour,
      size = publication_geom_text_size
    ) +
    scale_colour_manual(
      values = ranking_model_colours,
      guide = "none"
    ) +
    scale_x_continuous(
      name = NULL,
      limits = c(0, 1.1),
      breaks = NULL,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      limits = ranking_limits,
      breaks = NULL
    ) +
    ranking_base_theme +
    theme(
      plot.margin = margin(6, 80, 6, -90, unit = "pt")
    )
  
  (ranking_left | ranking_right) +
    plot_layout(widths = c(1.05, 1))
}

ranking_panel_all <- make_model_ranking_panel(
  model_rankings_all,
  "All situations"
)

ranking_panel_previous_no_conflict <- make_model_ranking_panel(
  model_rankings_previous_no_conflict,
  "No previous conflict"
)

ranking_panel_previous_conflict <- make_model_ranking_panel(
  model_rankings_previous_conflict,
  "Previous conflict"
)

views_triple_model_ranking_plot <- patchwork::wrap_plots(
  list(
    ranking_panel_all,
    ranking_panel_previous_no_conflict,
    ranking_panel_previous_conflict
  ),
  ncol = 1,
  heights = c(1, 1, 1)
)

save_plot_for_latex(
  filename = file.path(
    publication_plot_dir,
    paste0(
      "Plot_VIEWS_triple_model_ranking_",
      greq_id,
      ".png"
    )
  ),
  plot_object = views_triple_model_ranking_plot & theme_no_plot_title,
  height_px = triple_model_ranking_height_px
)


# ---- Completion message ------------------------------------------------------

generated_figures <- c(
  "Plot_VIEWS_fatalities.png",
  paste0("Plot_VIEWS_triple_country_", greq_id, ".png"),
  paste0(
    "Plot_VIEWS_triple_country_top5_per_measure_",
    greq_id,
    ".png"
  ),
  paste0("Plot_VIEWS_triple_situation_", greq_id, ".png"),
  paste0("Plot_VIEWS_double_thrBS_", greq_id, ".png"),
  paste0("Plot_VIEWS_double_thrBS_ge200_", greq_id, ".png"),
  paste0("Plot_VIEWS_ROC_", greq_id, ".png"),
  paste0("Plot_VIEWS_reliability_", greq_id, ".png"),
  paste0("Plot_VIEWS_triple_model_ranking_", greq_id, ".png")
)

message(
  "Generated publication figures:\n",
  paste(file.path(publication_plot_dir, generated_figures), collapse = "\n")
)



##### Additional analyses of contributions to mean CRPS for magnitude per conflict situations

# ==============================================================================
# Sustained conflict events ranked by their mean percentage of the dark red CRPS segments
# ==============================================================================

# An event is one country-target-month classified as "conflict" by the script.
sustained_conflict_scores <- score_long %>%
  filter(model %in% selected_models, situation_month == "conflict")

if (nrow(sustained_conflict_scores) == 0L) {
  stop("No country-target-months are classified as sustained conflict.")
}

assert_unique_key(
  sustained_conflict_scores,
  c("model", "country_id", "month_id"),
  "sustained_conflict_scores"
)

sustained_conflict_model_totals <- sustained_conflict_scores %>%
  group_by(model) %>%
  summarise(sustained_conflict_crps_sum = sum(crps), .groups = "drop")

if (any(!is.finite(sustained_conflict_model_totals$sustained_conflict_crps_sum) |
        sustained_conflict_model_totals$sustained_conflict_crps_sum <= 0)) {
  stop("Each model must have a finite, positive sustained conflict CRPS total.")
}

sustained_conflict_scores <- sustained_conflict_scores %>%
  left_join(sustained_conflict_model_totals, by = "model") %>%
  mutate(
    # The common country-month denominator cancels in this percentage.
    pct_of_dark_red = 100 * crps / sustained_conflict_crps_sum,
    CRPS_component = crps / n_country_months
  )

sustained_conflict_events <- sustained_conflict_scores %>%
  group_by(country_id, month_id) %>%
  summarise(
    n_models = n_distinct(model),
    avg_dark_red_share_pct = mean(pct_of_dark_red),
    avg_CRPS_component = mean(CRPS_component),
    .groups = "drop"
  )

if (any(sustained_conflict_events$n_models != length(selected_models))) {
  stop("Every sustained conflict event must have a score for every selected model.")
}

sustained_conflict_events <- sustained_conflict_events %>%
  left_join(country_id_list, by = "country_id") %>%
  left_join(
    observations_18_23 %>%
      transmute(country_id, month_id, target_fatalities = outcome),
    by = c("country_id", "month_id")
  ) %>%
  mutate(
    # Exactly the reference-October rule used in the situation plot.
    reference_october_id = first_reference_october_id +
      12L * ((month_id - first_forecast_month_id) %/% 12L)
  )

# One percentage column per model, in selected_models order.
sustained_conflict_model_percentages <- purrr::map(
  selected_models,
  function(model_name) {
    result <- sustained_conflict_scores %>%
      filter(model == model_name) %>%
      select(country_id, month_id, pct_of_dark_red)
    names(result)[3L] <- paste0("pct_", unname(model_labels[[model_name]]))
    result
  }
) %>%
  purrr::reduce(left_join, by = c("country_id", "month_id"))

# Retrieve October and its preceding 11 months for each country/reference year.
sustained_conflict_reference_keys <- sustained_conflict_events %>%
  distinct(country_id, reference_october_id)

sustained_conflict_history_keys <- purrr::map_dfr(
  seq.int(-11L, 0L),
  function(relative_month) {
    sustained_conflict_reference_keys %>%
      transmute(country_id, month_id = reference_october_id + relative_month)
  }
) %>%
  distinct(country_id, month_id)

# Reuse the plot's historical outcomes. Supplement earlier months from features:
# reference October 2017 requires November-December 2016, which the original
# historical_outcomes object does not contain.
sustained_conflict_history <- historical_outcomes %>%
  select(country_id, month_id, fatalities = outcome)

sustained_conflict_missing_history_keys <- sustained_conflict_history_keys %>%
  anti_join(sustained_conflict_history, by = c("country_id", "month_id"))

if (nrow(sustained_conflict_missing_history_keys) > 0L) {
  sustained_conflict_extra_history <- arrow::read_parquet(
    features_path,
    col_select = c("country_id", "month_id", "ged_sb")
  ) %>%
    transmute(country_id, month_id, fatalities = ged_sb) %>%
    semi_join(sustained_conflict_missing_history_keys, by = c("country_id", "month_id"))
  
  sustained_conflict_history <- bind_rows(sustained_conflict_history, sustained_conflict_extra_history)
}

assert_unique_key(
  sustained_conflict_history,
  c("country_id", "month_id"),
  "sustained_conflict_history"
)

# Columns are chronological: -11 = previous November, ..., -1 = September.
# Include reference October itself as an additional check on the classification.
sustained_conflict_history_columns <- c(
  sprintf("fatalities_Oct_minus_%02d", 11:1),
  "fatalities_reference_Oct"
)

sustained_conflict_history_wide <- purrr::map(
  seq.int(-11L, 0L),
  function(relative_month) {
    result <- sustained_conflict_reference_keys %>%
      mutate(month_id = reference_october_id + relative_month) %>%
      left_join(sustained_conflict_history, by = c("country_id", "month_id")) %>%
      select(country_id, reference_october_id, fatalities)
    
    names(result)[3L] <- if (relative_month == 0L) {
      "fatalities_reference_Oct"
    } else {
      sprintf("fatalities_Oct_minus_%02d", -relative_month)
    }
    result
  }
) %>%
  purrr::reduce(left_join, by = c("country_id", "reference_october_id"))

if (anyNA(sustained_conflict_history_wide[sustained_conflict_history_columns])) {
  warning("Some historical fatalities are unavailable; retained as NA, not zero.")
}

# Extend the existing calendar labels to earlier monthly IDs.
sustained_conflict_month_label <- function(month_id) {
  origin_date <- min(month_calendar$date)
  month_index <- 12L * as.integer(format(origin_date, "%Y")) +
    as.integer(format(origin_date, "%m")) - 1L +
    (month_id - min(month_calendar$month_id))
  sprintf("%04d-%02d", month_index %/% 12L, month_index %% 12L + 1L)
}

sustained_conflict_contribution_table <- sustained_conflict_events %>%
  left_join(sustained_conflict_model_percentages, by = c("country_id", "month_id")) %>%
  left_join(
    sustained_conflict_history_wide,
    by = c("country_id", "reference_october_id")
  ) %>%
  arrange(desc(avg_dark_red_share_pct), country_id, month_id) %>%
  mutate(
    rank = row_number(),
    cumulative_dark_red_share_pct = cumsum(avg_dark_red_share_pct),
    target_month = sustained_conflict_month_label(month_id),
    reference_october = sustained_conflict_month_label(reference_october_id)
  ) %>%
  select(
    rank,
    country_id,
    country = name,
    target_month,
    target_fatalities,
    reference_october,
    avg_dark_red_share_pct,
    cumulative_dark_red_share_pct,
    avg_CRPS_component,
    all_of(paste0("pct_", unname(model_labels[selected_models]))),
    all_of(sustained_conflict_history_columns)
  )

readr::write_csv(
  sustained_conflict_contribution_table,
  file.path(
    publication_plot_dir,
    paste0("Table_VIEWS_sustained_conflict_CRPS_contributions_", greq_id, ".csv")
  ),
  na = "NA"
)

print(sustained_conflict_contribution_table, n = Inf, width = Inf)

head(sustained_conflict_contribution_table, 20)

# ==============================================================================
# Ethiopia's contribution to the dark blue magnitude-CRPS components
# ==============================================================================

ethiopia_country_id <- country_id_list %>%
  filter(trimws(name) == "Ethiopia") %>%
  pull(country_id)

if (length(ethiopia_country_id) != 1L) {
  stop("Expected exactly one country_id for Ethiopia.")
}

ethiopia_onset_crps_contribution <- score_long %>%
  group_by(model) %>%
  summarise(
    # Full dark blue segment: all countries' onset CRPS / all country-months.
    dark_blue_CRPS_component =
      sum(crps[situation_month == "onset"]) / n_country_months,
    
    # Ethiopia's contribution, using the same denominator.
    ethiopia_CRPS_component = sum(
      crps[
        situation_month == "onset" &
          country_id %in% ethiopia_country_id
      ]
    ) / n_country_months,
    
    .groups = "drop"
  ) %>%
  mutate(
    Model = unname(model_labels[model]),
    
    # Ethiopia's percentage of the dark blue segment.
    ethiopia_share_pct = if_else(
      dark_blue_CRPS_component > 0,
      100 * ethiopia_CRPS_component / dark_blue_CRPS_component,
      NA_real_
    )
  ) %>%
  arrange(match(model, model_order_for_situation_plot)) %>%
  select(
    Model,
    dark_blue_CRPS_component,
    ethiopia_CRPS_component,
    ethiopia_share_pct
  )

print(
  ethiopia_onset_crps_contribution %>%
    mutate(ethiopia_share_pct = round(ethiopia_share_pct, 2)),
  n = Inf,
  width = Inf
)

# ==============================================================================
# Onset events ranked by their mean percentage of the dark blue CRPS segments
# ==============================================================================

# An event is one country-target-month classified as "onset" by the script.
onset_scores <- score_long %>%
  filter(model %in% selected_models, situation_month == "onset")

if (nrow(onset_scores) == 0L) {
  stop("No country-target-months are classified as conflict onsets.")
}

assert_unique_key(
  onset_scores,
  c("model", "country_id", "month_id"),
  "onset_scores"
)

onset_model_totals <- onset_scores %>%
  group_by(model) %>%
  summarise(onset_crps_sum = sum(crps), .groups = "drop")

if (any(!is.finite(onset_model_totals$onset_crps_sum) |
        onset_model_totals$onset_crps_sum <= 0)) {
  stop("Each model must have a finite, positive onset CRPS total.")
}

onset_scores <- onset_scores %>%
  left_join(onset_model_totals, by = "model") %>%
  mutate(
    # The common country-month denominator cancels in this percentage.
    pct_of_dark_blue = 100 * crps / onset_crps_sum,
    CRPS_component = crps / n_country_months
  )

onset_events <- onset_scores %>%
  group_by(country_id, month_id) %>%
  summarise(
    n_models = n_distinct(model),
    avg_dark_blue_share_pct = mean(pct_of_dark_blue),
    avg_CRPS_component = mean(CRPS_component),
    .groups = "drop"
  )

if (any(onset_events$n_models != length(selected_models))) {
  stop("Every onset must have a score for every selected model.")
}

onset_events <- onset_events %>%
  left_join(country_id_list, by = "country_id") %>%
  left_join(
    observations_18_23 %>%
      transmute(country_id, month_id, target_fatalities = outcome),
    by = c("country_id", "month_id")
  ) %>%
  mutate(
    # Exactly the reference-October rule used in the situation plot.
    reference_october_id = first_reference_october_id +
      12L * ((month_id - first_forecast_month_id) %/% 12L)
  )

# One percentage column per model, in selected_models order.
onset_model_percentages <- purrr::map(
  selected_models,
  function(model_name) {
    result <- onset_scores %>%
      filter(model == model_name) %>%
      select(country_id, month_id, pct_of_dark_blue)
    names(result)[3L] <- paste0("pct_", unname(model_labels[[model_name]]))
    result
  }
) %>%
  purrr::reduce(left_join, by = c("country_id", "month_id"))

# Retrieve October and its preceding 11 months for each country/reference year.
onset_reference_keys <- onset_events %>%
  distinct(country_id, reference_october_id)

onset_history_keys <- purrr::map_dfr(
  seq.int(-11L, 0L),
  function(relative_month) {
    onset_reference_keys %>%
      transmute(country_id, month_id = reference_october_id + relative_month)
  }
) %>%
  distinct(country_id, month_id)

# Reuse the plot's historical outcomes. Supplement earlier months from features:
# reference October 2017 requires November-December 2016, which the original
# historical_outcomes object does not contain.
onset_history <- historical_outcomes %>%
  select(country_id, month_id, fatalities = outcome)

onset_missing_history_keys <- onset_history_keys %>%
  anti_join(onset_history, by = c("country_id", "month_id"))

if (nrow(onset_missing_history_keys) > 0L) {
  onset_extra_history <- arrow::read_parquet(
    features_path,
    col_select = c("country_id", "month_id", "ged_sb")
  ) %>%
    transmute(country_id, month_id, fatalities = ged_sb) %>%
    semi_join(onset_missing_history_keys, by = c("country_id", "month_id"))
  
  onset_history <- bind_rows(onset_history, onset_extra_history)
}

assert_unique_key(
  onset_history,
  c("country_id", "month_id"),
  "onset_history"
)

# Columns are chronological: -11 = previous November, ..., -1 = September.
# Include reference October itself as an additional check on the classification.
onset_history_columns <- c(
  sprintf("fatalities_Oct_minus_%02d", 11:1),
  "fatalities_reference_Oct"
)

onset_history_wide <- purrr::map(
  seq.int(-11L, 0L),
  function(relative_month) {
    result <- onset_reference_keys %>%
      mutate(month_id = reference_october_id + relative_month) %>%
      left_join(onset_history, by = c("country_id", "month_id")) %>%
      select(country_id, reference_october_id, fatalities)
    
    names(result)[3L] <- if (relative_month == 0L) {
      "fatalities_reference_Oct"
    } else {
      sprintf("fatalities_Oct_minus_%02d", -relative_month)
    }
    result
  }
) %>%
  purrr::reduce(left_join, by = c("country_id", "reference_october_id"))

if (anyNA(onset_history_wide[onset_history_columns])) {
  warning("Some historical fatalities are unavailable; retained as NA, not zero.")
}

# Extend the existing calendar labels to earlier monthly IDs.
onset_month_label <- function(month_id) {
  origin_date <- min(month_calendar$date)
  month_index <- 12L * as.integer(format(origin_date, "%Y")) +# Count onset-labelled events with >= 25 fatalities in at least one
    # of the 11 months preceding reference October (excluding October).
    
    pre_october_columns <- sprintf("fatalities_Oct_minus_%02d", 11:1)
  
  onset_contribution_table <- onset_contribution_table %>%
    mutate(
      conflict_before_reference_october = if_any(
        all_of(pre_october_columns),
        ~ .x >= 25
      )
    )
  
  onset_prior_conflict_summary <- onset_contribution_table %>%
    summarise(
      total_onset_events = n(),
      n_with_prior_conflict = sum(
        conflict_before_reference_october,
        na.rm = TRUE
      ),
      n_without_prior_conflict = sum(
        !conflict_before_reference_october,
        na.rm = TRUE
      ),
      n_unknown = sum(is.na(conflict_before_reference_october)),
      percentage_with_prior_conflict =
        100 * n_with_prior_conflict / total_onset_events
    )
  
  print(onset_prior_conflict_summary, width = Inf)
  as.integer(format(origin_date, "%m")) - 1L +
    (month_id - min(month_calendar$month_id))
  sprintf("%04d-%02d", month_index %/% 12L, month_index %% 12L + 1L)
}

onset_contribution_table <- onset_events %>%
  left_join(onset_model_percentages, by = c("country_id", "month_id")) %>%
  left_join(
    onset_history_wide,
    by = c("country_id", "reference_october_id")
  ) %>%
  arrange(desc(avg_dark_blue_share_pct), country_id, month_id) %>%
  mutate(
    rank = row_number(),
    target_month = onset_month_label(month_id),
    reference_october = onset_month_label(reference_october_id)
  ) %>%
  select(
    rank,
    country_id,
    country = name,
    target_month,
    target_fatalities,
    reference_october,
    avg_dark_blue_share_pct,
    avg_CRPS_component,
    all_of(paste0("pct_", unname(model_labels[selected_models]))),
    all_of(onset_history_columns)
  )

onset_contribution_table <- onset_contribution_table %>%
  mutate(
    cumulative_dark_blue_share_pct = cumsum(avg_dark_blue_share_pct),
    .after = avg_dark_blue_share_pct
  )


onset_contribution_table %>% arrange(country) %>% filter(country == "Ethiopia")

pre_october_columns <- sprintf("fatalities_Oct_minus_%02d", 11:1)

onset_contribution_table <- onset_contribution_table %>%
  mutate(
    previous_conflict_months = rowSums(
      across(all_of(pre_october_columns), ~ .x >= 25)
    )
  )

onset_previous_conflict_distribution <- tibble::tibble(
  previous_conflict_months = 1:11
) %>%
  left_join(
    onset_contribution_table %>%
      filter(!is.na(previous_conflict_months)) %>%
      group_by(previous_conflict_months) %>%
      summarise(
        n_countries = n_distinct(country_id),
        n_country_reference_octobers = n_distinct(
          country_id,
          reference_october
        ),
        n_onset_events = n(),
        .groups = "drop"
      ),
    by = "previous_conflict_months"
  ) %>%
  mutate(
    across(starts_with("n_"), ~ coalesce(.x, 0L))
  )

print(onset_previous_conflict_distribution, n = Inf, width = Inf)

cumsum(onset_previous_conflict_distribution)

readr::write_csv(
  onset_contribution_table,
  file.path(
    publication_plot_dir,
    paste0("Table_VIEWS_onset_CRPS_contributions_", greq_id, ".csv")
  ),
  na = "NA"
)

print(onset_contribution_table, n = Inf, width = Inf)

head(onset_contribution_table, 20)
