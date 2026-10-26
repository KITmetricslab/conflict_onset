# UCDP GED v25.1 | Country-Month Fatalities | 2018--2023
# Autor: ChatGPT | Stand: 2026-10-26
# Aufruf (Verzeichnis mit GEDEvent_v25_1.csv und 6 Parquets):
# Rscript ucdp_multimonth_analysis.R "/PFAD/ZU/DATEIEN"
#
# Benötigte Pakete (einmalig installieren):
# install.packages(c("data.table", "arrow", "ggplot2", "clue", "scales"))
#
# WICHTIG:
# - Die Parquets enthalten offenbar type_of_violence == 1 (state-based).
# - month_id 457 = 2018-01, 528 = 2023-12 (Nullpunkt 1980-01).
# - parquet$country_id ist NICHT das UCDP-country_id-Feld! Ein Länderabgleich
#   wird über den monatlichen Verlauf und eine optimale 1:1-Zuordnung geschätzt.
# - Gleichmäßige Verteilung über Monate ist nur eine SENSITIVITÄTSANALYSE;
#   daraus lässt sich die tatsächliche monatliche Opferzahl nicht ableiten.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(clue)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[[1L]] else "."
output_dir <- file.path(input_dir, "ucdp_results_R")
plot_dir <- file.path(output_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

START_YEAR <- 2018L
END_YEAR <- 2023L
EXCLUDE_LAST_MULTIMONTH <- TRUE   # "alle bis auf das allerletzte" (nach date_end + id)
EXPORT_ALL_ORIGINAL_COLUMNS <- TRUE

csv_file <- file.path(input_dir, "GEDEvent_v25_1.csv")
parquet_files <- file.path(input_dir, paste0("cm_actuals_", START_YEAR:END_YEAR, ".parquet"))
stopifnot(file.exists(csv_file), all(file.exists(parquet_files)))

cat("1) Original-UCDP laden ...\n")
cols <- c("id", "relid", "year", "country", "country_id", "region",
          "type_of_violence", "date_prec", "event_clarity", "date_start",
          "date_end", "best", "low", "high", "conflict_name", "dyad_name",
          "conflict_new_id", "dyad_new_id", "where_coordinates")
events <- fread(csv_file, select = cols, showProgress = TRUE)
events[, date_start := as.Date(substr(date_start, 1L, 10L))]
events[, date_end   := as.Date(substr(date_end,   1L, 10L))]
if (anyNA(events$date_start) || anyNA(events$date_end)) {
  stop("Datum konnte nicht vollständig eingelesen werden.")
}
if (any(events$date_end < events$date_start)) stop("date_end liegt vor date_start")
events <- events[year >= START_YEAR & year <= END_YEAR]
events[, `:=`(
  month_start = as.Date(format(date_start, "%Y-%m-01")),
  month_end   = as.Date(format(date_end,   "%Y-%m-01")),
  days_spanned = as.integer(date_end - date_start) + 1L
)]
events[, `:=`(
  multi_month = month_start != month_end,
  months_spanned = 12L * (
    as.integer(format(month_end, "%Y")) - as.integer(format(month_start, "%Y"))
  ) + as.integer(format(month_end, "%m")) - as.integer(format(month_start, "%m")) + 1L
)]
events[, end_month_id := (as.integer(format(month_end, "%Y")) - 1980L) * 12L +
                           as.integer(format(month_end, "%m"))]

all_multimonth <- copy(events[multi_month == TRUE])
state_events <- copy(events[type_of_violence == 1L])
state_multimonth <- copy(state_events[multi_month == TRUE])

cat("Gesamte state-based Events:", nrow(state_events), "\n")
cat("Mehrmonatige state-based Events (vor Ausschluss):", nrow(state_multimonth), "\n")
cat("Fatalities in diesen mehrmonatigen Events:", sum(state_multimonth$best), "\n")

# Ursprüngliche, nicht ausgefilterte Zahlen immer archivieren.
annual_before <- state_events[, .(
  events = .N,
  multi_events = sum(multi_month),
  total_fatalities = sum(best),
  multi_fatalities = sum(best[multi_month])
), by = year][order(year)]
annual_before[, `:=`(
  multi_event_share = multi_events / events,
  multi_fatality_share = multi_fatalities / total_fatalities
)]
fwrite(annual_before, file.path(output_dir, "yearly_summary_unfiltered.csv"))

# Die letzte mehrmonatige Beobachtung wird nach Enddatum UND Event-ID gewählt.
# Für den bereitgestellten Datensatz lautet diese ID 561999, best=0.
excluded_id <- NA_integer_
if (EXCLUDE_LAST_MULTIMONTH && nrow(all_multimonth) > 0L) {
  setorder(all_multimonth, date_end, id)
  excluded_id <- tail(all_multimonth$id, 1L)
  cat("Ausgeschlossen: letztes mehrmonatiges Event, id =", excluded_id, "\n")
  all_multimonth <- all_multimonth[id != excluded_id]
  state_multimonth <- state_multimonth[id != excluded_id]
  state_events <- state_events[id != excluded_id]
}

# Gezielter Event-Export: gut für manuelle Plausibilitätskontrolle.
setorder(all_multimonth, date_end, id)
setorder(state_multimonth, date_end, id)
fwrite(all_multimonth, file.path(output_dir, "multimonth_events_all_types.csv"))
fwrite(state_multimonth, file.path(output_dir, "multimonth_events_state_based.csv"))

if (EXPORT_ALL_ORIGINAL_COLUMNS) {
  # Liest zusätzlich ALLE Originalspalten einschließlich Quellenfeldern.
  # Bei geringem RAM kann dies mit FALSE deaktiviert werden.
  cat("2) Volle Originalspalten der betroffenen Events exportieren ...\n")
  full_ucdp <- fread(csv_file, showProgress = TRUE)
  fwrite(full_ucdp[id %in% all_multimonth$id],
         file.path(output_dir, "original_full_rows_all_multimonth.csv"))
  fwrite(full_ucdp[id %in% state_multimonth$id],
         file.path(output_dir, "original_full_rows_state_based.csv"))
  rm(full_ucdp)
  gc(FALSE)
}

annual <- state_events[, .(
  events = .N,
  multi_events = sum(multi_month),
  total_fatalities = sum(best),
  multi_fatalities = sum(best[multi_month])
), by = year][order(year)]
annual[, `:=`(
  multi_event_share = multi_events / events,
  multi_fatality_share = multi_fatalities / total_fatalities
)]
fwrite(annual, file.path(output_dir, "yearly_summary.csv"))

by_country <- state_multimonth[, .(
  multi_events = .N,
  multi_fatalities = sum(best),
  max_months = max(months_spanned),
  max_deaths = max(best)
), by = country][order(-multi_fatalities)]
fwrite(by_country, file.path(output_dir, "country_summary_multimonth.csv"))
by_precision <- state_multimonth[, .(
  multi_events = .N,
  multi_fatalities = sum(best)
), by = date_prec][order(date_prec)]
fwrite(by_precision, file.path(output_dir, "date_precision_summary.csv"))

cat("3) 6 Parquet-Dateien laden und country_id zuordnen ...\n")
pq <- rbindlist(lapply(parquet_files, function(f) as.data.table(arrow::read_parquet(f))),
                use.names = TRUE)
pq[, `:=`(outcome = as.numeric(outcome),
          month_id = as.integer(month_id),
          country_id = as.integer(country_id))]
pq[, `:=`(
  year = 1980L + (month_id - 1L) %/% 12L,
  month = 1L + (month_id - 1L) %% 12L
)]
pq[, month_date := as.Date(sprintf("%04d-%02d-01", year, month))]
stopifnot(all(pq$month_id %in% 457:528))
if (anyDuplicated(pq, by = c("month_id", "country_id"))) {
  stop("Parquet: Doppelte country-month Schlüssel gefunden")
}

# Parquet-country_id ist kein UCDP-country_id! Landnamen werden durch
# Abgleich der gesamten 72-monatigen Todesfallserien ermittelt.
p_ids <- sort(pq[, .(total = sum(outcome)), by = country_id][total > 0, country_id])
u_countries <- sort(state_events[, .(total = sum(best)), by = country][total > 0, country])
if (length(p_ids) > length(u_countries)) {
  stop("Zu viele aktive Parquet-Länder für 1:1-Abgleich")
}
month_ids <- 457:528
P <- matrix(0, nrow = length(p_ids), ncol = length(month_ids),
            dimnames = list(as.character(p_ids), as.character(month_ids)))
p_active <- pq[country_id %in% p_ids]
P[cbind(match(p_active$country_id, p_ids), match(p_active$month_id, month_ids))] <-
  p_active$outcome

U <- matrix(0, nrow = length(u_countries), ncol = length(month_ids),
            dimnames = list(u_countries, as.character(month_ids)))
u_group <- state_events[end_month_id %in% month_ids,
                        .(deaths = sum(best)), by = .(country, end_month_id)]
U[cbind(match(u_group$country, u_countries),
        match(u_group$end_month_id, month_ids))] <- u_group$deaths

cost <- vapply(seq_len(nrow(U)), function(j) {
  rowSums(abs(sweep(P, 2, U[j, ], FUN = "-")))
}, FUN.VALUE = numeric(nrow(P)))
assignment <- as.integer(clue::solve_LSAP(cost))
country_map <- data.table(
  parquet_country_id = p_ids,
  country = u_countries[assignment],
  absolute_monthly_difference = cost[cbind(seq_along(p_ids), assignment)],
  parquet_total = rowSums(P),
  ucdp_total = rowSums(U)[assignment]
)
country_map[, relative_difference := absolute_monthly_difference /
              pmax(parquet_total, ucdp_total, 1)]
setorder(country_map, -parquet_total)
fwrite(country_map, file.path(output_dir, "country_id_mapping_inferred.csv"))
cat("Zuordnungen mit weniger als 5% relativer Differenz:",
    sum(country_map$relative_difference < 0.05), "/", nrow(country_map), "\n")

pq_named <- merge(pq, country_map[, .(parquet_country_id, country)],
                  by.x = "country_id", by.y = "parquet_country_id", all = FALSE)
u_end <- state_events[, .(ucdp_endmonth = sum(best)),
                      by = .(country, month_id = end_month_id)]
comparison <- merge(pq_named, u_end, by = c("country", "month_id"), all.x = TRUE)
comparison[is.na(ucdp_endmonth), ucdp_endmonth := 0]
comparison[, difference := outcome - ucdp_endmonth]
fwrite(comparison, file.path(output_dir, "country_month_comparison.csv"))

cat("4) Sensitivitätsanalyse: Monatsende vs. gleiche Verteilung ...\n")
# Einzelmonats-Events beibehalten, nur die mehrmonatigen verteilen.
single_month <- state_events[multi_month == FALSE,
                             .(fatalities_equal_split = sum(best)),
                             by = .(country, month_date = month_end)]
spanning_month <- state_events[multi_month == TRUE, {
  m <- seq.Date(month_start[1L], month_end[1L], by = "month")
  .(month_date = m, fatalities_equal_split = best[1L] / length(m))
}, by = .(id, country)]
equal_month <- rbindlist(list(single_month,
                               spanning_month[, .(country, month_date, fatalities_equal_split)]))[
  , .(fatalities_equal_split = sum(fatalities_equal_split)),
  by = .(country, month_date)]
end_month <- state_events[, .(fatalities_end_month = sum(best)),
                          by = .(country, month_date = month_end)]
observed <- pq_named[, .(parquet_outcome = sum(outcome)), by = .(country, month_date)]
monthly <- merge(end_month, equal_month, by = c("country", "month_date"), all = TRUE)
monthly <- merge(monthly, observed, by = c("country", "month_date"), all = TRUE)
for (col in c("fatalities_end_month", "fatalities_equal_split", "parquet_outcome")) {
  set(monthly, which(is.na(monthly[[col]])), col, 0)
}
monthly <- monthly[month_date >= as.Date("2018-01-01") &
                   month_date <= as.Date("2023-12-01")]
fwrite(monthly, file.path(output_dir, "monthly_end_vs_even_split_vs_parquet.csv"))

cat("5) Diagramme exportieren ...\n")
theme_set(theme_minimal(base_size = 12))
fmt <- scales::label_number(big.mark = ".", decimal.mark = ",")
saveplot <- function(plot, fname, w = 10, h = 5) {
  ggsave(filename = file.path(plot_dir, fname), plot = plot,
         width = w, height = h, dpi = 180, bg = "white")
}

p1 <- ggplot(annual, aes(x = factor(year), y = 100 * multi_fatality_share)) +
  geom_col(fill = "#3075a8") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Fatalities aus mehrmonatigen state-based Events",
       x = "Jahr", y = "Anteil aller Fatalities (%)")
saveplot(p1, "01_fatality_share_by_year.png", 9, 5)

p2 <- ggplot(head(by_country, 12),
             aes(x = reorder(country, multi_fatalities), y = multi_fatalities)) +
  geom_col(fill = "#3e819a") + coord_flip() +
  scale_y_continuous(labels = fmt) +
  labs(title = "Fatalities in mehrmonatigen Events nach Land",
       x = NULL, y = "Fatalities (best)")
saveplot(p2, "02_countries.png", 10, 6)

top_events <- head(copy(state_multimonth)[order(-best)], 12L)
top_events[, event_label := paste0(country, " | ", format(date_start, "%m/%Y"),
                                   "–", format(date_end, "%m/%Y"), " | ", id)]
p3 <- ggplot(top_events, aes(x = reorder(event_label, best), y = best)) +
  geom_col(fill = "#4e81a2") + coord_flip() +
  scale_y_continuous(labels = fmt) +
  labs(title = "12 Events mit den meisten Fatalities über mehrere Monate",
       x = NULL, y = "Fatalities (best)")
saveplot(p3, "03_top_events.png", 12, 7)

p4 <- ggplot(by_precision, aes(x = factor(date_prec), y = multi_fatalities)) +
  geom_col(fill = "#677b94") + scale_y_continuous(labels = fmt) +
  labs(title = "Fatalities nach UCDP date_prec", x = "date_prec", y = "Fatalities")
saveplot(p4, "04_date_precision.png", 8, 5)

ethiopia <- monthly[country == "Ethiopia" &
                     month_date >= as.Date("2021-01-01") &
                     month_date <= as.Date("2022-12-01")]
ethiopia_long <- melt(ethiopia, id.vars = "month_date",
                      measure.vars = c("parquet_outcome", "fatalities_end_month",
                                       "fatalities_equal_split"),
                      variable.name = "Serie", value.name = "Fatalities")
p5 <- ggplot(ethiopia_long, aes(month_date, Fatalities, colour = Serie)) +
  geom_line(linewidth = 0.85) +
  scale_y_continuous(labels = fmt) +
  labs(title = "Äthiopien 2021–2022: Endmonat vs. gleichmäßige Verteilung",
       x = "Monat", y = "Fatalities")
saveplot(p5, "05_ethiopia.png", 11, 5)

ukraine <- monthly[country == "Ukraine" &
                    month_date >= as.Date("2022-01-01") &
                    month_date <= as.Date("2023-12-01")]
ukraine_long <- melt(ukraine, id.vars = "month_date",
                     measure.vars = c("parquet_outcome", "fatalities_end_month"),
                     variable.name = "Serie", value.name = "Fatalities")
p6 <- ggplot(ukraine_long, aes(month_date, Fatalities, colour = Serie)) +
  geom_line(linewidth = 0.85) + scale_y_continuous(labels = fmt) +
  labs(title = "Ukraine: Parquet vs. UCDP im Endmonat",
       x = "Monat", y = "Fatalities")
saveplot(p6, "06_ukraine_validation.png", 11, 5)

cat("\nAnalyse fertig. Ergebnisse unter:\n", normalizePath(output_dir), "\n")
cat("Hinweis: Die letzte mehrmonatige Event-ID ist:", excluded_id, "\n")
cat("Ausgabe-Country-IDs sind anhand von Fatality-Zeitreihen geschätzt und zu prüfen.\n")
