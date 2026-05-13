## clear environment
rm(list = ls())

## set working directory
current_path <- rstudioapi::getActiveDocumentContext()$path # get the path of your current open file
setwd(dirname(current_path))



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
library(distfromq)
library(readxl)


# Erkenne das Betriebssystem
os <- Sys.info()["sysname"]

# Setze den Pfad abhängig vom Betriebssystem
# data_path <- "../Data/"
data_path <- ifelse(os == "Windows",
                    "//stat-meth-file1.stat.kit.edu/share-alle/Data/WNV forecasting challenge 2020/",  # Windows
                    "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/WNV forecasting challenge 2020/")  # macOS/Linux



df_actual <- read.csv(paste0(data_path, "West Nile virus human and non-human activity by area of residence for year selected below_.csv"))


### Surveillance data are reported by county of residence, not the location (county or state) of exposure.
# https://www.cdc.gov/west-nile-virus/data-maps/historic-data.html

### achtung nur 480 counties im Datensatz: ist es üblich die zeros einfach nicht zu reporten?
df_actual <- df_actual %>%
  filter(Year == "2020") %>%
  select(-Activity, 
         -Total.human.disease.cases,
         -X..Presumptive.viremic.blood.donors, 
         -Notes)

geocodes <- read_excel(paste0(data_path, "all-geocodes-v2020.xlsx"), skip = 4)

# all files that begin with "2020-04-30"
prediction_files <- list.files(path = data_path, 
                               pattern = "^2020-04-30.*\\.csv$", 
                               full.names = TRUE)

test_df <- df_predictions %>% 
  filter(team == "ARS",
         type == "Point")


# read all files and bind them into a single data frame
df_predictions <- do.call(rbind, lapply(prediction_files, read.csv))


model_names <- c("ARS", "LANL", "MHC", "MSSM", "NCSU", "NYSW", 
                       "NYSW-CVD", "Rutgers", "Standford", "UA", 
                       "UCD", "UI", "UI-NCSA", "WDH")

n_models <- length(model_names)












## ------------------------------------------------------------------------------------------------------------
## define lower threshold for binary event of a present outbreak
## -----
lower_infections_thresh = 0










# upper bounds of the bins (quantiles)
qs <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 150, 200, 999)











































### example

example_forecast_ARS <- df_predictions$value[17:31]
#example_forecast_ARS <- df_predictions$value[1:15]
bin_probs <- example_forecast_ARS



# cdf (ps)
ps <- cumsum(bin_probs)
# WICHTIG: Falls die Summe durch Rundungsfehler minimal über 1 liegt, auf 1 begrenzen
ps[ps > 1] <- 1


x <- seq(from = -0.5, to = 1000, length = 1000)

p_lognormal_approx <- make_p_fn(ps = ps,
                                qs = qs,
                                tail_dist = "lnorm")
cdf_lognormal_approx <- p_lognormal_approx(x)

data.frame(
  x = x,
  y = cdf_lognormal_approx
) %>%
  ggplot() +
  geom_line(
    mapping = aes(x = x, y = y),
    size = 0.8
  ) +
  geom_point(
    data = data.frame(q = qs, p = ps),
    mapping = aes(x = q, y = p),
    size = 1.2
  ) +
  scale_x_log10() + 
  ylim(0, 1) +
  ylab("Probability") +
  xlab("") +
  theme_bw() + 
  xlim(0,20)




x_counts <- 0:999  # Alternativ: seq(from = 0, to = 1000, by = 1)

# 2. Werte die CDF an genau diesen ganzzahligen Stellen aus
cdf_at_counts <- p_lognormal_approx(x_counts)

# Wenn du für weitere Berechnungen die Wahrscheinlichkeiten für JEDEN EINZELNEN Count brauchst
# (also P(X=0), P(X=1), P(X=2) etc. aus der Treppenfunktion ableiten willst):
# Das ist die Differenz zwischen dem aktuellen und dem vorherigen CDF-Wert.
prob_at_counts <- c(cdf_at_counts[1], diff(cdf_at_counts))

# 3. Speichere das Ergebnis übersichtlich in einem Data Frame
df_counts <- data.frame(
  Count = x_counts,
  CDF = cdf_at_counts,       # P(X <= Count)
  Probability = prob_at_counts # P(X == Count)
)


actual <- 0


indicator <- as.numeric(x_counts >= actual)
weight <- 1 / (x_counts + 1)

# RPS discrete version of CRPS
rps <- sum((df_counts$CDF - indicator)^2)
# twRPS discrete version of twCRPS
twrps <- sum(weight * (df_counts$CDF - indicator)^2)

# BRIER Score for Onset meaning > 0 WNVND infections




## für alle predictions aus der CDF (die für alle auf die selbe Art und Weise erzeugt wird) die WKs
# für die count-valued outcomes ablesen( 0-999 ). Das ist implizit die Treppenfunktion (und somit CRPS = RPS).
# Wie genau CDF erzeugt wird laut Johannes nicht so wichtig.





