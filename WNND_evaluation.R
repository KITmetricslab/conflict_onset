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


# Erkenne das Betriebssystem
os <- Sys.info()["sysname"]

# Setze den Pfad abhängig vom Betriebssystem
# data_path <- "../Data/"
data_path <- ifelse(os == "Windows",
                    "//stat-meth-file1.stat.kit.edu/share-alle/Data/WNV forecasting challenge 2020/",  # Windows
                    "smb://stat-meth-file1.stat.kit.edu/share-alle/Data/WNV forecasting challenge 2020/")  # macOS/Linux


predictions_ARS <- paste0(data_path, "2020-04-30_ARS.csv")
data_ARS <- read.csv(predictions_ARS)




prediction_files <- list.files(path = data_path, 
                               pattern = "^2020-.*\\.csv$", 
                               full.names = TRUE)


# read all files and bind them into a single data frame
df_predictions <- do.call(rbind, lapply(prediction_files, read.csv))

# upper bounds of the bins (quantiles)
qs <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 150, 200, 999)



### example

example_forecast_ARS <- df_predictions$value[17:31]
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




## für alle predictions aus der CDF (die für alle auf die selbe Art und Weise erzeugt wird) die WKs
# für die count-valued outcomes ablesen( 0-999 ). Das ist implizit die Treppenfunktion (und somit CRPS = RPS).
# Wie genau CDF erzeugt wird laut Johannes nicht so wichtig.





