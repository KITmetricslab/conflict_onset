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





















# example prediction
bin_probs <- c(
  0.572377761,    # Bin 0
  0.426117879,    # Bin 1-5
  0.001459964,    # Bin 6-10
  4.11e-05,       # Bin 11-15
  2.88e-06,       # Bin 16-20
  3.45e-07,       # Bin 21-25
  5.68e-08,       # Bin 26-30
  1.17e-08,       # Bin 31-35
  2.82e-09,       # Bin 36-40
  1.77e-10        # Bin 41-45
)

# 2. upper bounds of the bins (quantiles)
qs <- c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)

# 3. cumulativ distribution (ps)
ps <- cumsum(bin_probs)

# WICHTIG: Falls die Summe durch Rundungsfehler minimal über 1 liegt, auf 1 begrenzen
ps[ps > 1] <- 1


x <- seq(from = -1, to = 50, length = 1000)

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
  ylim(0, 1) +
  ylab("Probability") +
  xlab("") +
  theme_bw()



plot_ps <- seq(0.001, 0.999, length.out = 1000)

length(plot_ps)

q_lognormal_approx <- make_q_fn(ps = ps,
                                qs = qs,
                                tail_dist = "lnorm")
qf_lognormal_approx <- q_lognormal_approx(plot_ps)

# Rundung problematisch? damit Wk für x=0 nicht mehr korrekt sondern höher
qf_lognormal_approx <- round(qf_lognormal_approx, 0)

data.frame(
  x = plot_ps,
  y = qf_lognormal_approx
) %>%
  ggplot() +
  geom_line(
    mapping = aes(x = x, y = y),
    size = 0.8
  ) +
  geom_point(
    data = data.frame(q = ps, p = qs),
    mapping = aes(x = q, y = p),
    size = 1.2
  ) +
  xlim(0, 1) +
  xlab("Probability") +
  ylab("") +
  theme_bw()



