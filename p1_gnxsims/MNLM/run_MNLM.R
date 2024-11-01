library(raster)
library(here)
#install.packages("devtools")
#devtools::install_github("ropensci/NLMR")
library(NLMR)
library(tidyverse)
source(here("p1_gnxsims", "MNLM", "MNLM_functions.R"))

combos <- expand.grid(seed = c(1, 2, 3), H = c(0.05, 0.5), r = c(0.3, 0.6))
mnlms <- pmap(combos, mnlm_create, .progress = TRUE)
walk(mnlms, mnlm_write)