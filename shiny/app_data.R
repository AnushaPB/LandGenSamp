library(shiny)
library(tidyverse)
library(here)
source(here("p4_analysis", "analysis_functions.R"))


combos <- expand_grid(method = c("lfmm", "rda", "gdm", "mmrr", "mmrr2", "gdm2"), sampling = c("individual", "site"))
clean_data <- pmap(combos, \(method, sampling){
  df <- format_data(method, sampling, p_filter = FALSE)
  write_csv(df, here("shiny", "app_data", paste0(method, "_", sampling, ".csv")))
  })
