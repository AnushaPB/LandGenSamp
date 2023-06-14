# Parallel processing libraries
library(furrr)
library(dplyr)
library(tidyr)
library(here)
library(vcfR)
library(lfmm)
library(gdm)
library(vegan)
library(purrr)

# Read in general functions and objects
source(here("general_functions.R"))
source(here("p3_methods", "general_run_functions.R"))
source(here("p3_methods", "GEA_functions.R"))

# run analysis for individual sampling
method = "rda"
ns = nsites
strats = sitestrats
site = TRUE
site_results <- run_analysis(params, ns = ns, strats = strats, method = method, full_result = NULL, site = site, ncores = 10)
path <- here("p3_methods", "outputs", paste0(method, "_sitesampling_results.csv"))
write.csv(site_results, path, row.names = FALSE)
