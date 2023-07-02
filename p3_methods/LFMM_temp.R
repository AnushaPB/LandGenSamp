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
library(doParallel)
library(foreach)
# Read in general functions and objects
source(here("general_functions.R"))
source(here("p3_methods", "general_run_functions.R"))
source(here("p3_methods", "GEA_functions.R"))

# run analysis for individual sampling
method = "lfmm"
ns = nsamps
strats = sampstrats
site = FALSE
ind_results <- run_analysis(params, ns = ns, strats = strats, method = method, full_result = NULL, site = site, ncores = 25)
path <- here("p3_methods", "outputs", paste0(method, "_indsampling_results.csv"))write.csv(ind_results, path, row.names = FALSE)
