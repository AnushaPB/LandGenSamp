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
source(here("p3_methods", "IBDIBE_functions.R"))

# run analysis for individual sampling
method = "mmrr"
ns = nsamps
strats = sampstrats
site = FALSE
safe_save <- safely(save)
full_result <- run_full(params, method = method, ncores = 20, n = 2000) 
safe_save(full_result, file = "full_result.rda")
ind_results <- run_analysis(params, ns = nsamps, strats = sampstrats, method = method, full_result = full_result, site = FALSE, ncores = 25)
path <- here("p3_methods", "outputs", paste0(method, "_indsampling_results.csv"))
write.csv(ind_results, path, row.names = FALSE)
