# in bash:
# Set the variables
#m="lfmm"
#s="individual"
# Call the Rscript and pass the variables as command-line arguments
#Rscript my_script.R "$m" "$s"

library(furrr)
library(dplyr)
library(tidyr)
library(here)
library(vcfR)
library(lfmm)
library(gdm)
library(vegan)

# Read in general functions and objects
source(here("general_functions.R"))
source(here("p3_methods", "general_run_functions.R"))
source(here("p3_methods", "GEA_functions.R"))
source(here("p3_methods", "IBDIBE_functions.R"))

# Access the variables passed from the command line
args <- commandArgs(trailingOnly = TRUE)
method <- args[1]
sampling <- args[2]

# Print the values of the variables
print(method)
print(sampling)

# Run method
run_method(method, sampling = c("individual", "site"), ncores = 27)