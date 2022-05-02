library("here") #paths
library("gdm") #GDM
library("vcfR")
#parallel
library(foreach)
library(doParallel)
library(tidyverse)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")
source("GDM_functions.R")

set.seed(42)


###########
#   GDM   #
###########
#register cores
cores <- 2
cl <- makeCluster(cores) 
registerDoParallel(cl)


res_gdm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("vcfR")
  library("gdm")
  library("adegenet")
  
  results <- run_gdm_params(i, params, sampstrats, npts, "outputs/GDM/gdm_sitesampling_results", mode = "site")
  
  return(results)
}

#stop cluster
stopCluster(cl)

write.csv(res_gdm, "sitesampling/outputs/gdm_sitesampling_results_dps.csv", row.names = FALSE)

