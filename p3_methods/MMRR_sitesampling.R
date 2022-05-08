set.seed(42)

library(here) #paths

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")

############
#   MMRR   #
############
#register cores
cores <- 5
cl <- makeCluster(cores) 
registerDoParallel(cl)


res_mmrr <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("vcfR")
  library("mmrr")
  library("adegenet")
  library("here")
  library("tidyverse")
  results <- run_mmrr_params(i, params, sampstrats, npts, "outputs/MMRR/mmrr_sitesampling_results", mode = "site")
  
  return(results)
}

#stop cluster
stopCluster(cl)

write.csv(res_mmrr, "outputs/mmrr_sitesampling_results.csv", row.names = FALSE)

