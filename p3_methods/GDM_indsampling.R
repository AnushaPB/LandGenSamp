set.seed(42)

library("here") #paths
library("gdm") #GDM
library("vcfR")
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("GDM_functions.R")


###########
#   GDM   #
###########

#register cores
cores <- 20
cl <- makeCluster(cores) 
registerDoParallel(cl)


res_gdm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  #vcfR
  library("vcfR")
  library("gdm")
  library("adegenet")
  
  results <- run_gdm_params(i, params, "outputs/GDM/gdm_results", mode = "ind")
  
  return(results)
  
}

#stop cluster
stopCluster(cl)

write.csv(res_gdm, "outputs/gdm_results_euc.csv", row.names = FALSE)

