set.seed(42)
#paths
library("here") 
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")
library("AssocTests")
library("purrr")

#read in general functions and objects
source("general_functions.R")

##########
#  LFMM  #
##########

#register cores
cores <- 25
cl <- makeCluster(cores)
#not to overload your computer
registerDoParallel(cl)

system.time(
res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  library("stringr")
  
  run_lfmm_params(i, params, path = "outputs/LFMM/lfmm_results", mode = "ind")
  
  gc()
}
)

#stop cluster
stopCluster(cl)

write.csv(res_lfmm, "outputs/LFMM/lfmm_results.csv", row.names = FALSE)

