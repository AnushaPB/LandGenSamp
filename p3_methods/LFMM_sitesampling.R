set.seed(42)

library("here") #paths
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")

##########
#  LFMM  #
##########

#register cores
cores <- 25
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)


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
    
    run_lfmm_params(i, params, path = "outputs/LFMM/lfmm_sitesampling_results", mode = "site")
    
    gc()
  }
)

#stop cluster
stopCluster(cl)

#stop cluster
stopCluster(cl)

write.csv(res_lfmm, "outputs/LFMM_sitesampling/lfmm_sitesampling_results.csv", row.names = FALSE)

