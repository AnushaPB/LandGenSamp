set.seed(42)
#paths
library("here") 
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")
source("LFMM_functions.R")

#register cores
cores <- 30
cl <- makeCluster(cores)
registerDoParallel(cl)

system.time(
  res_lfmm <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("here", "vcfR", "lfmm", "stringr", "AssocTests", "adegenet", "purrr", "dplyr")) %dopar% {
    
    #skip iteration if files do not exist
    skip_to_next <- skip_check(i, params)
    if(skip_to_next) { result <- NA } 
    
    #run LFMM
    if(skip_to_next == FALSE){
      gen <- get_data(i, params = params, "gen")
      gsd_df <- get_data(i, params = params, "gsd")
      loci_df <- get_data(i, params = params, "loci")
      
      # make data.frame
      result <- data.frame()
      
      for(nsamp in npts){
        for(sampstrat in sampstrats){
          #subsample from data based on sampling strategy and number of samples
          subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
          subgen <- gen[subIDs,]
          subgsd_df <- gsd_df[subIDs,]
          
          #run analysis using subsample
          sub_result <- 
            cross(list(K_selection = c("tracy.widom", "find.clusters", "quick.elbow"), method = c("lasso", "ridge"))) %>%
            map_dfr(run_lfmm_helper, gen = subgen, gsd_df = subgsd_df, loci_df = loci_df)
          
          #save and format new result
          sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsamp, sub_result)
          
          #bind results
          result <- bind_rows(result, sub_result)
          
        }
      }
    }
    
    return(result)
    
    gc()
  }
)


#stop cluster
stopCluster(cl)

write.csv(res_lfmm, "outputs/lfmm_indsampling_results.csv", row.names = FALSE)

