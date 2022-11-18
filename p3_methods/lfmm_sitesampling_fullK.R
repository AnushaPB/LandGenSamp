set.seed(42)

library("here") #paths
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")
library("AssocTests")

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")
source("lfmm_functions.R")

#register cores
cores <- 30
cl <- makeCluster(cores)
registerDoParallel(cl)

system.time(
res_lfmm <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("here", "vcfR", "lfmm", "stringr", "AssocTests", "adegenet", "purrr", "dplyr")) %dopar% {
  
  #skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if(skip_to_next) result <- NA 
  
  #run LFMM
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    # subset and get K
    s <- sample(nrow(gen), 1000)
    gen2k <- gen[s,]
    gsd2k <- gsd_df[s,]
    K <- get_K(gen2k, coords = gsd2k[,c("x", "y")], method = "tess") 

    # make data.frame
    result <- data.frame()
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        # subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsite)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        # get sites
        siteIDs <- get_sites(params[i,], sampstrat, nsite)
        # confirm that number of sites matches number of sample IDs
        stopifnot(length(subIDs) == length(siteIDs))
        # calculate allele frequency by site (average)
        sitegen <- data.frame(aggregate(subgen, list(siteIDs), FUN = mean)[,-1])
        # calculate env values by site
        sitegsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN = mean)[,-1]) 
        
        # run analysis using subsample
        sub_result <- 
          cross(list(K_selection = "full", method = c("ridge"))) %>%
          map_dfr(run_lfmm_helper, gen = sitegen, gsd_df = sitegsd_df, loci_df = loci_df, K = K)
        
        # save and format new result
        sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsite, sub_result)
        
        # bind results
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


write.csv(res_lfmm, "outputs/lfmm_sitesampling_fullK_results.csv", row.names = FALSE)

