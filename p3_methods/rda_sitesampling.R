
set.seed(42)

library(here) #paths
library(vegan) #RDA
library(vcfR)  #read VCF files
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source(here("general_functions.R"))
source("sitesampling_functions.R")
source("rda_functions.R")

#register cores
cores <- 25
cl <- makeCluster(cores) 
registerDoParallel(cl)

res_rda <- foreach(i = 1:nrow(params), .combine=rbind, .packages = c("vcfR", "vegan", "here", "stringr", "tidyverse", "qvalue", "robust", "dplyr")) %dopar% {

  #skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if(skip_to_next) { result <- NA } 
  
  #run RDA
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    # make data.frame
    result <- data.frame()
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsite)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #get sites
        siteIDs <- get_sites(params[i,], sampstrat, nsite)
        #confirm that number of sites matches number of sample IDs
        stopifnot(length(subIDs) == length(siteIDs))
        #calculate allele frequency by site (average)
        sitegen <- data.frame(aggregate(subgen, list(siteIDs), FUN=mean)[,-1])
        #calculate env values by site
        sitegsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN=mean)[,-1])
        
        # run pRDA
        # run analysis using subsample
        sub_result_pRDA <- run_rda(sitegen, sitegsd_df, loci_df, correctPC = TRUE)
        # save and format new result
        sub_result_pRDA <- data.frame(params[i, ], sampstrat = sampstrat, nsamp = nsite, correctPC = TRUE, sub_result_pRDA)
        
        # run regular RDA
        # run analysis using subsample
        sub_result_RDA <- run_rda(sitegen, sitegsd_df, loci_df, correctPC = FALSE)
        # save and format new result
        sub_result_RDA <- data.frame(params[i, ], sampstrat = sampstrat, nsamp = nsite, correctPC = FALSE, sub_result_RDA)
        
        # bind results
        result <- bind_rows(result, sub_result_RDA, sub_result_pRDA)
      }
    }
  }
  
  return(result)
  
  gc()
  
}

#stop cluster
stopCluster(cl)

write.csv(res_rda, "outputs/rda_sitesampling_results.csv", row.names = FALSE)
