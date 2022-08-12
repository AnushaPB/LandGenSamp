set.seed(42)

library(here) #paths
library(vegan) #RDA
library(vcfR)  #read VCF files
library(robust)
library(qvalue)
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("RAD_functions.R")

#register cores
cores <- 20
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

res_rda <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("vcfR", "vegan", "here", "stringr", "tidyverse", "qvalue", "robust", "dplyr")) %dopar% {
  #skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if(skip_to_next) { result <- NA } 
  
  #run RDA
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #make data.frame
    result <- data.frame()
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_rda(subgen, subgsd_df, loci_df)
        
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

#stop cluster
stopCluster(cl)

write.csv(res_rda, "outputs/rda_indsampling_results.csv", row.names = FALSE)
