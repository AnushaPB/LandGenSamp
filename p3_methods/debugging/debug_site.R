set.seed(42)

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")
source("MMRR_functions.R")

#register cores
cores <- 20
cl <- makeCluster(cores)
registerDoParallel(cl)

res_mmrr <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("here", "vcfR", "adegenet", "stringr")) %dopar% {

  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #skip iteration if files do not exist
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  
  #run MMRR
    gsd_df <- get_data(i, params = params, "gsd")
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsite)
        subgsd_df <- gsd_df[subIDs,]
        
        #get sites
        siteIDs <- get_sites(params[i,], params, sampstrat, nsite)
        #confirm that number of sites matches number of sample IDs
        stopifnot(length(subIDs) == length(siteIDs))
        #calculate env values by site
        sitegsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN=mean)[,-1]) 
        #confirm that number of values matches number of sites
        stopifnot(nsite == nrow(sitegsd_df))
        
      }
  }
  
  return(NULL)
  
}

#stop cluster
stopCluster(cl)
