library("foreach")
library("doParallel")
library("here")
source(here("general_functions.R"))
source(here("p2_sampling", "sampling_functions.R"))

set.seed(42)

cores <- 20
cl <- makeCluster(cores) 
registerDoParallel(cl)

for (n in nsites) {
  samples <- foreach(i = 1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("raster")
    library("sf")
    library("dplyr")    
    
    #create file path
    gsd_filepath <- create_filepath(i, params = params, "gsd")
    
    #skip iteration if file does not exist
    skip_to_next <- FALSE
    if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
    if(skip_to_next) { print("File does not exist:")
      print(params[i,]) } 
    if(skip_to_next) { samples <- NA } 
    
    #run sampling
    if(skip_to_next == FALSE){
      #get data
      gsd_df <- get_gsd(gsd_filepath)
      #sample
      set.seed(4)
      samples <- SiteSample(gsd_df, nsite = n, npts = global_npts, site_method = "rand", sample_method = "near")
    }
    
    #return vector of sample IDs
    return(samples)
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("rand", 1:ncol(samples))
  #sample IDs
  sampleIDs <- gsub("\\_.*","",samples)
  #site IDs
  siteIDs <- gsub("^.*\\_","", samples)
  
  #create df of sample IDs
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, sampleIDs[,i])}
  colnames(samp_out) <- c(colnames(params), colnames(samples))
  write.csv(samp_out, paste0("outputs/site_samples_rand",n,".csv"), row.names = FALSE)
  
  #create df of site IDs
  site_out <- params
  for(i in 1:ncol(samples)){site_out <- cbind.data.frame(site_out, siteIDs[,i])}
  colnames(site_out) <- c(colnames(params), colnames(samples))
  write.csv(site_out, paste0("outputs/site_ids_rand",n,".csv"), row.names = FALSE)
  
  message(paste(n, "complete"))
}

#stop cluster
stopCluster(cl)
