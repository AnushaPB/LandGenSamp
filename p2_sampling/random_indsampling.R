library("here")
library("foreach")
library("doParallel")

source(here("general_functions.R"))
source(here("p2_sampling", "sampling_functions.R"))

cl <- makeCluster(25) 
registerDoParallel(cl)

for (n in nsamps){
  samples <- foreach(i = 1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("tidyverse")
    library("dplyr")
    #create file path
    gsd_filepath <- create_filepath(i, params = params, "gsd")
    
    #skip iteration if file does not exist
    skip_to_next <- FALSE
    if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
    if(skip_to_next) { print("File does not exist:")
      print(params[i,]) } 
    if(skip_to_next) { return(NA) } 
    
    #run sampling
    if (skip_to_next == FALSE){
      gsd_df <- get_gsd(gsd_filepath)
      #take random sample
      set.seed(1)
      samples <- gsd_df[sample(1:nrow(gsd_df), n), "idx"]
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("rand", 1:ncol(samples))
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, samples[,i])}
  colnames(samp_out) <- c(colnames(params), colnames(samples))
  write.csv(samp_out, paste0("outputs/samples_rand",n,".csv"), row.names = FALSE)
}


#stop cluster
stopCluster(cl)
  


