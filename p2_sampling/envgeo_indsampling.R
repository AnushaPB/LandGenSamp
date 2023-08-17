
library("here")
library("foreach")
library("doParallel")
library("vegan")
library("tidyverse")
source(here("general_functions.R"))
source(here("p2_sampling", "sampling_functions.R"))

#register cores
cores <- 20
cl <- makeCluster(cores) 
registerDoParallel(cl)

for(n in nsamps){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("vegan")
    library("tidyverse")
    #create file path
    gsd_filepath <- create_filepath(i, params = params, "gsd")
    
    #skip iteration if file does not exist
    skip_to_next <- FALSE
    if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
    if(skip_to_next) { print("File does not exist:")
      print(params[i,]) } 
    if(skip_to_next) { result <- NA } 
    
    #run sampling
    if(skip_to_next == FALSE){
      gsd_df <- get_gsd(gsd_filepath)
      set.seed(2)
      samples <- envgeo_indsamp(gsd_df, npts = n, Nreps = 1000)
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("envgeo",1:ncol(samples))
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, samples[,i])}
  colnames(samp_out) <- c(colnames(params),colnames(samples))
  write.csv(samp_out, paste0("outputs/samples_envgeo",n,".csv"), row.names = FALSE)
}


#stop cluster
stopCluster(cl)
  


