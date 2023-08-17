library("here")
library("foreach")
library("doParallel")
library("raster")

source(here("p2_sampling", "site_functions.R"))
source(here("general_functions.R"))

set.seed(42)

#register cores
cores <- 10
cl <- makeCluster(cores) 
registerDoParallel(cl)

#confirm correct ldim
print(ldim)

# confirm correct nsites
print(nsites)

for (n in nsites){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("raster")
    library("rgeos")
    library("dplyr")
    
    source(here("p2_sampling", "site_functions.R"))
    source(here("general_functions.R"))
    
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
      gsd_df <- get_gsd(gsd_filepath)
      set.seed(5)
      samples <- SiteSample(gsd_df, nsite = n, npts = global_npts, site_method = "equi", sample_method = "near", ldim = ldim)
      # print("nsite", n)
      # print("nsamp", length(samples))
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("equi", 1:ncol(samples))
  #sample IDs
  sampleIDs <- gsub("\\_.*","",samples)
  #site IDs
  siteIDs <- gsub("^.*\\_","", samples)
  
  #create df of sample IDs
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, sampleIDs[,i])}
  colnames(samp_out) <- c(colnames(params), colnames(samples))
  write.csv(samp_out, paste0("outputs/site_samples_equi",n,".csv"), row.names = FALSE)
  
  #create df of site IDs
  site_out <- params
  for(i in 1:ncol(samples)){site_out <- cbind.data.frame(site_out, siteIDs[,i])}
  colnames(site_out) <- c(colnames(params), colnames(samples))
  write.csv(site_out, paste0("outputs/site_ids_equi",n,".csv"), row.names = FALSE)
  
  message(n, "complete")
}


#stop cluster
stopCluster(cl)
  


