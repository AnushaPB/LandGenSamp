source("general_functions.R")
source("site_functions.R")
library("here")
library("foreach")
library("doParallel")

set.seed(42)

cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

nistes <- c(9, 16, 25)
for(n in nsites){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("raster")
    library("sp")
    
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
      #get data
      gsd_df <- get_gsd(gsd_filepath)
      coords <- gsd_df[,c("idx", "x","y")]
      coordinates(coords) <- ~x+y
      
      
      #buffer away from edges
      coords_buffer <- crop(coords, extent(5,35,5,35))
      #plot(coords)
      #points(coords_buffer, col="red")
      #randomly select points to act as sites
      sample_sites <- coords_buffer[sample(1:length(coords_buffer), n),]

      #sample from around sites based on a buffer
      #400000 chosen arbitrarily, 300000 was too small/not enough points in buffer for smaller sample sizes
      site_samples <- SiteSample(sample_sites, coords, npts = 10, buffer_size = 400000)
      
      #plot (for debugging)
      plot(sample_sites, xlim = c(0,40), ylim = c(0,40))
      points(gsd_df[,c("x","y")], col = "gray")
      points(site_samples[,c("x","y")], col = "red")
      #points(site_samples[,c("xsite","ysite")], col = "blue", pch = 19)
      
      samples <- paste0(site_samples$idx, "_", site_samples$site)
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
  
}

#stop cluster
stopCluster(cl)
