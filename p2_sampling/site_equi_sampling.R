source("general_functions.R")
source("site_functions.R")
library("here")
library("foreach")
library("doParallel")
library("raster")

set.seed(42)

#function to make equidistant sampling sites
equi_samp <- function(nsite, ldim){
  #pts - coords to sample from
  #nsite - number of points (or sites) to sample (should be a perfect square)
  #ldim - landscape dimension of one side (landscape should be a square)
  inc <- ldim/(sqrt(nsite)+1)
  xgrid <- ygrid <- seq(0+inc, ldim-inc, inc) 
  cgrid <- expand.grid(xgrid, ygrid)
  colnames(cgrid) <- c("x","y")
  coordinates(cgrid) <- ~x+y

  plot(cgrid)
  return(cgrid)
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

#confirm correct ldim
print(ldim)

for(n in nsites){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("raster")
    
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
      pts <- gsd_df[,c("idx","x","y")]
      coords <- pts
      coordinates(coords) <- ~x+y
 
      #equidistant sample sites
      sample_sites <- equi_samp(nsite = n, ldim = ldim)
      
      #sample from around sites based on a buffer
      #chosen arbitrarily, lower was too small/not enough points in buffer for smaller sample sizes
      site_samples <- SiteSample(sample_sites, coords, npts = 10, buffer_size = 600000)
      
      #plot (for debugging)
      #par(pty="s")
      #plot(gsd_df[,c("x","y")], xlim = c(0,ldim), ylim = c(0,ldim), col = "gray")
      #points(sample_sites, pch = 3)
      #points(site_samples[,c("x","y")], col = "red")
      #points(site_samples[,c("xsite","ysite")], col = "blue", pch = 19)
      
      samples <- paste0(site_samples$idx, "_", site_samples$site)
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
  
}


#stop cluster
stopCluster(cl)
  


