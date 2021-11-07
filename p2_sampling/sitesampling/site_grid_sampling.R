source("general_functions.R")
source("sitesampling/site_functions.R")
library("here")
library("foreach")
library("doParallel")
library("raster")

set.seed(42)

grid_samp <- function(pts, nsite, ldim, buffer = 0){
  #pts - coords to sample from
  #nsite - number of points (or sites) to sample (should be a perfect square)
  #ldim - landscape dimension of one side (landscape should be a square)
  inc <- (ldim - 2*buffer)/sqrt(nsite)
  xgrid <- ygrid <- seq(0 + buffer, ldim - buffer, inc) 
  subs <- c()
  #first round of sampling: entire grid
  for(i in 1:(length(xgrid)-1)){ 
    for(j in 1:(length(ygrid)-1)){ 
      gridsq = subset(pts, y > ygrid[j] & y < ygrid[j+1] & x > xgrid[i] & x < xgrid[i+1]) 
      if(dim(gridsq)[1]>0){ subs = rbind(subs, gridsq[sample(1:dim(gridsq)[1],1 ), ]) }
    } 
  }
  #reset grid indices
  i=1
  j=1
  #second round of sampling: cycle through gridcells again until number of desired samples is reached
  while(nrow(subs) != nsite & i < (length(xgrid)-1) & j < (length(ygrid)-1)){
    i = i+1
    j = j+1
    gridsq = subset(pts, y > ygrid[j] & y < ygrid[j+1] & x > xgrid[i] & x < xgrid[i+1]) 
    if(dim(gridsq)[1]>0){subs = rbind(subs, gridsq[sample(1:dim(gridsq)[1],1 ), ])}
  }
  
  #save IDs to vector
  samples <- as.character(subs$idx)
  
  return(samples)
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

nsites <- c(9,16,25)

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
 
      #grid sample sites
      #note - buffer 5 from ldim so sites aren't sampled close to the edge
      sample_sites <- grid_samp(pts, nsite = n, ldim = 40, buffer = 5)
      #overwrite sample sites with coordinates for sample sites using indexes
      sample_sites <- gsd_df[sample_sites, c("x","y")]
      #convert to coordinates
      coordinates(sample_sites) <- ~x+y
      
      #sample from around sites based on a buffer
      #400000 chosen arbitrarily, 300000 was too small/not enough points in buffer for smaller sample sizes
      site_samples <- SiteSample(sample_sites, coords, npts = 10, buffer_size = 400000)
      
      #plot (for debugging)
      par(pty="s")
      plot(gsd_df[,c("x","y")], xlim = c(0,40), ylim = c(0,40), col = "gray")
      points(sample_sites, pch = 3)
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
  write.csv(samp_out, paste0("outputs/site_samples_grid",n,".csv"), row.names = FALSE)
  
  #create df of site IDs
  site_out <- params
  for(i in 1:ncol(samples)){site_out <- cbind.data.frame(site_out, siteIDs[,i])}
  colnames(site_out) <- c(colnames(params), colnames(samples))
  write.csv(site_out, paste0("outputs/site_ids_grid",n,".csv"), row.names = FALSE)
  
}


#stop cluster
stopCluster(cl)
  


