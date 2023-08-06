library("here")
library("foreach")
library("doParallel")

source(here("general_functions.R"))

set.seed(42)

grid_samp <- function(pts, npts, ldim){
  # switch coordinates back to positive y to sample grid
  pts$y <- -pts$y
  inc <- ldim/sqrt(npts)
  xgrid <- ygrid <- seq(0, ldim, inc) 
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
  while(nrow(subs) != npts & i < (length(xgrid)-1) & j < (length(ygrid)-1)){
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
cl <- makeCluster(5) 
registerDoParallel(cl)

for(n in nsamps){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
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
      pts <- gsd_df[,c("idx","x","y")]
      set.seed(7)
      samples <- grid_samp(pts, npts = n, ldim = ldim)
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("grid",1:ncol(samples))
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, samples[,i])}
  colnames(samp_out) <- c(colnames(params),colnames(samples))
  write.csv(samp_out, paste0("outputs/samples_grid",n,".csv"), row.names = FALSE)
}


#stop cluster
stopCluster(cl)
  


