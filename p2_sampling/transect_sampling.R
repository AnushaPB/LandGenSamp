source("general_functions.R")
library("here")
library("foreach")
library("doParallel")

transect_samp <- function(pts, npts, ytsct, buffer){
  #pts - dataframe with IDs and coords
  #npts - total number of points to sample (evenly split across transects)
  #buffer - buffer around transects within which points are sampled 
  
  #divide number of samples evenly among the transects
  npts_tsct <- npts/length(ytsct)
    
  #plot all points (gray) (for debugging, comment out later)
  par(pty="s")
  plot(gsd_df$x, gsd_df$y, pch=19, cex=0.2, col="gray", main = npts)
    
  #create empty vector to store IDs
  samples <- c()
  for(i in 1:length(ytsct)){ 
    #subset points around transect based on buffer
    tsctsq <- subset(pts, y > (ytsct[i] - buffer) & y < (ytsct[i] + buffer))
    #randomly sample subset of transect points to match number of samples needed for each transect
    tsctsq <- tsctsq[sample(nrow(tsctsq), npts_tsct),]
    #plot points sampled (for debugging, comment out later)
    points(tsctsq$x, tsctsq$y, col=i+1)
    #store IDs in list
    samples <- c(samples, tsctsq$idx)
  }
  
  #confirm correct number of samples were subsetted
  stopifnot(npts == length(samples))
  
  #return vec of sample IDs
  return(samples)
}

#vec of number of sample points
npts_vec <- c(36, 81, 144, 225, 324)
#horizontal transects (y-coords)
ytsct <- c(10, 20, 30)
#buffer around transects
buffer <- 1


#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

for(n in npts){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    
    #create file path
    gsd_filepath <- create_filepath(i, "gsd")
    
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
      samples <- transect_samp(pts, npts_vec[n], ytsct, buffer)
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("trans",1:ncol(samples))
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, samples[,i])}
  colnames(samp_out) <- c(colnames(params),colnames(samples))
  write.csv(samp_out, paste0("outputs/samples_trans",n,".csv"))
}


#stop cluster
stopCluster(cl)


