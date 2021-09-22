source("general_functions.R")
source("site_functions.R")
library("here")
library("foreach")
library("doParallel")

set.seed(42)

cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

nsite <- c(9, 18, 36)
npts <- 10


for(n in nsite){
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
      
      #randomly select points to act as sites
      #buffer away from edges
      coords_buffer <- crop(coords, extent(5,35,5,35))
      #plot(coords)
      #points(coords_buffer, col="red")
      sample_sites <- coords_buffer[sample(1:length(coords_buffer), n),]

      site_samples <- data.frame()
      for(s in 1:n){
        #create buffer around sites from which to sample points
        site_buffers <- buffer(sample_sites[s,], 300000)
        #subset out only coordinates falling within site buffers
        buffer_samples <- coords[site_buffers,]
        #convert from SPDF to df
        buffer_samples_df <- data.frame(buffer_samples)
        #randomly sample samples from coordinates
        buffer_samples_df <- buffer_samples_df[sample(nrow(buffer_samples_df), npts),]
        
        #THINK ABOUT THIS CODE
        #calculate the mean x coord (psuedo-site)
        #buffer_samples_df$xsite <- mean(buffer_samples_df$x)
        #calculate the mean y coord (psuedo-site)
        #buffer_samples_df$ysite <- mean(buffer_samples_df$y)
        
        #bind samples
        site_samples <- rbind(site_samples, buffer_samples_df)
        site_samples$xsite
      }
      
      #plot (for debugging)
      plot(sample_sites, xlim = c(0,40), ylim = c(0,40))
      points(gsd_df[,c("x","y")], col = "gray")
      points(site_samples[,c("x","y")], col = "red")
      #points(site_samples[,c("xsite","ysite")], col = "blue", pch = 19)
      
      samples <- site_samples$idx 
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("rand", 1:ncol(samples))
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, samples[,i])}
  colnames(samp_out) <- c(colnames(params), colnames(samples))
  write.csv(samp_out, paste0("outputs/site_samples_rand",n,".csv"), row.names = FALSE)
}

#stop cluster
stopCluster(cl)

