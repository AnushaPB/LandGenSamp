
SiteSample <- function(sample_sites, coords, npts, buffer_size = 300000){
  #sample_sites - coordinates of sampling sites
  #coords - coordinates of all individuals in data set
  #npts - number of points to sample from each site
  #buffer - buffer around site from which to draw samples randomly
  #edge_buffer - buffer from landscape edges to prevent sampling of sites
  site_samples <- data.frame()
  for(s in 1:nrow(sample_sites)){
    #create buffer around sites from which to sample points
    site_buffers <- buffer(sample_sites[s,], buffer_size)
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
    
    #save site IDs
    buffer_samples_df$site <- s
    
    #bind samples
    site_samples <- rbind(site_samples, buffer_samples_df)
    site_samples$xsite
  }
  return(site_samples)
}

nsite <- c(9, 14, 25)
npts <- 10

