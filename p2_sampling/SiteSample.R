
SiteSample <- function(sample_sites, coords, npts, buffer_size = NULL, method = "buffer"){
  if(method == "near"){site_samples <- SiteSampleNear(sample_sites, coords, npts)}
  if(method == "buffer"){site_samples <- SiteSampleBuffer(sample_sites, coords, npts, buffer_size = buffer_size)}
}

SiteSampleBuffer <- function(sample_sites, coords, npts, buffer_size = 5){
  #sample_sites - coordinates of sampling sites
  #coords - coordinates of all individuals in data set
  #npts - number of points to sample from each site
  #buffer - buffer around site from which to draw samples randomly
  #edge_buffer - buffer from landscape edges to prevent sampling of sites
  site_samples <- data.frame()
  
  #sample sites can be provided as either sp or coords in vector format (nloop just counts how many sites)
  if(class(sample_sites)[1] == "SpatialPoints"){nloop <- length(sample_sites)} else {nloop <- nrow(sample_sites)}
  
  for(s in 1:nloop){
    #create buffer around sites from which to sample points
    site_buffers <- gBuffer(sample_sites[s,], width = buffer_size)
    #subset out only coordinates falling within site buffers
    buffer_samples <- coords[site_buffers,]
    #convert from SPDF to df
    buffer_samples_df <- data.frame(buffer_samples)
    #randomly sample samples from coordinates
    randsamp <- sample(nrow(buffer_samples_df), npts)
    #subset df
    buffer_samples_df <- buffer_samples_df[randsamp,]
    #IMPORTANT: remove sample from coords so they are not sampled twice
    coords <- coords[!coords$idx %in% buffer_samples_df$idx,]
    
    #THINK ABOUT THIS CODE
    #calculate the mean x coord (psuedo-site)
    #buffer_samples_df$xsite <- mean(buffer_samples_df$x)
    #calculate the mean y coord (psuedo-site)
    #buffer_samples_df$ysite <- mean(buffer_samples_df$y)
    
    #save site IDs
    buffer_samples_df$site <- s
    
    #bind samples
    site_samples <- rbind(site_samples, buffer_samples_df)
  }
  return(site_samples)
}

SiteSampleNear <- function(sample_sites, coords, npts){
  #sample_sites - coordinates of sampling sites
  #coords - coordinates of all individuals in data set
  #npts - number of points to sample from each site
  site_samples <- data.frame()
  
  #sample sites can be provided as either sp or coords in vector format (nloop just counts how many sites)
  if(class(sample_sites)[1] == "SpatialPoints"){nloop <- length(sample_sites)} else {nloop <- nrow(sample_sites)}
  
  for(s in 1:nloop){
    # get site coords
    site_coords <- sample_sites[s,]
    # creates a vector of distances between the site and all other points in the dataset
    dist_vec <- sqrt((coords$x - site_coords$x)^2 + (coords$y - site_coords$y)^2)
    # remove the distance of 0 (same site)
    dist_vec <- dist_vec[dist_vec != 0]
    dist_vec <- dist_vec[order(dist_vec)]
    sample_idx <- names(dist_vec)[1:npts]
    #IMPORTANT: remove sample from coords so they are not sampled twice
    coords <- coords[!coords$idx %in% sample_idx,]
    
    #save site IDs
    sample_df <- data.frame(site = s, sample_idx)
    
    #bind samples
    site_samples <- rbind(site_samples, sample_df)
  }
  return(site_samples)
}
