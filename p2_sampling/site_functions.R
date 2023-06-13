# INDIVIDUAL SAMPLING ----------------------------------------------------------------------------------------------------------------------------------------



# SITE SAMPLING ---------------------------------------------------------------------------------------------------------------------------------------------------

SiteSample <- function(gsd_df, nsite, npts, site_method, sample_method = "near",  buffer_size = NULL, edge_buffer = NULL, ldim = NULL, Nreps = 1000){
  # make coords
  coords <- gsd_df[,c("idx","x","y")]
  coordinates(coords) <- ~x+y
  
  # sample sites
  if(site_method == "rand"){sample_sites <- rand_samp(coords = coords, nsite = nsite, buffer_size = buffer_size, edge_buffer = edge_buffer, ldim = ldim)}
  if(site_method == "envgeo"){sample_sites <- envgeo_samp(gsd_df, nsite = nsite, Nreps = Nreps, edge_buffer = global_edge_buffer, ldim = ldim)}
  if(site_method == "equi"){sample_sites <- equi_samp(nsite = nsite, ldim = ldim)}
  
  # sample points around sites 
  if(sample_method == "near"){site_samples <- SiteSampleNear(sample_sites, coords, npts = npts)}
  if(sample_method == "buffer"){site_samples <- SiteSampleBuffer(sample_sites, coords, npts = npts, buffer_size = buffer_size)}
  
  # make a vector
  samples <- paste0(site_samples$idx, "_", site_samples$site)
  
  return(samples)
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
    #check if there are enough samples
    if(length(buffer_samples) < npts){stop(paste0("Less than npts found in buffer around site, try choosing a different site or increasing the buffer size"))}
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
    sample_df <- data.frame(site = s, idx = sample_idx)
    
    #bind samples
    site_samples <- rbind(site_samples, sample_df)
  }
  return(site_samples)
}


rand_samp <- function(coords, nsite, buffer_size = 5, edge_buffer = NULL, ldim = NULL){
  #buffer away from edges if ldim and edge_buffer provided
  if(is.null(ldim) | is.null(edge_buffer)){coords_buffer <- coords} else {coords_buffer <- crop(coords, extent(edge_buffer, ldim-edge_buffer, edge_buffer, ldim-edge_buffer))}
  
  #randomly select points to act as sites
  sample_sites <- coords_buffer[sample(1:length(coords_buffer), nsite),]
  
  return(sample_sites)
}


#function to make equidistant sampling sites
equi_samp <- function(nsite, ldim = 100, buffer = 10){
  #nsite - number of points (or sites) to sample (should be a perfect square)
  #ldim - landscape dimension of one side (landscape should be a square)
  inc <- (ldim - buffer*2)/(sqrt(nsite) - 1)
  xgrid <- ygrid <- seq(0+buffer, ldim-buffer, inc) 
  cgrid <- expand.grid(xgrid, ygrid)
  
  par(pty = "s")
  plot(cgrid, xlim = c(0,ldim), ylim = c(0,ldim))
  
  colnames(cgrid) <- c("x","y")
  sp::coordinates(cgrid) <- ~x+y
  
  
  return(cgrid)
}

# function to perform envgeo sampling
envgeo_samp <- function(gsd_df, nsite, Nreps = 1000, edge_buffer = NULL, ldim = NULL){
  sample.sets <- matrix(nrow=Nreps, ncol=nsite)
  results <- data.frame(env1.var=numeric(Nreps), env2.var=numeric(Nreps),
                        Mantel.r=numeric(Nreps), Mantel.p=numeric(Nreps),
                        mean.dist=numeric(Nreps))
  if(is.null(ldim) | is.null(edge_buffer)){
    gsd_df <- gsd_df
  } else {
    #define buffer
    buffmin <- edge_buffer
    buffmax <- ldim - edge_buffer
    #buffer coordinates away from edge
    gsd_df <- gsd_df[gsd_df$x > buffmin & gsd_df$x < buffmax & gsd_df$y > buffmin & gsd_df$y < buffmax,]
  }

  env.df <- gsd_df[,c("env1","env2")]
  
  e.dist <-  as.matrix(dist(gsd_df[,c("env1","env2")], diag = TRUE, upper = TRUE)) 
  g.dist <- as.matrix(dist(gsd_df[,c("x","y")], diag = TRUE, upper = TRUE))
  for(i in 1:Nreps){
    NN <- sample(1:nrow(gsd_df), nsite, replace = FALSE)
    sample.sets[i,] <- NN
    
    env.sub <- env.df[NN,]
    g.dist.sub <- g.dist[NN, NN]
    e.dist.sub <- e.dist[NN, NN]
    
    results$mean.dist[i] <- mean(g.dist.sub)
    
    results$env1.var[i] <- var(env.sub$env1)
    results$env2.var[i] <- var(env.sub$env2)
    
    DxE <- mantel(g.dist.sub, e.dist.sub, permutations = 99)
    results$Mantel.r[i] <- DxE$statistic
    results$Mantel.p[i] <- DxE$signif
  }
  
  score <- scale(1-results$Mantel.r) + scale(results$env1.var) + scale(results$env2.var)
  best_sample <- sample.sets[which.max(score),]
  sub_df <- gsd_df[best_sample,]
  
  #overwrite sample sites with coordinates for sample sites using indexes
  sample_sites <- gsd_df[as.character(sub_df$idx), c("x","y")]
  
  #convert to coordinates
  coordinates(sample_sites) <- ~x+y
  
  return(sample_sites)
}

#method to use for site sampling (can be "near" or "buffer")
site_method <- "near"
#number of individuals to sample per site
global_npts <- 10
#buffer away from edge for site selection (if using buffer method)
global_edge_buffer <- 5
#buffer size around sites from which individuals (pts) are sampled (if using buffer method)
global_buffer_size <- 5
