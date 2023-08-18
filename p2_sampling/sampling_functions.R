# INDIVIDUAL SAMPLING ---------------------------------------------------------------------------------------------------------------------------------------------------

# perform EG sampling
envgeo_indsamp <- function(gsd_df, npts, Nreps = 1000){
  Nreps <- 1000
  sample.sets <- matrix(nrow=Nreps, ncol=npts)
  results <- data.frame(env1.var=numeric(Nreps), env2.var=numeric(Nreps),
                        Mantel.r=numeric(Nreps), Mantel.p=numeric(Nreps),
                        mean.dist=numeric(Nreps))
  
  env.df <- gsd_df[,c("env1","env2")]
  e.dist <-  as.matrix(dist(gsd_df[,c("env1","env2")], diag = TRUE, upper = TRUE)) 
  g.dist <- as.matrix(dist(gsd_df[,c("x","y")], diag = TRUE, upper = TRUE))
  for(i in 1:Nreps){
    NN <- sample(1:nrow(gsd_df), npts, replace = FALSE)
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
  
  #save IDs to vector
  samples <- as.character(sub_df$idx)
  
  return(samples)
}

# perform grid sampling
grid_indsamp <- function(pts, npts, ldim){
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

# perform transect sampling
transect_indsamp <- function(pts, npts){
  #pts - dataframe with IDs and coords
  #npts - total number of points to sample (evenly split across transects)

  #horizontal transects (y-coords)
  ytsct <- c(ldim/2 - ldim/4, ldim/2, ldim/2 + ldim/4)
  
  #buffer around transects
  #NOTE: changed from 2 to 3 because a buffer of 2 did not include enough points
  buffer <- 3
  
  #convert y coords back to positive for transect sampling
  pts$y <- -pts$y
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

# SITE SAMPLING ---------------------------------------------------------------------------------------------------------------------------------------------------

SiteSample <- function(gsd_df, nsite, npts, site_method, sample_method = "near",  buffer_size = NULL, edge_buffer = NULL, ldim = 100, Nreps = 1000){
  # make coords
  coords <- gsd_df[,c("idx","x","y")]
  # correct y coords if not corrected
  if (all(gsd_df$y > 0)) gsd_df$y <- -gsd_df$y
  coords <- sf::st_as_sf(coords, coords = c("x","y"))
  
  # sample sites
  if (site_method == "rand") sample_sites <- rand_sitesamp(coords = coords, nsite = nsite, edge_buffer = edge_buffer, ldim = ldim)
  if (site_method == "envgeo") sample_sites <- envgeo_sitesamp(gsd_df, nsite = nsite, Nreps = Nreps, edge_buffer = global_edge_buffer, ldim = ldim)
  if (site_method == "equi") sample_sites <- equi_sitesamp(nsite = nsite, ldim = ldim)
  
  # sample points around sites 
  if (sample_method == "near") site_samples <- SiteSampleNear(sample_sites, coords, npts = npts)
  if (sample_method == "buffer") site_samples <- SiteSampleBuffer(sample_sites, coords, npts = npts, buffer_size = buffer_size)
  
  # make a vector
  samples <- paste0(site_samples$idx, "_", site_samples$site)
  
  return(samples)
}

# note: buffer was not used in the final analyses
SiteSampleBuffer <- function(sample_sites, coords, npts, buffer_size = 5){
  #sample_sites - coordinates of sampling sites
  #coords - coordinates of all individuals in data set
  #npts - number of points to sample from each site
  #buffer - buffer around site from which to draw samples randomly
  #edge_buffer - buffer from landscape edges to prevent sampling of sites
  
  # convert coords to sp
  coords <- sf::as_Spatial(coords)
  
  # make df
  site_samples <- data.frame()
  
  #sample sites can be provided as either sf or coords in vector format (nloop just counts how many sites)
  nloop <- nrow(sample_sites)
  
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
  
  #sample sites can be provided as either sf or coords in vector format (nloop just counts how many sites)
  nloop <- nrow(sample_sites)
  
  for(s in 1:nloop){
    # get site coords
    site_coords <- sample_sites[s,]
    # creates a vector of distances between the site and all other points in the dataset
    dist_vec <- as.vector(sf::st_distance(coords, site_coords))
    names(dist_vec) <- coords$idx
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


# randomly sample sites
rand_sitesamp <- function(coords, nsite, edge_buffer = NULL, ldim = 100){
  
  # buffer away from edges if ldim and edge_buffer provided
  # negatives are to deal with -y values in coords
  xmin = 0
  xmax = ldim
  ymin = -ldim
  ymax = 0
  if(is.null(ldim) | is.null(edge_buffer)){coords_buffer <- coords} else {coords_buffer <- sf::st_crop(coords, xmin = xmin + edge_buffer, ymin = ymin + edge_buffer, xmax = xmax - edge_buffer, ymax = ymax - edge_buffer)}
  
  #randomly select points to act as sites
  sample_sites <- coords_buffer[sample(1:nrow(coords_buffer), nsite),]
  
  return(sample_sites)
}


#make equidistant sampling sites
equi_sitesamp <- function(nsite, ldim = 100, buffer = 10){
  
  #nsite - number of points (or sites) to sample (should be a perfect square)
  #ldim - landscape dimension of one side (landscape should be a square)
  inc <- (ldim - buffer*2)/(sqrt(nsite) - 1)
  xgrid <- ygrid <- seq(0+buffer, ldim-buffer, inc) 
  cgrid <- expand.grid(xgrid, ygrid)
  colnames(cgrid) <- c("x","y")
  
  #flip y coordinates to match simulation coords
  cgrid$y <- -cgrid$y
  
  #par(pty = "s")
  #plot(cgrid, xlim = c(0,ldim), ylim = c(-ldim,0))
  
  cgrid <- sf::st_as_sf(cgrid, coords = c("x", "y"))
  
  return(cgrid)
}

#perform envgeo sitesampling
envgeo_sitesamp <- function(gsd_df, nsite, Nreps = 1000, edge_buffer = NULL, ldim = 100){
  
  sample.sets <- matrix(nrow=Nreps, ncol=nsite)
  results <- data.frame(env1.var=numeric(Nreps), env2.var=numeric(Nreps),
                        Mantel.r=numeric(Nreps), Mantel.p=numeric(Nreps),
                        mean.dist=numeric(Nreps))
  
  if(is.null(ldim) | is.null(edge_buffer)){
    gsd_df <- gsd_df
  } else {
    #define buffer
    xmin = 0 + edge_buffer
    xmax = ldim - edge_buffer
    ymin = -ldim + edge_buffer
    ymax = 0 - edge_buffer
    
    #buffer coordinates away from edge
    gsd_df <- gsd_df[gsd_df$x > xmin & gsd_df$x < xmax & gsd_df$y > ymin & gsd_df$y < ymax,]
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
  sample_sites <- gsd_df[as.character(sub_df$idx), c("x","y", "idx")]
  
  #convert to coordinates
  sample_sites <- sf::st_as_sf(sample_sites, coords = c("x","y"))
  
  return(sample_sites)
}

# GLOBAL PARAMETERS

#method to use for site sampling (can be "near" or "buffer")
site_method <- "near"
#number of individuals to sample per site
global_npts <- 10
#buffer away from edge for site selection (if using buffer method)
global_edge_buffer <- 5
#buffer size around sites from which individuals (pts) are sampled (if using buffer method)
global_buffer_size <- 5
