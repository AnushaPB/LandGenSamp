
#MODIFIED FOR SITE SAMPLING - neaten this up later
#get list of sampling IDs that correspond with parameter set, sampling strategy, and number of samples
get_samples <- function(param_set, params = params, sampstrat, nsamp, outdir = here(dirname(getwd()), "p2_sampling", "outputs")){
  #param_set - vector of one set of parameters (e.g. params[i,])
  #sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  #nsamp - number of samples
  
  #Check if files for parameter exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  print(gen_filepath)
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  print(gsd_filepath)
  loci_filepath <- create_filepath(i, params = params, "loci")
  print(loci_filepath)
  file_exists <- TRUE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){file_exists <- FALSE}
  if(!file_exists) { 
    print("File does not exist:")
    print(params[i,]) 
  } 
  stopifnot(file_exists)
  
  #TO DO - only thing changed is "samples" to "site_samples" so fix this later to make it cleaner
  subIDs <- read.csv(paste0(outdir, "/site_samples_", sampstrat, nsamp, ".csv"))
  
  subIDs <- subIDs[subIDs$K == param_set$K 
                   & subIDs$phi == param_set$phi
                   & subIDs$m == param_set$m 
                   & subIDs$seed == param_set$seed
                   & subIDs$H == param_set$H
                   & subIDs$r == param_set$r
                   & subIDs$it == param_set$it,]
  
  #confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)
  
  #remove parameter columnds and convert to vector of IDs
  subIDs <- subIDs[,!names(subIDs) %in% colnames(params)]
  subIDs <- unlist(subIDs)
  
  #confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))
  
  return(as.character(subIDs))
}


SiteSample <- function(sample_sites, coords, npts, buffer_size = 300000){
  #sample_sites - coordinates of sampling sites
  #coords - coordinates of all individuals in data set
  #npts - number of points to sample from each site
  #buffer - buffer around site from which to draw samples randomly
  #edge_buffer - buffer from landscape edges to prevent sampling of sites
  site_samples <- data.frame()
  
  #if sp
  if(class(sample_sites)[1] == "SpatialPoints"){nloop <- length(sample_sites)} else {nloop <- nrow(sample_sites)}
  
  for(s in 1:nloop){
    #create buffer around sites from which to sample points
    site_buffers <- buffer(sample_sites[s,], buffer_size)
    #subset out only coordinates falling within site buffers
    buffer_samples <- coords[site_buffers,]
    #convert from SPDF to df
    buffer_samples_df <- data.frame(buffer_samples)
    #randomly sample samples from coordinates
    randsamp <- sample(nrow(buffer_samples_df), npts)
    #subset df
    buffer_samples_df <- buffer_samples_df[randsamp,]
    #IMPORTANT: remove sample from coords so they are not sampled twice
    coords <- coords[-c(coords$idx %in% buffer_samples_df$idx),]
    
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

nsites <- c(9, 16, 25)
sampstrats <- c("rand", "equi", "envgeo")
npts <- 10

