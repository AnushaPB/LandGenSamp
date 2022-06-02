
source("general_functions.R")
source("site_functions.R")
library("here")
library("foreach")
library("doParallel")
library("vegan")

set.seed(42)

envgeo_samp <- function(gsd_df, nsite, Nreps = 1000, edge_buffer, ldim){
  Nreps <- 1000
  sample.sets <- matrix(nrow=Nreps, ncol=nsite)
  results <- data.frame(env1.var=numeric(Nreps), env2.var=numeric(Nreps),
                        Mantel.r=numeric(Nreps), Mantel.p=numeric(Nreps),
                        mean.dist=numeric(Nreps))
  #define buffer
  buffmin <- edge_buffer
  buffmax <- ldim - edge_buffer
  #buffer coordinates away from edge
  gsd_df <- gsd_df[gsd_df$x > buffmin & gsd_df$x < buffmax & gsd_df$y > buffmin & gsd_df$y < buffmax,]
  
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
  
  
  #save IDs to vector
  samples <- as.character(sub_df$idx)
  
  return(samples)
}

#register cores
#these calculations are RAM intensive so only run two at a time
cores <- 2
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

for(n in nsites){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("vegan")
    library("raster")
    library("rgeos")
    
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
      #note - add an edge buffer from ldim so sites aren't sampled close to the edge
      sample_sites <- envgeo_samp(gsd_df, nsite = n, Nreps = 1000, edge_buffer = global_edge_buffer, ldim = ldim)
      #overwrite sample sites with coordinates for sample sites using indexes
      sample_sites <- gsd_df[sample_sites, c("x","y")]
      #convert to coordinates
      coordinates(sample_sites) <- ~x+y
      
      #sample from around sites based on a buffer
      #buffer size chosen arbitrarily, other size was too small/not enough points in buffer for smaller sample sizes
      site_samples <- SiteSample(sample_sites, coords, npts = global_npts, buffer_size = global_buffer_size)
      
      #plot (for debugging)
      plot(sample_sites, xlim = c(0,ldim), ylim = c(0,ldim))
      points(gsd_df[,c("x","y")], col = "gray")
      points(site_samples[,c("x","y")], col = "red")
      #points(site_samples[,c("xsite","ysite")], col = "blue", pch = 19)
      
      samples <- paste0(site_samples$idx, "_", site_samples$site)
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("envgeo", 1:ncol(samples))
  #sample IDs
  sampleIDs <- gsub("\\_.*","",samples)
  #site IDs
  siteIDs <- gsub("^.*\\_","", samples)
  
  #create df of sample IDs
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, sampleIDs[,i])}
  colnames(samp_out) <- c(colnames(params), colnames(samples))
  write.csv(samp_out, paste0("outputs/site_samples_envgeo",n,".csv"), row.names = FALSE)
  
  #create df of site IDs
  site_out <- params
  for(i in 1:ncol(samples)){site_out <- cbind.data.frame(site_out, siteIDs[,i])}
  colnames(site_out) <- c(colnames(params), colnames(samples))
  write.csv(site_out, paste0("outputs/site_ids_envgeo",n,".csv"), row.names = FALSE)
  
}


#stop cluster
stopCluster(cl)

