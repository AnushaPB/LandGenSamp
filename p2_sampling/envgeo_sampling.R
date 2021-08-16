source("general_functions.R")
library("here")
library("foreach")
library("doParallel")
library("vegan")

envgeo_samp <- function(pts, npts, Nreps = 1000){
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

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

for(n in npts){
  samples <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    library("here")
    library("vegan")
    
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
      samples <- envgeo_samp(pts, npts = n, Nreps = 1000)
    }
    
    #return vector of sample IDs
    return(samples)
    
  }
  
  #bind sample IDs together and export (rows are parameter sets/columns are individual IDs)
  colnames(samples) <- paste0("envgeo",1:ncol(samples))
  samp_out <- params
  for(i in 1:ncol(samples)){samp_out <- cbind.data.frame(samp_out, samples[,i])}
  colnames(samp_out) <- c(colnames(params),colnames(samples))
  write.csv(samp_out, paste0("outputs/samples_envgeo",n,".csv"), row.names = FALSE)
}


#stop cluster
stopCluster(cl)
  


