# CODE TO CALCULATE AND COMPARE DIFFERENT GENETIC DISTANCE MEASURES

set.seed(42)

library("here") #paths
library("gdm") #GDM
library("vcfR")
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")


###########
#   GDM   #
###########


gendist <- function(gen, distmeasure = "dps"){
  
  if(distmeasure == "bray"){
    K <- nrow(gen)
    nloc <- ncol(gen)
    ret <- matrix(0,K,K)
    rownames(ret) <- colnames(ret) <- rownames(gen)
    for( i in 1:K){
      for( j in 1:i){
        if( i != j){
          ret[i,j] <- ret[j,i] <- sum(apply(gen[ c(i,j), ], 2, min)) / nloc
        }
      }
    }
    gendist <- ret
  } else if(distmeasure == "dps"){
    #DPS GENETIC DISTANCE
    gen[gen == 0] <- "11"
    gen[gen == 1] <- "12"
    gen[gen == 2] <- "22"
    
    genindobj <- df2genind(gen, ploidy=2, ncode=1)
    psh <- propShared(genindobj)
    dps <- 1 - psh
    gendist <- dps
  } else if(distmeasure == "pca"){
    #perform PCA
    pc <- prcomp(gen)
    #Calculate PC distance based on  PCs (?MODIFY?)
    #use npcs based on sample size
    npcs <- round(nrow(gen)*0.5,0)
    
    pc_dist <- as.matrix(dist(pc$x[,1:npcs], diag = TRUE, upper = TRUE))
    gendist <- pc_dist
  } else if(distmeasure == "euc"){
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  } else {
    print("appropriate gen dist measure not specified")
  }
  
  return(gendist)
}


range01 <- function(x){(x-min(x))/(max(x)-min(x))}


plot_dist <- function(gen_dist, xlim = c(0,1), range01 = TRUE, xlab = "Gen Dist", main = ""){
  if(range01){gen_dist <- range01(gen_dist)}
  d <- density(gen_dist)
  plot(d, xlab = xlab, main = main, col = "#CCCCFF", lwd=2)
  polygon(d, col=rgb(	204, 204, 255, maxColorValue = 255, alpha = 200), border = "#7c83bc", lwd = 2)
}


#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


for(i in 1:nrow(params)){
  #vcfR
  library("vcfR")
  library("gdm")
  library("adegenet")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  skip_to_next <- FALSE
  if(file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run GDM
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    par(mfrow=c(1,1))
    #get gen dist for full data set
    full_result_dps <- gen_dist(gen_2k, distmeasure = "dps")
    plot_dist(full_result_dps, xlab = "dps", main = "2000 | full", range01 = FALSE)
    #full_result_pca <- gen_dist(gen_2k, distmeasure = "pca")
    #plot_dist(full_result_dps, xlab = "pca", main = "2000 | full")
    #env
    #plot_dist(dist(gsd_df[,c("env1","env2")]), xlab = "envdist", main = "2000 | full")
    
    
    par(mfrow=c(4,4))
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result_dps <- gen_dist(subgen, distmeasure = "dps")
        plot_dist(sub_result_dps, xlab = "dps", main = paste0(nsamp," | ", sampstrat), range01 = FALSE)
        #sub_result_pca <- gen_dist(subgen, distmeasure = "pca")
        #plot_dist(sub_result_dps, xlab = "pca", main = paste0(nsamp," | ", sampstrat))
        
        #plot_dist(dist(subgsd_df[,c("env1","env2")]), xlab = "envdist",  main = paste0(nsamp," | ", sampstrat))
        
      }
    }
    
    par(mfrow=c(1,1))
    plot_dist(dist(gsd_df[,c("env1","env2")]), xlab = "envdist", main = "2000 | full")
    
    par(mfrow=c(4,4))
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        plot_dist(dist(subgsd_df[,c("env1","env2")]), xlab = "envdist",  main = paste0(nsamp," | ", sampstrat))
        
      }
    }
    
    
    par(mfrow=c(1,1))
    plot_dist(dist(gsd_df[,c("x","y")]), xlab = "geodist", main = "2000 | full")
    
    par(mfrow=c(4,4))
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        plot_dist(dist(subgsd_df[,c("x","y")]), xlab = "geodist",  main = paste0(nsamp," | ", sampstrat))
        
      }
    }
    
  }
  
  return(result)
  
  gc()
  
}

#stop cluster
stopCluster(cl)

write.csv(res_gdm, "outputs/gdm_results_dps.csv", row.names = FALSE)

