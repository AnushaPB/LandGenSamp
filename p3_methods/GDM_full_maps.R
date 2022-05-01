# CODE TO PLOT RGB MAPS FROM MODELS TO VISUALIZE GDM

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


#for scaling genetic distances from 0 to 1 for GDM
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


run_gdm_map <- function(gen, gsd_df, envlayers, distmeasure = "euc"){
  
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
    #Calculate PC distance based on  PCs (MODIFY to make based on % var explained)
    #use npcs based on sample size
    npcs <- round(nrow(gen)*0.5,0)
    pc_dist <- as.matrix(dist(pc$x[,1:npcs], diag = TRUE, upper = TRUE))
    gendist <- range01(pc_dist)
  } else if(distmeasure == "euc"){
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  } else {
    print("appropriate gen dist measure not specified, defaulting to euclidean")
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  }
  
  #Format gdm dataframe
  site <- 1:nrow(gendist) #vector of sites
  gdmGen <- cbind(site, gendist) #bind vector of sites with gen distances
  gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, env1 = gsd_df$env1, env2 = gsd_df$env2)
  gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
  
  #scale distance from 01
  #!THINK THIS THROUGH!
  gdmData$distance <- range01(gdmData$distance) 
  
  #run GDM
  gdm_model <- gdm(gdmData, geo = TRUE)
  
  if(is.null(gdm_model)){
    #create and return empty raster
    pcaRast <- raster(matrix(nrow=100,ncol=100))
    extent(pcaRast) <- c(0,100,0,100)
  } else {
    # CREATE MAP
    # Transform GIS layers
    rastTrans <- gdm.transform(gdm_model, envlayers)
    #remove na values
    rastDat <- na.omit(getValues(rastTrans))
    #run pca
    pcaSamp <- prcomp(rastDat)
    
    # make PCA raster
    nl <- nlayers(rastTrans)
    pcaRast <- predict(rastTrans, pcaSamp, index=1:nl)
    
    # scale rasters to get colors
    for(i in 1:nl){
      pcaRast[[i]] <- (pcaRast[[i]]-pcaRast[[i]]@data@min) /
        (pcaRast[[i]]@data@max-pcaRast[[i]]@data@min)*255
    }
    if(nl == 2){
      pcaRast <- stack(pcaRast, pcaRast[[2]]*0+255)
    } else if (nl == 1){
      pcaRast <- stack(pcaRast, pcaRast[[2]]*0+255, pcaRast[[2]]*0+225)
    }
  }
  
  
  
  #return raster
  return(pcaRast)
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


res_gdm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
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
    
    #get envlayers 
    #get path to dir with layers
    layer_path <- here(dirname(getwd()), "p1_gnxsims", "MNLM", "layers")
    #get path for each env layer
    layer_name1 <- paste0("seed", params[i,"seed"],"_env1_H",params[i,"H"]*100,"_r",params[i,"r"]*100,".csv")
    layer_name2 <- paste0("seed", params[i,"seed"],"_env2_H",params[i,"H"]*100,"_r",params[i,"r"]*100,".csv")
    #read in layers
    env1 <- read.csv(here(layer_path, layer_name1), header = FALSE)
    env2 <- read.csv(here(layer_path, layer_name1), header = FALSE)
    #convert to rasters
    env1 <- raster(as.matrix(env1))
    env2 <- raster(as.matrix(env2))
    #fix extents
    extent(env1) <- extent(env2) <- c(0,100,0,100)
    envlayers <- stack(env1,env2)
    
    #run model on full data set
    full_map <- run_gdm_map(gen_2k, gsd_df_2k, envlayers, distmeasure = "euc")
    plotRGB(full_map, r=1, g=2, b=3)
    
    #write full datafile (temp)
    file_name <- paste0("outputs/GDM_maps/gdm_map_",paramset,".tif")
    writeRaster(full_map, file_name)
    
    #save as png
    file_name <- paste0("outputs/GDM_maps/gdm_map_",paramset,".png")
    png(file_name)
    par(mfrow=c(1,1))
    plotRGB(full_map, r=1, g=2, b=3)
    dev.off()
    
    #save as png
    file_name <- paste0("outputs/GDM_maps/gdm_submaps_",paramset,".png")
    png(file_name)
    par(mfrow=c(length(npts),length(sampstrats)), mar = rep(1.5,4))
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_map <- run_gdm_map(subgen, subgsd_df, envlayers, distmeasure = "euc")
        plotRGB(sub_map, r=1, g=2, b=3, margins=TRUE, main=paste0(nsamp," | ",sampstrat))
        
      }
    }
    dev.off()
    
  }
  
  gc()
  
}