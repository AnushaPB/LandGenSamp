library("here") #paths
library("gdm") #GDM
library("viridis")
library("vcfR")

###########
#   GDM   #
###########

#CALCULATING PC BASED GENETIC DISTANCE
#convert gen data to matrix
gen <- as.matrix(gen[s,])
#perform PCA
pc <- prcomp(gen)
#Calculate PC distance based on first three PCs (?MODIFY?)
pc_dist <- as.matrix(dist(pc$x[,1:100], diag = TRUE, upper = TRUE))

#Format gdm dataframe
site <- 1:nrow(pc_dist) #vector of sites
gdmGen <- cbind(site, pc_dist) #bind vector of sites with gen distances
gdmPred <- data.frame(site = site, Longitude = gea_df$x, Latitude = gea_df$y, env1 = gea_df$env1, env2 = gea_df$env2)
gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
#SCALE DISTANCE FROM 0 to 1 if max(distance) >1 (gdm only works for 0<vals<1) (?MODIFY?)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
gdmData$distance <- range01(gdmData$distance) 

#run GDM
gdm.model <- gdm(gdmData, geo = TRUE)

# Sum coefficients for each predictor (each has 3 splines)
coeffs <- function(gdm.model){
  coefSums <- c()
  for (i in 1:length(gdm.model$predictors)){
    j <- (i * 3) - 2
    coefSums[i] <- sum(gdm.model$coefficients[j:(j+2)])
  }
  
  # Add those values to a simple data frame
  coeffs <- data.frame(predictor = gdm.model$predictors, coefficient = coefSums)
  return(coeffs)
}

predictors <- coeffs(gdm.model)
predictors

env <- read.csv("data/seed1_env1_H50_r60.csv", header = F)
#env <- read.csv("data/env_nlm1.csv", header = F)
envM <- as.matrix(env)
r <- raster(envM)
extent(r) <- c(0,41,0,41) #change extent to match coords
env1 <- flip(r, direction='y') #idk why but you have to flip the raster to match the coords

env <- read.csv("data/seed1_env2_H50_r60.csv", header = F)
#env <- read.csv("data/env_nlm1.csv", header = F)
envM <- as.matrix(env)
r <- raster(envM)
extent(r) <- c(0,41,0,41) #change extent to match coords
env2 <- flip(r, direction='y') #idk why but you have to flip the raster to match the coords


# Transform GIS layers
rastTrans <- gdm.transform(gdm.model, stack(env1,env2))
rastDat <- na.omit(getValues(rastTrans))
pcaSamp <- prcomp(rastDat)

# note the use of the 'index' argument
pcaRast <- predict(rastTrans, pcaSamp, index=1:3)

# scale rasters
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255

par(mfrow=c(3,2), mar = rep(0,4), oma = rep(0,4), pty="s")
plotRGB(pcaRast, r=1, g=2, b=3)
plotRGB(pcaRast, r=1, g=3, b=2)
plotRGB(pcaRast, r=2, g=3, b=1)
plotRGB(pcaRast, r=2, g=1, b=3)
plotRGB(pcaRast, r=3, g=2, b=1)
plotRGB(pcaRast, r=3, g=1, b=2)

par(mfrow=c(1,3), mar = rep(1,4), oma = rep(0,4), pty="s")
plot(env1, axes=FALSE, box=FALSE, col=magma(100), legend=FALSE, main = "env1")
plot(env2, axes=FALSE, box=FALSE, col=magma(100), legend=FALSE, main = "env2")
plotRGB(pcaRast, r=1, g=3, b=2, main = "GDM")

