
library(raster)
library(here)
#install.packages("devtools")
#devtools::install_github("ropensci/NLMR")
library(NLMR)
#options(repos = c(
#  jeffreyevans = 'https://jeffreyevans.r-universe.dev',
#  CRAN = 'https://cloud.r-project.org'))
# Download and install spatialEco in R
# install.packages('spatialEco')
library(spatialEco)
library(viridis)
library(terra)

wdir <- here("p1_gnxsims", "MNLM")

# Function to generate MNLMs:
mnlm <- function(ncol, nrow, H, sd = 1, r){
  if(ncol %% 2 != 0){
    ncol <- ncol+2
    nrow <- nrow+2
  }
  
  theta <- acos(r) # Corresponding angle for desired correlation
  
  # Generate first NLM
  nlm.1 <- nlm_mpd(ncol, nrow, resolution = 1, roughness = 1-H, rand_dev = sd, rescale = TRUE, verbose = TRUE)
  nlm.1.v <- getValues(nlm.1)
  
  # Generate values for generating second NLM with same spatial structure
  xvals.r <- nlm_mpd(ncol, nrow, resolution = 1, roughness = 1-H, rand_dev = sd, rescale = TRUE, verbose = TRUE)
  xvals.v <- getValues(xvals.r)
  m <- cbind(nlm.1.v, xvals.v)
  m.ctr  <- scale(m, center = TRUE, scale = FALSE)   
  
  # Prepare matrix projection
  Id.m <- diag(length(nlm.1.v))  # Identity matrix
  Q <- qr.Q(qr(m.ctr[ , 1, drop = FALSE]))  # QR decomposition
  P <- tcrossprod(Q)  # Cross-product of matrix for projection
  xvals.o <- (Id.m-P) %*% m.ctr[ , 2]  # Make xvals orthogonal to m.ctr
  m.2  <- cbind(m.ctr[ , 1], xvals.o)
  Y <- m.2 %*% diag(1/sqrt(colSums(m.2^2)))  # Scale columns to length 1
  
  # Generate second NLM
  nlm.2.v <- Y[ , 2] + (1/tan(theta)) * Y[ , 1]  # Vector of new values
  nlm.2.v <- (nlm.2.v - min(nlm.2.v)) / (max(nlm.2.v) - min(nlm.2.v))  # Scale new values from 0 to 1
  nlm.2.m <- matrix(data = nlm.2.v, nrow = nrow(xvals.r), ncol = ncol(xvals.r), byrow = TRUE)
  nlm.2 <- raster(nlm.2.m)
  
  # Match extents and stack NLMs
  extent(nlm.2) <- extent(nlm.1)  
  s <- stack(nlm.1, nlm.2)
  return(s)
}

f <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)

dim <- 100
n <- 3
d <- (n^2-1)/2
f <- matrix(c(rep(1,d),0,rep(1,d)), n, n)

par(mfrow=c(4,3), mar=rep(2,4), oma=rep(1,4))
for(seed in c(1, 2, 3)){
  for(H in c(0.05,0.5)){
    for(r in c(0.3, 0.6)){
      set.seed(seed)
      
      nlm.s <- mnlm(101, 101, H = H, sd = 1, r = r)
      
      env1 <- nlm.s[[1]]
      #crop from 101x101 to 100x100 so that it is easier to partition the landscape for grid sampling later
      env1 <- crop(env1, c(0,100,0,100))
      MI1 <- round(Moran(env1, w=f),2)
      
      env2 <- nlm.s[[2]]
      env2 <- crop(env2, c(0,100,0,100))
      MI2 <- round(Moran(env2,  w=f),2)
      
      pearcor <- round(layerStats(stack(env1,env2),"pearson")$`pearson correlation coefficient`,1)[1,2]
      
      plot(env1, box=FALSE, axes=FALSE, col=inferno(100), legend.width = 2, main = paste0("env1 (Moran's I = ",MI1,")"))
      plot(env2, box=FALSE, axes=FALSE, col=inferno(100), legend.width = 2, main = paste0("env2 (Moran's I = ",MI2,")"))
      
      rascor <- rasterCorrelation(rast(env1),rast(env2),5)
      MIr <- round(Moran(raster(rascor)),2)
      plot(rascor, col=inferno(100), box=FALSE, axes=FALSE, zlim = c(-1,1), 
           main = paste0("H = ",H," | r = ",pearcor,"\n(Moran's I = ",MIr,")"), legend.width=2)
      
      write.table(as.matrix(env1), row.names = FALSE, col.names=FALSE, sep=",",
                  here(wdir, "layers", paste0("seed",seed,
                                              "_env1",
                                              "_H", gsub("\\.", "", as.character(H*100)),
                                              "_r", gsub("\\.", "", as.character(r*100)),
                                              ".csv"
                  )))
      
      write.table(as.matrix(env2), row.names = FALSE, col.names=FALSE, sep=",",
                  here(wdir, "layers", paste0("seed",seed,
                                              "_env2",
                                              "_H", gsub("\\.", "", as.character(H*100)),
                                              "_r", gsub("\\.", "", as.character(r*100)),
                                              ".csv"
                  )))
      
    }
  }
}
