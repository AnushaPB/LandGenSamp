#Transect Sampling

transect_samp <- function(pts, npts, ytsct, buffer){
  #pts - dataframe with IDs and coords
  #npts - total number of points to sample (evenly split across transects)
  #buffer - buffer around transects within which points are sampled 
  
  #divide number of samples evenly among the transects
  npts_tsct <- npts/length(ytsct)
    
  #plot all points (gray) (for debugging, comment out later)
  par(pty="s")
  plot(gsd_df$x, gsd_df$y, pch=19, cex=0.2, col="gray", main = npts)
    
  #create empty vector to store IDs
  IDvec <- c()
  for(i in 1:length(ytsct)){ 
    #subset points around transect based on buffer
    tsctsq <- subset(pts, y > (ytsct[i] - buffer) & y < (ytsct[i]+buffer))
    #randomly sample subset of transect points to match number of samples needed for each transect
    tsctsq <- tsctsq[sample(nrow(tsctsq), npts_tsct),]
    #plot points sampled (for debugging, comment out later)
    points(tsctsq$x, tsctsq$y, col=i+1)
    #store IDs in list
    IDvec <- c(IDvec, tsctsq$idx)
  }
  
  #confirm correct number of samples were subsetted
  stopifnot(npts == length(IDvec))
  
  #return vec of sample IDs
  return(IDvec)
}

#data (EXAMPLE)
gsd_df <- read.csv("mod-10k_K1_phi50_m100_seed1_H50_r60_it--1_t-0_spp-spp_0.csv")
#subset out necessary columns from df
pts <- gsd_df[,c("idx","x","y")]

#define vec of number of sample points
npts_vec <- c(36, 81, 144, 225, 324)
#define horizontal transects (y-coords)
ytsct <- c(10, 20, 30)
#define buffer around transects
buffer <- 3

#plot of transects and buffer
par(pty="s")
plot(pts$x, pts$y, pch=19, col="gray", cex=0.2)
abline(h = ytsct, col="black", lty="dashed", lwd=2)
abline(h = ytsct + buffer, col="red", lty="dashed")
abline(h = ytsct - buffer, col="red", lty="dashed")

par(pty="s", mfrow = c(2,3))
#create list to store IDs of samples
subs <- vector(mode = "list", length = length(npts_vec))
for(npts in npts_vec){
  subs[[n]] <- transect_samp(pts, npts, ytsct, buffer)
}

