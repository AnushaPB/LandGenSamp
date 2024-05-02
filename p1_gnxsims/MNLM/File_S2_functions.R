# Create bootstrapped distribution of correlations and Moran's I
# x - raster stack of two layers
# ldim - subraster dimensions (creates a 100x100 subsampled raster)
# percent.na - maximum percentage of NA values allowed in the subsampled raster
# nsamp - how many subsampled rasters to generate
boot_rnM <- function(x,
                     ldim = 100,
                     nsamp = 100, 
                     percent.na = 0){
  
  #convert from pixels to coordinate units (assumes x and y units are the same)
  ldimconv <- ldim*res(x)[1]
  
  #Create buffered raster based on ldim (dimensions of subrasters)
  #so that when the subrasters are created in the next step they don't extend outside of the extent 
  #(only needed for xmax/ymax since the box is created by adding to an initial coordinate
  buffrast <- crop(x, ext(ext(x)[1],
                          ext(x)[2] - ldimconv,
                          ext(x)[3],
                          ext(x)[4] - ldimconv))
  
  #vectors for saving results
  r <- c()
  Mtavg <- c()
  Mprec <- c()
  
  #vector to save initial cell values so the same subraster isn't sampled twice
  initvec <- c()
  
  #get non na cell numbers
  notna <- which(complete.cases(values(buffrast)))
  
  #loop to create subrasters
  for (i in 1:nsamp){
    
    #Initial conditions for while loops
    pna <- 1
    tn <- 0
    init <- notna[1]
    
    #Loop through subrasters while the percent of NAs is above a threshold 
    #AND while the number of attempts to get the raster to meet this condition is below 1000 (just so it doesn't run forever)
    testn = 1000
    while (pna > percent.na & tn <= testn){
      #Loop through initial points for the subraster corner that have not been sampled before
      while (init %in% initvec){
        init <- sample(notna, size = 1)
      }
      
      #add point to vector so it isn't sampled again
      initvec <- c(initvec, init)
      #get coordinates of the initial corner (lower left)
      coord1 <- xyFromCell(buffrast, init)
      #add ldim to corner to get the upper right corner of the subraster
      coord2 <- coord1 + ldimconv
      #create coords object for cropping
      coords <- c(coord1[1], coord2[1], coord1[2], coord2[2])
      #crop raster with coords to get subraster
      newx <- crop(x, ext(coords))
      #calculate the percentage of NAs
      pna <- sum(is.na(values(newx)))/(ncell(newx)*nlyr(x))
      #tracker for number of loops
      tn <- tn + 1
    }
    if (tn == testn){"failed to find a subraster below the NA threshold in testn attempts"}
    #calculate and save pearson's correlation for each subraster
    #calc correlation
    lyrcor <- layerCor(newx, "pearson", na.rm=T)$`pearson`[1,2]
    
    Mt <- Moran(raster(newx[[1]]))
    Mp <- Moran(raster(newx[[2]]))
    
    r <- append(r, lyrcor)
    Mtavg <- append(Mtavg, Mt)
    Mprec <- append(Mprec, Mp)
  }
  
  return(list(r = r, Moran_prec = Mprec, Moran_tavg = Mtavg))
}


# Create bootstrapped distribution of Moran's I for tree canopy
# x - raster stack of one layer (tree canopy)
# ldim - subraster dimensions (creates a 100x100 subsampled raster)
# percent.na - maximum percentage of NA values allowed in the subsampled raster
# nsamp - how many subsampled rasters to generate
boot_tree <- function(x,
                      ldim = 100,
                      nsamp = 100, 
                      percent.na = 0){
  
  #convert from pixels to coordinate units (assumes x and y units are the same)
  ldimconv <- ldim*res(x)[1]
  
  #Create buffered raster based on ldim (dimensions of subrasters)
  #so that when the subrasters are created in the next step they don't extend outside of the extent 
  #(only needed for xmax/ymax since the box is created by adding to an initial coordinate
  buffrast <- crop(x, ext(ext(x)[1],
                          ext(x)[2] - ldimconv,
                          ext(x)[3],
                          ext(x)[4] - ldimconv))
  
  #vectors for saving results
  Mtree <- c()
  
  #vector to save initial cell values so the same subraster isn't sampled twice
  initvec <- c()
  
  #get non na cell numbers
  notna <- which(complete.cases(values(buffrast)))
  
  #loop to create subrasters
  for(i in 1:nsamp){
    
    #Initial conditions for while loops
    pna <- 1
    tn <- 0
    init <- notna[1]
    
    #Loop through subrasters while the percent of NAs is above a threshold 
    #AND while the number of attempts to get the raster to meet this condition is below 1000 (just so it doesn't run forever)
    testn = 1000
    while (pna > percent.na & tn <= testn){
      #Loop through initial points for the subraster corner that have not been sampled before
      while (init %in% initvec){
        init <- sample(notna, size = 1)
      }
      
      #add point to vector so it isn't sampled again
      initvec <- c(initvec, init)
      #get coordinates of the initial corner (lower left)
      coord1 <- xyFromCell(buffrast, init)
      #add ldim to corner to get the upper right corner of the subraster
      coord2 <- coord1 + ldimconv
      #create coords object for cropping
      coords <- c(coord1[1], coord2[1], coord1[2], coord2[2])
      #crop raster with coords to get subraster
      newx <- crop(x, ext(coords))
      #calculate the percentage of NAs
      pna <- sum(is.na(values(newx)))/(ncell(newx)*nlyr(x))
      #tracker for number of loops
      tn <- tn + 1
    }
    if (tn == testn){"failed to find a subraster below the NA threshold in testn attempts"}
    
    Mtc <- Moran(raster(newx))
    Mtree <- append(Mtree, Mtc)
  }
  return(list(Moran_tree = Mtree))
}