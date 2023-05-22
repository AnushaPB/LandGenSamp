library(NLMR)
library(raster)

mnlm <- function(ncol, nrow, H, sd = 1, r){
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
  Y <- m.2 %*% diag(1/sqrt(colSums(m.2^2)))  # scale columns to length 1
  
  nlm.2.v <- Y[ , 2] + (1/tan(theta)) * Y[ , 1]  # Vector of new values
  nlm.2.v <- (nlm.2.v - min(nlm.2.v)) / (max(nlm.2.v) - min(nlm.2.v))  # Scale new values from 0 to 1
  nlm.2.m <- matrix(data = nlm.2.v, nrow = nrow, ncol = ncol, byrow = TRUE)
  nlm.2 <- raster(nlm.2.m)
  
  extent(nlm.2) <- extent(nlm.1)  # Match extents
  s <- stack(nlm.1, nlm.2)
  return(s)
}

nlm.s <- mnlm(51, 51, H = 0.9, sd = 1, r = 0.6)

layerStats(nlm.s, stat = 'pearson')
Moran(nlm.s[[1]])
Moran(nlm.s[[2]])
plot(nlm.s)
