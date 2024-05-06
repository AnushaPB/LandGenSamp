
#' Multiple neutral landscape models (fractal Brownian motion)
#'
#' @param nlayers The number of NLMs to generate (numeriCal)
#' @param r The correlation coefficient between the first NLM and each successive NLM (numerical or vector)
#' @param nrow The number of rows in the rasters (numerical)
#' @param ncol The number of columns in the rasters (numerical)
#' @param reslution The resolution of the rasters (numerical)
#' @param roughness FILL
#' @param rand_dev FILL
#' @param torus FILL
#' @param user_seed Set seed for simulation (numerical)
#' @param rescale If TRUE (default), raster values are scaled from 0 to 1 (logical)
#' @param ... character, logical, or numeric (optional). Additional arguments to be passed to RandomFields::RFoptions (n.b. if using a fractal dimension between ~ 1.6 and 1.9, one must set the option modus_operandi = "sloppy".
#' @details
#' Generates multiple neutral landscape models using fractional Brownian motion, an extension of Brownian motion in which the correlation between steps is controlled by frac_dim.
#' Higher values of frac_dim produce smoother, more correlated surfaces, while lower value produce rougher, less correlated surfaces.
#' The r argument can accept either a single value, in which case all NLMs produced will have the same correlation with the first layer, or a vector containing the desired correlation coefficients for each layer.
#' @examples
#' NLMs <- mnlm_mpd(nlayers = 3, r = c(0.3, 0.6), ncol = 20, nrow = 20)
#' ## Not run:
#' layerStats(NLMs, stat = "pearson")
#' ## End(**Not run**)
#' @importFrom NLMR nlm_mpd
#' @importFrom raster stack layerStats
#' @export
mnlm_mpd <- function(nlayers = 2, r, ncol, nrow, resolution = 1, roughness = 0.5, rand_dev = 1, torus = FALSE, user_seed = NULL, rescale = TRUE, ...){
  if(length(r) == 1) r <- rep(r, nlayers - 1)

  nlm.s <- stack()
  for(i in 1:nlayers){
    nlm <- nlm_mpd(ncol = ncol, nrow = nrow, resolution = resolution, rescale = rescale, roughness = roughness, rand_dev = rand_dev, torus = torus)
    nlm.s <- stack(nlm.s, nlm)
  }

  for(j in 2:nlayers){
    s <- stack(nlm.s[[1]], nlm.s[[j]])
    newNLM <- rasterQR(s, r[j-1])
    extent(newNLM) <- extent(nlm.s[[1]]) # Match extents
    nlm.s[[j]] <- newNLM
  }

  return(nlm.s)

}
#' Raster QR decomposition
#'
#' @param s A RasterStack
#' @param r A correlation coefficient (numeric)
#' @importFrom raster getValues raster
#' @export
rasterQR <- function(s, r){
  # Extract values and bind
  vals.1 <- getValues(s[[1]])
  vals.2 <- getValues(s[[2]])
  m <- cbind(vals.1, vals.2)
  m.ctr  <- scale(m, center = TRUE, scale = FALSE)
  
  # Prepare matrix projection
  Id.m <- diag(length(m.ctr[, 1]))  # Identity matrix
  Q <- qr.Q(qr(m.ctr[ , 1, drop = FALSE]))  # QR decomposition
  P <- tcrossprod(Q)  # Cross-product of matrix for projection
  xvals.o <- (Id.m - P) %*% m.ctr[ , 2]  # Make xvals orthogonal to m.ctr
  m.2  <- cbind(m.ctr[ , 1], xvals.o)
  Y <- m.2 %*% diag(1 / sqrt(colSums(m.2^2)))  # Scale columns to length 1
  
  # Generate second NLM raster
  theta <- acos(r) # Calculate corresponding angle for desired correlation
  newNLM.v <- Y[ , 2] + (1/tan(theta)) * Y[ , 1]  # Vector of new values
  newNLM.v <- (newNLM.v - min(newNLM.v)) / (max(newNLM.v) - min(newNLM.v))  # Scale new values from 0 to 1
  newNLM.m <- matrix(data = newNLM.v, nrow = nrow(s[[2]]), ncol = ncol(s[[2]]), byrow = TRUE)
  newNLM <- raster(newNLM.m)
  return(newNLM)
}

# Create MNLMs with a fixed H and r
# Iteratively find a pair of MNLMs where the difference in autocorrelation is <0.01
mnlm_create <- function(seed, H, r){
  set.seed(seed)
  
  # Find a pair of MNLMs where the difference in autocorrelation is <0.01
  Hdif = TRUE
  it = 0
  while (Hdif & it < 100){
    # note: 103x103 produces a 101x101 raster
    NLMs <- mnlm_mpd(nlayers = 2, r = r, ncol = 103, nrow = 103, roughness = 1 - H)
    
    # Crop to desired extent
    NLM1 <- crop(NLMs[[1]], c(0, 100, 0, 100))
    NLM2 <- crop(NLMs[[2]], c(0, 100, 0, 100))
    NLMs <- stack(NLM1, NLM2)
    
    # Calculate Moran's I
    MI1 <- Moran(NLMs[[1]])
    MI2 <- Moran(NLMs[[2]])
    
    # Determine difference in Moran's I 
    Hdif <- abs(MI1 - MI2) 
    
    print(paste("it:", it, "| H dif:", round(Hdif, 3)))
    
    # Continue if difference is greater than 0.01 and it < 100
    Hdif <- Hdif > 0.01
    it <- it + 1
    
    if (it == 100) warning("No landscapes found")
  }
  
  # Name NLMs
  names(NLMs) <- paste0("seed", seed, "_", c("env1", "env2"), "_H", gsub("\\.", "", as.character(H*100)), "_r", gsub("\\.", "", as.character(r*100)))
  return(NLMs)
}

# Write out MNLMs
mnlm_write <- function(x){
  path <- here("p1_gnxsims", "MNLM", "layers")
  rls <- as.list(x)
  
  walk(rls, ~write.table(
    as.matrix(.x),
    row.names = FALSE, 
    col.names = FALSE, 
    sep = ",", 
    here(path, paste0(names(.x), ".csv"))
    ))
}

# Get all MNLMs
mnlm_get <- function(){
  combos <- expand.grid(seed = c(1, 2, 3), H = c(0.05, 0.5), r = c(0.3, 0.6))
  pmap(combos, mnlm_read)
}


# Read in a pair of MNLMs
mnlm_read <- function(seed, H, r){
  folder_path <- here("p1_gnxsims", "MNLM", "layers")
  file_paths <- paste0("seed", seed, "_", c("env1", "env2"), "_H", gsub("\\.", "", as.character(H*100)), "_r", gsub("\\.", "", as.character(r*100)), ".csv")
  env <- 
    map(file_paths, ~raster(as.matrix(read.csv(here(folder_path, .x), header = FALSE)))) %>%
    stack() 
  extent(env) <- c(0, 100, 0, 100)
  names(env) <- paste0("seed", seed, "_", c("env1", "env2"), "_H", gsub("\\.", "", as.character(H*100)), "_r", gsub("\\.", "", as.character(r*100)))
  return(env)
}

# Calculate Moran's I and raster correlation for MLNMs
mnlm_stats <- function(env){
  MI1 <- round(Moran(env[[1]]), 2)
  MI2 <- round(Moran(env[[2]]), 2)
  r <- round(layerStats(env, "pearson")$`pearson correlation coefficient`, 1)[1, 2]
  data.frame(env = names(env), M = c(MI1, MI2), r = r)
}
