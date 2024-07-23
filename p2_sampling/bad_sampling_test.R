library(here)
library(tidyverse)
library(raster)
library(terra)
library(sf)
source(here("p1_gnxsims", "MNLM", "MNLM_functions.R"))

# FAKE GSD DF
gsd_df <- data.frame(x = sample(seq(0,100,0.01), 8000, replace = TRUE),
                     y = sample(seq(0,100,0.01), 8000, replace = TRUE),
                     ID = 1:1000)

coords <- st_as_sf(gsd_df, coords = c("x", "y"))

# Load the landscapes
nlms <- mnlm_get()

# Convert to terra rast
nlms <- map(nlms, rast)

# Mask out extremes from raster and coords
nlms_mask <- 
  nlms %>%
  map(~{
    mask_values <- 0.125 * 2
    mask <- .x < mask_values | .x > (1 - mask_values)
    .x[mask] <- NA
    .x
  })

# Get values from landscapes for coords to make a set of fake gsds
fake_gsds <- 
  nlms %>% 
  map(~extract(.x, coords, xy = TRUE)) %>%
  map(~{
    names(.x) <- c("ID", "env1", "env2", "x", "y")
    return(.x)
  })

bad_sampling <- function(x, nsamp){
  # Remove extreme environmental values
  extremes <- 0.125 * 2
  mask_coords <- 
    x %>%
    filter(env1 > extremes, env1 < 1 -extremes, env2 > extremes, env2 < 1 - extremes) %>%
    st_as_sf(coords = c("x", "y"))

  # Sample 3 points from the masked coordinates
  # buffer by 10 coordinate units so sample isn't on the edge
  buffer_coords <-
   mask_coords %>% 
   st_crop(ext(10, 90, 10, 90))
  # Check: ggplot() + geom_sf(data = mask_coords, col = "black") + geom_sf(data = buffer_coords, col = "red")
  s <- sample(1:nrow(buffer_coords), 3, replace = FALSE)
  init_coords <- buffer_coords[s,]

  # Sample additional n points based on an inverse distance weighted probability distribution
  idw_coords <- map(1:nrow(init_coords), ~idw_sampling(init_coords[.x,], mask_coords, nsamp = nsamp)) %>% bind_rows()

  return(idw_coords)
}

# Sample additional 10 points from each landscape based on an inverse distance weighted probability distribution
idw_sampling <- function(init_sample, coords, nsamp){
  # Calculate the distances from the given coordinate to the other coordinates
  distances <- as.vector(sf::st_distance(init_sample, coords, method = "euclidean"))

  # Calculate the inverse distances
  inverse_distances <- 1 / distances

  # Transform
  inverse_distances_2 <- inverse_distances ^ 2

  # Normalize the inverse distances to get probabilities
  # Remove Inf values (Inf = distance of 0 = the original coordinate/any overlapping coordinates)
  probabilities <- inverse_distances_2 / sum(inverse_distances_2[!is.infinite(inverse_distances_2)], na.rm = TRUE)

  # Give the Inf inverse distance value the max probability value
  # This increase the probability that the original site is included
  probabilities[is.infinite(probabilities)] <- max(probabilities[!is.infinite(probabilities)], na.rm = TRUE)

  # Check: 
  #coords$probabilities <-as.vector(probabilities)
  #ggplot() +
  #  geom_sf(data = coords, aes(color = probabilities)) +
  #  geom_sf(data = init_sample, color = "red") +
  # theme_void()

  # Sample coordinates based on the probabilities
  sampled_coords <- sample(1:nrow(coords), size = nsamp/3, prob = probabilities, replace = FALSE)

  # Get the sampled coordinates
  sampled_coords <- coords[sampled_coords, ]

  return(sampled_coords)
}

example_coords <- map(c(36, 81, 144, 225), ~map(fake_gsds, bad_sampling, nsamp = .x))

plts <- 
  map(example_coords, ~{
    map2(.x, nlms, ~{
      ggdf <- terra::as.data.frame(.y, xy = TRUE) %>% pivot_longer(-c(x, y))
      ggplot() +
        geom_raster(data = ggdf, aes(x = x, y = y, fill = value)) +
        geom_sf(data = .x, color = "tomato") +
        facet_wrap(~name) +
        theme_void() +
        scale_fill_viridis_c(option = "mako") +
        ggtitle(paste("n = ", nrow(.x)))
    })
  })
  

ge <- map(plts, ~do.call(gridExtra::grid.arrange, .x[c(1,4)]))
do.call(gridExtra::grid.arrange, ge)
