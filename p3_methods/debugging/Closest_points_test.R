# 
sample_site <- gsd_df[2,"idx"]
nsamp = 10

# get site coords
site_coords <- gsd_df[as.character(sample_site), c("x","y")]

# creates a vector of distances between the site and all other points in the dataset
dist_vec <- sqrt((gsd_df$x - site_coords$x)^2 + (gsd_df$y - site_coords$y)^2)

names(dist_vec) <- gsd_df$idx


# remove the distance of 0 (same site)
dist_vec <- dist_vec[dist_vec != 0]
dist_vec <- dist_vec[order(dist_vec)]
idx <- names(dist_vec)[1:npts]

# order distances
# find the 5 nearest stations and create new data frame
near <- dist_vec %>%  
  data.frame() %>%
  gather('idx','dist') %>% 
  filter(!is.na(dist)) %>% 
  arrange(dist) %>% 
  slice(1:nsamp)

# return IDs
return(as.character(near$idx))
