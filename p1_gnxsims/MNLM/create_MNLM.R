
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

combos <- expand.grid(seed = c(1, 2, 3), H = c(0.05, 0.5), r = c(0.3, 0.6))

mnlms <- pmap(combos, mnlm_create, .progress = TRUE)

walk(mnlms, mnlm_write)

mnlms <- pmap(combos, mnlm_read)

stats <- map(mnlms, mnlm_stats) %>% bind_rows()

mnlm_df <- 
  map(mnlms, ~as.data.frame(.x, xy = TRUE, ID = FALSE)) %>%
  map(~pivot_longer(.x, -c(x, y), names_to = "env")) %>%
  bind_rows() %>%
  left_join(stats) %>%
  mutate(name = paste0(env, " | M = ", M, " | r = ", r))

ggplot(mnlm_df) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  facet_wrap(~name, ncol = 4) +
  theme_void() +
  scale_fill_viridis_c(option = "magma") +
  theme(legend.position = "none")

