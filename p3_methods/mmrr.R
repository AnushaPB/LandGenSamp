set.seed(42)

# Parallel processing libraries
library(furrr)

# Read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")
source("mmrr_functions.R")

perform_mmrr_analysis <- function(params, ns, strats, full, site = FALSE) {
  results <- future_map(1:nrow(params), \(i) {
    # Skip iteration if files do not exist
    skip_to_next <- skip_check(i, params)
  
    if (skip_to_next) {
      return(NA)
    } else {
      # Get full data for that set of params
      full_i <- full[[i]]
      
      # Create a data frame of all combinations of n and strats
      combinations <- expand.grid(n = ns, strat = strats)
      
      # Iterate over each combination using pmap
      sub_results <- pmap(combinations, \(n, strat) return(run_mmrr_subsampled(i, params, n = n, strat = strat, full = full_i, site)))
      
      # Combine with full result
      results <- bind_rows(sub_results, full_i$full_result)
      
      return(results)
    }
  }, .options = furrr_options(seed = TRUE, packages = c("here", "vcfR", "adegenet", "stringr", "dplyr", "tidyr", "purrr")))
  
  return(bind_rows(results))
}

run_mmrr_full <- function(i, params, n = 2000) {
  # Skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  
  if (skip_to_next) {
    return(NA)
  } else {
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    
    # Subsample full data randomly
    s <- sample(nrow(gsd_df), n, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    # Run model on full data set
    result <- run_mmrr(gen_2k, gsd_df_2k)
    
    # save and format result
    full <- data.frame(params[i,], 
                       sampstrat = "full", 
                       nsamp = nrow(gsd_df_2k), 
                       result)
  }
  
  return(list(gen = gen, gsd_df = gsd_df, full_result = full))
}

run_mmrr_subsampled <- function(i, params, n, strat, full, site) {
  gen <- full$gen
  gsd_df <- full$gsd_df
  full_result <- full$full_result
  
  subIDs <- get_samples(params[i,], sampstrat = strat, nsamp = n, site = site)
  subgen <- gen[subIDs,]
  subgsd_df <- gsd_df[subIDs,]
  
  if (site){
    # Get sites
    siteIDs <- get_sites(params[i,], sampstrat = strat, nsamp = n)
    # Confirm that the number of sites matches the number of sample IDs
    stopifnot(length(subIDs) == length(siteIDs))
    # Calculate allele frequency by site (average)
    subgen <- data.frame(aggregate(subgen, list(siteIDs), FUN = mean)[-1])
    # Calculate env values by site
    subgsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN = mean)[-1])
  }
  
  # Run model on sub data set
  sub_stats <- run_mmrr(subgen, subgsd_df)
  
  # Calculate stats
  full_stats <- full_result %>% select(-K, -m, -phi, -H, -r, -sampstrat, -nsamp, -seed, -it)
  stats <- mmrr_stats(sub_stats, full_stats)

  # Save and format new result
  sub_result <- data.frame(params[i,], 
                           sampstrat = strat, 
                           nsamp = n, 
                           sub_stats, 
                           stats)
  
  return(sub_result)
}

# setup parallel session
future::plan(future::multisession, workers = 2)

# Run common MMRR operations
full <- future_map(1:nrow(params), \(i) run_mmrr_full(i, params, n = 2000),
                   .options = furrr_options(seed = TRUE, packages = c("here", "vcfR", "adegenet", "stringr", "dplyr", "tidyr", "purrr")))

# Perform MMRR analysis for site sampling
site_results <- perform_mmrr_analysis(params, ns = nsites, strats = sitestrats, full, site = TRUE)

# Perform MMRR analysis for individual sampling
ind_results <- perform_mmrr_analysis(params, ns = nsamp, strats = sampstrats, full = full, site = FALSE)

#stop cluster
stopCluster(cl)

write.csv(site_results, "outputs/2mmrr_sitesampling_results.csv", row.names = FALSE)
write.csv(ind_results, "outputs/2mmrr_indsampling_results.csv", row.names = FALSE)
