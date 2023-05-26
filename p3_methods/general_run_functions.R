
run_method <- function(method, sampling = c("individual", "site"), ncores = NULL){
  # Parallel processing libraries
  library(furrr)
  library(dplyr)
  library(tidyr)
  library(here)
  library(vcfR)
  library(lfmm)
  library(gdm)
  library(vegan)
  
  # Read in general functions and objects
  source(here("general_functions.R"))
  source(here("p3_methods", "general_run_functions.R"))
  source(here("p3_methods", "mmrr_functions.R"))
  source(here("p3_methods", "lfmm_functions.R"))
  source(here("p3_methods", "rda_functions.R"))
  source(here("p3_methods", "gdm_functions.R"))
  
  # setup parallel session
  if (is.null(ncores)) ncores <- 20
  
  # Run common operations
  if (method == "mmrr" | method == "gdm") full_result <- run_full(params, method = method, ncores = ncores) else full_result <- NULL
  
  # run analysis for individual sampling
  if (any(sampling == "individual")){
    ind_results <- run_analysis(params, ns = nsamps, strats = sampstrats, method = method, full_result = full_result, site = FALSE, ncores = ncores)
    path <- here("p3_methods", "outputs", paste0(method, "_indsampling_results.csv"))
    write.csv(ind_results, path, row.names = FALSE)
  }
  
  # run analysis for site sampling
  if (any(sampling == "site")){
    site_results <- run_analysis(params, ns = nsites, strats = sitestrats, method = method, full_result = full_result, site = TRUE, ncores = ncores)
    path <- here("p3_methods", "outputs", paste0(method, "_sitesampling_results.csv"))
    write.csv(site_results, path, row.names = FALSE)
  }
}

run_analysis <- function(params, ns, strats, method, full_result = NULL, site = FALSE, ncores = 10) {
  
  future::plan(future::multisession, workers = ncores)
  
  results <- future_map(1:nrow(params), \(i) {
    # Skip iteration if files do not exist
    skip_to_next <- skip_check(i, params)
    if (skip_to_next) return(NA)

    # Get full result
    if (!is.null(full_result)) full_result_i <- full_result[[i]] else full_result_i <- NULL
    
    # Get full data
    gen <- get_data(i, params = params, "dos")
    gsd_df <- get_data(i, params = params, "gsd")
    
    # Create a data frame of all combinations of n and strats
    combinations <- expand.grid(n = ns, strat = strats)
    
    # Iterate over each combination using pmap
    results <-
      pmap(combinations, \(n, strat) return(
        run_subsampled(
          i,
          params = params,
          n = n,
          strat = strat,
          gen = gen, 
          gsd_df = gsd_df,
          full_result = full_result_i,
          method = method,
          site = site
        )
      ))
    
    # Combine with full result if mmrr/gdm
    if (!is.null(full_result_i)) results <- bind_rows(results, full_result_i)
    
    return(results)

  }, .options = furrr_options(seed = TRUE, packages = get_packages()), .progress = TRUE)
  
  results <- bind_rows(results)
  print(results)
  
  ## Shut down parallel workers
  future::plan("sequential")
  
  return(results)
}

run_full <- function(params, method, n = 2000, ncores = 10){
  
  future::plan(future::multisession, workers = ncores)
  
  future_map(
    1:nrow(params),
    \(i) run_full_helper(
      i,
      params = params,
      method = method,
      n = n
    ),
    .options = furrr_options(seed = TRUE, packages = get_packages())
  )
  
  ## Shut down parallel workers
  future::plan("sequential")
}

run_full_helper <- function(i, params, method, n = 2000) {
  # Skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if (skip_to_next) return(NA)

  gen <- get_data(i, params = params, "dos")
  gsd_df <- get_data(i, params = params, "gsd")
  
  # Subsample full data randomly
  s <- sample(nrow(gsd_df), n, replace = FALSE)
  gen_2k <- gen[s,]
  gsd_df_2k <- gsd_df[s,]
  
  # Run model on full data set
  run_method <- get_method(method, type = "run")
  result <- run_method(gen_2k, gsd_df_2k)
  
  # save and format result
  full_result <- data.frame(params[i, ],
                            sampstrat = "full",
                            nsamp = nrow(gsd_df_2k),
                            result)
  
  return(full_result)
}


run_subsampled <- function(i, params, n, strat, gen, gsd_df, full_result, method, site) {
  
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
  run_method <- get_method(method, type = "run")
  
  if (method == "mmrr" | method == "gdm") {
    sub_stats <- run_method(subgen, subgsd_df)
    # Calculate stats
    full_stats <- full_result %>% select(-K, -m, -phi, -H, -r, -sampstrat, -nsamp, -seed, -it)
    method_stat <- get_method(method, type = "stat")
    stats <- method_stat(sub_stats, full_stats)
    stats <- bind_cols(sub_stats, stats)
  } else {
    stats <- run_method(subgen, subgsd_df)
  }
    
  # Save and format new result
  sub_result <- data.frame(params[i,], 
                           sampstrat = strat, 
                           nsamp = n, 
                           stats)
  
  return(sub_result)
}

get_method <- function(method, type = "run"){
  if (type == "run") {
    if (method == "mmrr") return(run_mmrr)
    if (method == "gdm") return(run_gdm)
    if (method == "lfmm") return(run_lfmm)
    if (method == "rda") return(run_rda)
  }
  
  if (type == "stat") {
    if (method == "mmrr") return(stat_mmrr)
    if (method == "gdm") return(stat_gdm)
    if (method == "lfmm") return(stat_lfmm)
    if (method == "rda") return(stat_rda)
  }
  
  stop("invalid input")
}

get_packages <- function(){
  c("here", "vcfR", "adegenet", "stringr", "dplyr", "tidyr", "purrr", "lfmm", "AssocTests", "gdm", "vegan")
}