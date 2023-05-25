
make_dosage <- function(params){
  future::plan(future::multisession, workers = 20)
  future_map(
    1:nrow(params),
    \(i) {
      gen <- get_data(i, params = params, "gen")
      file_path <- create_filepath(i, params, type = "gen")
      new_file_path <- gsub("mod-(.*?)_", "dos-\\1_", file_path)
      new_file_path <- gsub(".vcf", ".csv", new_file_path)
      write.csv(gen, new_file_path, row.names = FALSE)
      }, .options = furrr_options(seed = TRUE, packages = get_packages())
  )
}

run_method <- function(method, sampling = c("individual", "site"), ncores = NULL){
  # Parallel processing libraries
  library(furrr)
  library(dplyr)
  
  # Read in general functions and objects
  source(here("general_functions.R"))
  source(here("p3_methods", "mmrr_functions.R"))
  source(here("p3_methods", "lfmm_functions.R"))
  source(here("p3_methods", "rda_functions.R"))
  source(here("p3_methods", "gdm_functions.R"))
  
  # setup parallel session
  if (is.null(ncores)) ncores <- 20
  future::plan(future::multisession, workers = ncores)
  
  # Run common operations
  full <- run_full(params, method = method)
  
  # run analysis for individual sampling
  if (any(sampling == "individual")){
    ind_results <- run_analysis(params, ns = nsamps, strats = sampstrats, full = full, method = method, site = FALSE)
    path <- here("p3_methods", "outputs", paste0(method, "_indsampling_results.csv"))
    write.csv(ind_results, path, row.names = FALSE)
  }
  
  # run analysis for site sampling
  if (any(sampling == "site")){
    site_results <- run_analysis(params, ns = nsites, strats = sitestrats, full = full, method = method, site = TRUE)
    path <- here("p3_methods", "outputs", paste0(method, "_sitesampling_results.csv"))
    write.csv(site_results, path, row.names = FALSE)
  }
}

run_analysis <- function(params, ns, strats, full, method, site = FALSE) {
  results <- future_map(1:nrow(params), \(i) {
    # Skip iteration if files do not exist
    skip_to_next <- skip_check(i, params)
    if (skip_to_next) return(NA)

    # Get full data for that set of params
    full_i <- full[[i]]
    
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
          full = full_i,
          method = method,
          site = site
        )
      ))
    
    # Combine with full result if mmrr/gdm
    if (method == "mmrr" | method == "gdm") results <- bind_rows(results, full_i$full_result)
    
    return(results)

  }, .options = furrr_options(seed = TRUE, packages = get_packages()))
  
  return(bind_rows(results))
}

run_full <- function(params, method, n = 2000){
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
}

run_full_helper <- function(i, params, method, n = 2000) {
  # Skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if (skip_to_next) return(NA)

  gen <- get_data(i, params = params, "dos")
  gsd_df <- get_data(i, params = params, "gsd")
  
  if (method == "mmrr" | method == "gdm"){
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
  } else full_result <- NULL
  
  return(list(gen = gen, gsd_df = gsd_df, full_result = full_result))
}


run_subsampled <- function(i, params, n, strat, full, method, site) {
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