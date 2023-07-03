
run_method <- function(method, sampling = c("individual", "site"), ncores = NULL){
  # Read in general functions and objects
  source(here::here("general_functions.R"))
  source(here::here("p3_methods", "general_run_functions.R"))
  source(here::here("p3_methods", "GEA_functions.R"))
  source(here::here("p3_methods", "IBDIBE_functions.R"))
  
  # set cores
  if (is.null(ncores)) ncores <- 25
  
  # make cluster
  cl <- parallel::makeCluster(ncores) 
  # using doSNOW to get progress bar
  doSNOW::registerDoSNOW(cl)
  
  # Run common operations
  if (method %in% c("mmrr", "mmrr2", "gdm", "gdm2", "lfmm_fullK"))
    full_result <- run_full(params, method = method, ncores = ncores, n = 1000) 
  else
    full_result <- NULL
  
  # run analysis for individual sampling
  if (any(sampling == "individual")){
    ind_results <- run_analysis(params, ns = nsamps, strats = sampstrats, method = method, full_result = full_result, site = FALSE, ncores = ncores)
    path <- here::here("p3_methods", "outputs", paste0(method, "_indsampling_results.csv"))
    write.csv(ind_results, path, row.names = FALSE)
  }
  
  # run analysis for site sampling
  if (any(sampling == "site")){
    site_results <- run_analysis(params, ns = nsites, strats = sitestrats, method = method, full_result = full_result, site = TRUE, ncores = ncores)
    path <- here::here("p3_methods", "outputs", paste0(method, "_sitesampling_results.csv"))
    write.csv(site_results, path, row.names = FALSE)
  }
  
  ## Shut down parallel workers
  snow::stopCluster(cl)
  
}


run_analysis2 <- function(params, ns, strats, method, full_result = NULL, site = FALSE, ncores = 25) {
  
  future::plan(future::multisession, workers = ncores)
  
  results <-
    furrr::future_map(
      1:nrow(params),
      ~ run_analysis_helper(
        i = .x,
        params = params,
        ns = ns,
        strats = strats,
        method = method,
        full_result = full_result,
        site = site
      ),
      .options = furrr::furrr_options(seed = TRUE, packages = get_packages()),
      .progress = TRUE
    )
  
  results <- dplyr::bind_rows(results)
  
  future::plan("sequential")
  
  return(results)
}


run_analysis <- function(params, ns, strats, method, full_result = NULL, site = FALSE) {
  pb <- txtProgressBar(max = nrow(params), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  results <-
    foreach::foreach(
      i = 1:nrow(params),
      .combine = dplyr::bind_rows,
      .packages = get_packages(),
      .options.snow = list(progress = progress)
    ) %dopar% {
      # Read in general functions and objects
      source(here::here("general_functions.R"))
      source(here::here("p3_methods", "general_run_functions.R"))
      source(here::here("p3_methods", "GEA_functions.R"))
      source(here::here("p3_methods", "IBDIBE_functions.R"))
      
      return(run_analysis_helper(
        i = i,
        params = params,
        ns = ns,
        strats = strats,
        method = method,
        full_result = full_result,
        site = site
      ))
    }
  
  # end progress
  close(pb)
  
  return(results)
}

run_analysis_helper <- function(i, params, ns, strats, method, full_result = NULL, site = FALSE){
  # Skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if (skip_to_next) return(NA)
  
  # Get full result
  if (!is.null(full_result)) full_result_i <- full_result[i,] else full_result_i <- NULL
  
  # Get full data
  gen <- get_data(i, params = params, "dos")
  gsd_df <- get_data(i, params = params, "gsd")
  
  # Create a data frame of all combinations of n and strats
  combinations <- expand.grid(n = ns, strat = strats)
  
  # Iterate over each combination using pmap
  results <-
    purrr::pmap(combinations, \(n, strat) return(
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
  
  # bind rows
  results <- results %>% dplyr::bind_rows()
  
  # Combine with full result if mmrr/gdm
  if (!is.null(full_result_i)) results <- dplyr::bind_rows(results, full_result_i)
  
  # remove large objects
  rm("gen")
  rm("gsd_df")
  gc()
  
  return(results)
}


run_full2 <- function(params, method, n = 1000, ncores = 10){
  
  future::plan(future::multisession, workers = ncores)
  
  results <- future_map(
    1:nrow(params),
    \(i) run_full_helper(
      i,
      params = params,
      method = method,
      n = n
    ),
    .options = furrr::furrr_options(seed = TRUE, packages = get_packages()),
    .progress = TRUE
  ) %>%
    dplyr::bind_rows()
  
  ## Shut down parallel workers
  future::plan("dsequential")
  
  return(results)
}


run_full <- function(params, method, n = 2000){
  pb <- txtProgressBar(max = nrow(params), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  
  results <- foreach(i = 1:nrow(params),
                     .combine = dplyr::bind_rows,
                     .packages = get_packages(),
                     .options.snow = list(progress = progress)) %dopar% {
                       # Read in general functions and objects
                       source(here::here("general_functions.R"))
                       source(here::here("p3_methods", "general_run_functions.R"))
                       source(here::here("p3_methods", "GEA_functions.R"))
                       source(here::here("p3_methods", "IBDIBE_functions.R"))
                       
                       return(run_full_helper(i,
                                       params = params,
                                       method = method,
                                       n = n))
                     }

  # end progress
  close(pb)
  
  return(results)
}

run_full_helper <- function(i, params, method, n = 2000) {
  # Skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if (skip_to_next) return(NA)

  gen <- get_data(i, params = params, "dos")
  gsd_df <- get_data(i, params = params, "gsd")
  
  # Subsample full data randomly
  set.seed(124)
  s <- sample(nrow(gsd_df), n, replace = FALSE)
  gen <- gen[s,]
  gsd_df <- gsd_df[s,]
  
  # Run model on full data set
  if (method == "lfmm_fullK") {
    K <- get_K_tess(gen, coords = gsd[,c("x", "y")], Kvals = 1:9)
    result <- data.frame(K_factor = K)
  } else {
    run_method <- get_method(method, type = "run")
    result <- run_method(gen, gsd_df)
  }
  
  # save and format result
  full_result <- data.frame(params[i, ],
                            sampstrat = "full",
                            nsamp = nrow(gsd_df),
                            result)
  
  
  # remove large objects
  rm("gen")
  rm("gsd_df")
  gc()
  
  return(full_result)
}


run_subsampled <- function(i, params, n, strat, gen, gsd_df, full_result, method, site, full_K = FALSE) {
  
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
  
  if (method %in% c("mmrr", "mmrr2", "gdm", "gdm2")) {
    sub_stats <- run_method(subgen, subgsd_df)
    # Calculate stats
    full_stats <- full_result %>% dplyr::select(-K, -m, -phi, -H, -r, -sampstrat, -nsamp, -seed, -it)
    method_stat <- get_method(method, type = "stat")
    stats <- method_stat(sub_stats, full_stats)
    stats <- dplyr::bind_cols(sub_stats, stats)
  } 
  
  if (method == "lfmm_fullK") stats <- run_lfmm(subgen, subgsd_df, K = full_result$K_factor)
  
  if (method == "lfmm" | method == "rda")  stats <- run_method(subgen, subgsd_df)
    
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
    if (method == "mmrr2") return(run_mmrr2)
    if (method == "gdm") return(run_gdm)
    if (method == "gdm2") return(run_gdm2)
    if (method == "lfmm") return(run_lfmm)
    if (method == "rda") return(run_rda)
  }
  
  if (type == "stat") {
    if (method %in% c("mmrr", "mmrr2", "gdm", "gdm2")) return(stat_ibdibe)
  }
  
  stop("invalid input")
}


