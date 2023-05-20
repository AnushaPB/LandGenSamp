set.seed(42)

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")
source("mmrr_functions.R")

#register cores
cores <- 20
cl <- makeCluster(cores)
registerDoParallel(cl)

res_mmrr <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("here", "vcfR", "adegenet", "stringr")) %dopar% {
  
  #skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if(skip_to_next) { result <- NA } 
  
  #run MMRR
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_mmrr(gen_2k, gsd_df_2k)
    fullratio <- (abs(full_result$env1_coeff) + abs(full_result$env2_coeff))/abs(full_result$geo_coeff)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = nrow(gsd_df_2k), full_result, 
                         ratio = fullratio, 
                         env1_err = NA, env2_err = NA, geo_err = NA, ratio_err = NA)
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsite)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #get sites
        siteIDs <- get_sites(params[i,], params, sampstrat, nsite)
        #confirm that number of sites matches number of sample IDs
        stopifnot(length(subIDs) == length(siteIDs))
        #calculate allele frequency by site (average)
        sitegen <- data.frame(aggregate(subgen, list(siteIDs), FUN=mean)[,-1])
        #calculate env values by site
        sitegsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN=mean)[,-1]) 
        #confirm that number of values matches number of sites
        stopifnot(nsite == nrow(sitegsd_df) & nsite == nrow(sitegen))
        
        #run analysis using subsample
        sub_result <- run_mmrr(sitegen, sitegsd_df)
        
        #calculate err
        env1_err <- err_coeff(full_result$env1_coeff, sub_result$env1_coeff)
        env2_err <- err_coeff(full_result$env2_coeff, sub_result$env2_coeff)
        geo_err <- err_coeff(full_result$geo_coeff, sub_result$geo_coeff)
        
        subratio <- (abs(sub_result$env1_coeff) + abs(sub_result$env2_coeff))/abs(sub_result$geo_coeff)
        ratio_err <- err_coeff(fullratio, subratio)
        
        #save and format new result
        sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsite, sub_result, 
                                 ratio = subratio,
                                 env1_err = env1_err, env2_err = env2_err, geo_err = geo_err,
                                 ratio_err = ratio_err)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
        
      }
    }
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

write.csv(res_mmrr, "outputs/mmrr_sitesampling_results.csv", row.names = FALSE)
