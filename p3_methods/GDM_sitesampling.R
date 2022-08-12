set.seed(42)

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")
source("GDM_functions.R")

#register cores
cores <- 25
cl <- makeCluster(cores) 
registerDoParallel(cl)

res_gdm <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("vcfR", "gdm", "adegenet", "stringr", "dplyr", "here")) %dopar% {
 
  #skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if(skip_to_next) { result <- NA } 
  
  #run GDM
  tryCatch({
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run analysis using subsample
    
    #run model on full data set
    full_result <- run_gdm(gen_2k, gsd_df_2k, distmeasure = "euc")
    fullratio <- (full_result$env1_coeff + full_result$env2_coeff)/full_result$geo_coeff
    result <- data.frame(params[i,], 
                         sampstrat = "full", 
                         nsamp = 2000, 
                         full_result, 
                         ratio = fullratio,
                         env1_err = NA, 
                         env2_err = NA, 
                         geo_err = NA,
                         ratio_err = NA)
    result <- sapply(result, as.character)
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsite)
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
        sub_result <- run_gdm(sitegen, sitegsd_df, distmeasure = "euc")
        
        #calculate err if not null
        if(sub_result$env1_coeff == "NULL" | sub_result$env2_coeff == "NULL" | sub_result$geo_coeff == "NULL"){
          subratio <- "NULL"
          env1_err <- "NULL"
          env2_err <- "NULL"
          geo_err <- "NULL"
          ratio_err <- "NULL"
        } else {
          subratio <- (sub_result$env1_coeff + sub_result$env2_coeff)/sub_result$geo_coeff
          env1_err <- err_coeff(full_result$env1_coeff, sub_result$env1_coeff)
          env2_err <- err_coeff(full_result$env2_coeff, sub_result$env2_coeff)
          geo_err <- err_coeff(full_result$geo_coeff, sub_result$geo_coeff)
          ratio_err <- err_coeff(fullratio, subratio)
        }
        
        
        #save and format new result
        sub_result <- data.frame(params[i,], 
                                 sampstrat = sampstrat, 
                                 nsamp = nsamp, 
                                 sub_result, 
                                 ratio = subratio,
                                 env1_err = env1_err, 
                                 env2_err = env2_err, 
                                 geo_err = geo_err, 
                                 ratio_err = ratio_err)
        
        #bind results 
        sub_result <- sapply(sub_result, as.character)
        result <- bind_rows(result, sub_result)
      }
    }
  }
    # end trycatch
  }, error = function(e) {
    err <<- conditionMessage(e)
    write.table(err, "error_msg.txt")
    write.csv(gen_2k, "error_gen_2k.csv", row.names = FALSE)
    write.csv(gsd_df_2k, "error_gsd_df_2k.csv", row.names = FALSE)
    write.csv(sitegen, "error_sitegen.csv", row.names = FALSE)
    write.csv(sitegsd_df, "error_sitegsd_df.csv", row.names = FALSE)
    write.csv(gen, "error_gen.csv", row.names = FALSE)
    write.csv(gsd_df, "error_gsd_df.csv", row.names = FALSE)
    
    message(err)
    
    stop(err)})
  
  return(result)
  
  gc()
  
}
  
  

#stop cluster
stopCluster(cl)

write.csv(res_gdm, "outputs/gdm_sitesampling_results.csv", row.names = FALSE)

