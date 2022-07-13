


for(i in 1:nrow(params)){
  #vcfR
  library("vcfR")
  library("gdm")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100)
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  
  skip_to_next <- FALSE
  if(file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run gdm
  if(skip_to_next == FALSE){
    
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_gdm(gen_2k, gsd_df_2k)
    result <- data.frame(sampstrat = "full", nsamp = nrow(gsd_df), full_result, env1_rmse = NA, env2_rmse = NA, geo_rmse = NA)
    
    #write full datafile (temp)
    csv_file <- paste0("gdm_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_gdm(subgen, subgsd_df)
        
        #calculate RMSE
        env1_rmse <- rmse_coeff(full_result$env1_coeff, sub_result$env1_coeff)
        env2_rmse <- rmse_coeff(full_result$env2_coeff, sub_result$env2_coeff)
        geo_rmse <- rmse_coeff(full_result$geo_coeff, sub_result$geo_coeff)
        
        #save and format new result
        sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, sub_result, 
                                 env1_rmse = env1_rmse, env2_rmse = env2_rmse, geo_rmse = geo_rmse)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  #return(result)
  
  gc()
  
}
