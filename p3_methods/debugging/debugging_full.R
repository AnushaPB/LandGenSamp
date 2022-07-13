
#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

system.time(
res_test <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  
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
    loci_filepath <- create_filepath(i, "loci")
    skip_to_next <- FALSE
    if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
    if(skip_to_next) { print("File does not exist:")
      print(params[i,]) } 
    if(skip_to_next) { result <- NA } 
    
    #run LFMM
    if(skip_to_next == FALSE){
      #run model on full data set
      full_result <- paramset
      result <- data.frame(sampstrat = "full", nsamp = 2000, test = full_result)
      
      #write full datafile (temp)
      #csv_file <- paste0("outputs/LFMM/LFMM_results_",paramset,".csv")
      #write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
      
      for(nsamp in npts){
        for(sampstrat in sampstrats){
          #subsample from data based on sampling strategy and number of samples

          #run analysis using subsample
          sub_result <- paramset
          
          #save and format new result
          sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, test = sub_result)
          
          ##export data to csv (temp)
          #csv_df <- read.csv(csv_file)
          #csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
          #write.csv(csv_df, csv_file, row.names = FALSE)
          
          #bind results
          result <- rbind.data.frame(result, sub_result)
        }
      }
    }
    
    
    return(result)
    
    gc()
  }
)

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_test)
dim(res_test)
dim(params)

res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {

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
  loci_filepath <- create_filepath(i, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    loci_df <- get_data(i, "loci")
    
    #run model on full data set
    result <- data.frame(sampstrat = "full", nsamp = 2000, full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LFMM/LFMM_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_lfmm(subgen, subgsd_df, loci_df, K = full_result$K)
        
        #save and format new result
        sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, sub_result)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  #end pdf()
  dev.off()
  
  return(result)
  
  gc()
}
)

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)