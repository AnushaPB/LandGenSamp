set.seed(42)

library("here") #paths
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")

##############
#  MISMATCH  #
##############


#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-3) #not to overload your computer
registerDoParallel(cl)

system.time(
res_mismatch <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  
  #skip iteration if files do not exist
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  skip_to_next <- FALSE
  if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
                      print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  if(skip_to_next == FALSE){
    gsd_df <- get_data(i, params = params, "gsd")
    
    #calculate phenotypic mistmatch
    z1mis <- gsd_df$z1 - gsd_df$env1
    z2mis <- gsd_df$z2 - gsd_df$env2
    mis <- c(z1mis, z2mis)
    
    result <- data.frame(params[i,], sampstrat = "full", nsamp = nrow(gsd_df), mismatch = mean(mis))
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgsd_df <- gsd_df[subIDs,]
        
        #calculate phenotypic mistmatch
        z1mis <- subgsd_df$z1 - subgsd_df$env1
        z2mis <- subgsd_df$z2 - subgsd_df$env2
        mis <- c(z1mis, z2mis)
        
        sub_result <- data.frame(params[i,],  sampstrat = sampstrat, nsamp = nsamp, mismatch = mean(mis))
        
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

write.csv(res_mismatch, "outputs/mismatch_results.csv", row.names = FALSE)

