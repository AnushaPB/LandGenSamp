set.seed(42)

library("here") #paths
library("foreach")
library("doParallel")
library("dplyr")
library("tidyr")
#read in general functions and objects
source(here("general_functions.R"))

##############
#  MISMATCH  #
##############

#register cores
cl <- makeCluster(20)
registerDoParallel(cl)

system.time(
res_mismatch <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("dplyr")
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
    z1mis <- abs(gsd_df$z1 - gsd_df$env1)
    z2mis <- abs(gsd_df$z2 - gsd_df$env2)
    mis <- z1mis + z2mis

    #calculate phenotype-env correlation
    mod1 <- summary(lm(gsd_df$z1 ~ gsd_df$env1))
    mod2 <- summary(lm(gsd_df$z2 ~ gsd_df$env2))
    p1 <- mod1$coefficients[-1,"Pr(>|t|)"]
    p2 <- mod2$coefficients[-1,"Pr(>|t|)"]
    coeff1 <- mod1$coefficients[-1,"Estimate"]
    coeff2 <- mod2$coefficients[-1,"Estimate"]

    result <- 
      data.frame(params[i,], 
      sampstrat = "full", 
      nsamp = nrow(gsd_df), 
      mismatch_mean = mean(mis), 
      mismatch_max = max(mis),
      mod1_p = p1,
      mod2_p = p2,
      mod1_coeff = coeff1,
      mod2_coeff = coeff2)

    for(nsamp in nsamps){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgsd_df <- gsd_df[subIDs,]
        
        #calculate phenotypic mistmatch
        z1mis <- abs(subgsd_df$z1 - subgsd_df$env1)
        z2mis <- abs(subgsd_df$z2 - subgsd_df$env2)
        mis <- z1mis + z2mis
        
        sub_result <- data.frame(params[i,],  sampstrat = sampstrat, nsamp = nsamp, mismatch_mean = mean(mis), mismatch_max = max(mis))
        
        #bind results
        result <- bind_rows(result, sub_result)
      }
    }
    
    result <-
      result %>%
      mutate(mod1_sig = mod1_p < 0.05, mod2_sig = mod2_p < 0.05, mod_sig = (mod1_sig + mod2_sig)/2)
  }
  
  return(result)
  
  gc()
}
)

#stop cluster
stopCluster(cl)

write.csv(res_mismatch, here("p3_methods", "outputs" , "mismatch_results.csv"), row.names = FALSE)

