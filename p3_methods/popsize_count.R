set.seed(42)

library("here") 
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")

#register cores
cores <- 4
cl <- makeCluster(cores)
#not to overload your computer
registerDoParallel(cl)

system.time(
res_popsize <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
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
  
  #get pop size
  if(skip_to_next == FALSE){
    gsd_df <- get_data(i, params = params, "gsd")
    result <- data.frame(params[i,], popsize = nrow(gsd_df))
  }
  
  #end pdf()
  # dev.off()
  
  return(result)
  
  gc()
}
)

#stop cluster
stopCluster(cl)

write.csv(res_popsize, "outputs/popsize_results.csv", row.names = FALSE)
