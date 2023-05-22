library("here") #paths
library("gdm") #GDM
library("vcfR")
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")

set.seed(42)


for(i in 1:nrow(params)) {

  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  skip_to_next <- FALSE
  if(file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run GDM
  if(skip_to_next == FALSE){
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsite)
      }
    }
  }
  
}


