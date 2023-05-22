source("general_functions.R")
library("here")
library("foreach")
library("doParallel")

set.seed(42)


missing <- data.frame()
for(i in 1:nrow(params)) {
    library("here")
    
    #create file path
    gsd_filepath <- create_filepath(i, params = params, "gsd")
    
    #skip iteration if file does not exist
    skip_to_next <- FALSE
    if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
    if(skip_to_next) { print("File does not exist:")
      print(params[i,]) 
      missing <- rbind(missing, params[i,])} 
    if(skip_to_next) { result <- NA } 

}


missing <- apply(missing, 2, as.factor)


summary(unique(missing))