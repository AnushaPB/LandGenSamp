

#get list of sampling IDs that correspond with parameter set, sampling strategy, and number of samples
get_samples <- function(param_set, params, sampstrat, nsamp, datadir =  here(dirname(getwd()), "p2_sampling", "outputs")){
  #param_set - vector of one set of parameters (e.g. params[i,])
  #params - full set of parameters
  #sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  #nsamp - number of samples
  
  subIDs <- read.csv(paste0(datadir, "/samples_", sampstrat, nsamp, ".csv"))
    
  subIDs <- subIDs[subIDs$K == param_set$K 
                   & subIDs$phi == param_set$phi
                   & subIDs$m == param_set$m 
                   & subIDs$seed == param_set$seed
                   & subIDs$H == param_set$H
                   & subIDs$r == param_set$r
                   & subIDs$it == param_set$it,]
  
  #confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)
  
  #remove parameter columns and convert to vector of IDs
  subIDs <- subIDs[,!names(subIDs) %in% colnames(params)]
  subIDs <- unlist(subIDs)
  
  #confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))
  
  return(as.character(subIDs))
}

#get list of sampling IDs that correspond with parameter set, sampling strategy, and number of samples
get_samples <- function(param_set, params, sampstrat, nsamp,  datadir =  here(dirname(getwd()), "p2_sampling", "outputs")){
  #param_set - vector of one set of parameters (e.g. params[i,])
  #params - full set of parameters
  #sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  #nsamp - number of samples
  
  subIDs <- read.csv(paste0(datadir, "/site_samples_", sampstrat, nsamp, ".csv"))
  
  subIDs <- subIDs[subIDs$K == param_set$K 
                   & subIDs$phi == param_set$phi
                   & subIDs$m == param_set$m 
                   & subIDs$seed == param_set$seed
                   & subIDs$H == param_set$H
                   & subIDs$r == param_set$r
                   & subIDs$it == param_set$it,]
  
  #confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)
  
  #remove parameter columns and convert to vector of IDs
  subIDs <- subIDs[,!names(subIDs) %in% colnames(params)]
  subIDs <- unlist(subIDs)
  
  #confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))
  
  return(as.character(subIDs))
}

#get list of site IDs that correspond with parameter set, sampling strategy, and number of samples (and sample IDs)
get_sites <- function(param_set, params, sampstrat, nsamp,  datadir =  here(dirname(getwd()), "p2_sampling", "outputs")){
  #param_set - vector of one set of parameters (e.g. params[i,])
  #params - full set of parameters
  #sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  #nsamp - number of samples
  
  subIDs <- read.csv(paste0(datadir, "/site_ids_", sampstrat, nsamp, ".csv"))
  
  subIDs <- subIDs[subIDs$K == param_set$K 
                   & subIDs$phi == param_set$phi
                   & subIDs$m == param_set$m 
                   & subIDs$seed == param_set$seed
                   & subIDs$H == param_set$H
                   & subIDs$r == param_set$r
                   & subIDs$it == param_set$it,]
  
  #confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)
  
  #remove parameter columns and convert to vector of IDs
  subIDs <- subIDs[,!names(subIDs) %in% colnames(params)]
  subIDs <- unlist(subIDs)
  
  #confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))
  
  return(subIDs)
}


######################################################
# GENERAL OBJECTS (objects used in multiple scripts) #
######################################################

#number of sites
npts <- c(9, 16, 25)
#sampling strategies
sampstrats <- c("rand", "equi", "envgeo")

