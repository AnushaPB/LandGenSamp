
#create filepath based on params index and data type (e.g. genetic data = gen, geospatial data = gsd, and adaptive loci = loci)
create_filepath <- function(i, type){
  #directory of data
  datadir <- "data/" 
  
  #set of parameter names in filepath form
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100)
  
  #different file patterns for different data types
  if(type == "gen"){filepath <- paste0(datadir, "GNX_mod-", paramset, "/it--",params[i,"it"],"/spp-spp_0/mod-", paramset,
                                       "_it--", params[i,"it"], "_t-1000_spp-spp_0.vcf")}
  if(type == "gsd"){filepath <- paste0(datadir, "GNX_mod-", paramset, "/it--",params[i,"it"],"/spp-spp_0/mod-", paramset,
                                       "_it--",params[i,"it"], "_t-1000_spp-spp_0.csv")}
  if(type == "loci"){filepath <- paste0(datadir, "nnloci_", paramset, ".csv")}
  
  return(filepath)
}


#Get gen data
get_gen <- function(filepath){
  #read vcf
  vcf <- read.vcfR(filepath)
  #convert to genlight from vcf
  genlight <- vcfR2genlight(vcf) #CHECK THIS
  #convert to matrix
  genmat <- as.matrix(genlight)
  #assign IDs from genlight to matrix rownames
  rownames(genmat) <- genlight@ind.names
  return(genmat)
}

#Get geospatial data
get_gsd <- function(filepath){
  gsd_df <- read.csv(filepath)
  gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
  gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
  rownames(gsd_df) <- gsd_df$idx
  return(gsd_df)
}

#general function to get data
get_data <- function(i, type){
  #different file patterns for different data types
  if(type == "gen"){
    filepath <- create_filepath(i, type)
    df <- get_gen(filepath)
  }
  
  if(type == "gsd"){
    filepath <- create_filepath(i, type)
    df <- get_gsd(filepath)
  }
  
  if(type == "loci"){
    filepath <- create_filepath(i, type)
    df <- read.csv(filepath)
  }
  
  return(df)
}

#get list of sampling IDs that correspond with parameter set, sampling strategy, and number of samples
get_samples <- function(param_set, sampstrat, nsamp){
  #param_set - vector of one set of parameters (e.g. params[i,])
  #sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  #nsamp - number of samples
  
  #Check if files for parameter exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  loci_filepath <- create_filepath(i, "loci")
  file_exists <- TRUE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){file_exists <- FALSE}
  if(!file_exists) { 
    print("File does not exist:")
    print(params[i,]) 
    } 
  stopifnot(file_exists)
  
  #directory of sample ID csvs (CHANGE)
  datadir <- "outputs/" 
  
  subIDs <- read.csv(paste0(datadir, "samples_", sampstrat, nsamp, ".csv"))
    
  subIDs <- subIDs[subIDs$K == param_set$K 
                   & subIDs$phi == param_set$phi
                   & subIDs$m == param_set$m 
                   & subIDs$seed == param_set$seed
                   & subIDs$H == param_set$H
                   & subIDs$r == param_set$r
                   & subIDs$it == param_set$it,]
  
  #confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)
  
  #remove parameter columnds and convert to vector of IDs
  subIDs <- subIDs[,!names(subIDs) %in% colnames(params)]
  subIDs <- unlist(subIDs)
  
  #confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))
  
  return(as.character(subIDs))
}

#function to calculate RMSE
#currently two functions because I am not sure whether to calculate the RMSE when the coeffs aren't signif
#on the one hand insignificant coeffs aren't meaningful, but it is hard to calculate averages for stats with NAs
#e.g. if the coeffs aren't signif they are likely to be far from the true values/should be represented in the mean stat
#currently using the first function for this reason (rmse_coeff)
rmse_coeff <- function(full_coeff, sub_coeff, full_p, sub_p, alpha = 0.05){
  sqerr <- (full_coeff - sub_coeff)^2
  res <- sqrt(mean(sqerr))
  return(res)
}

rmse_coeff_p <- function(full_coeff, sub_coeff, alpha = 0.05){
  sqerr <- (full_coeff - sub_coeff)^2
  res <- sqrt(mean(sqerr))
  res[full_p > alpha] <- NA
  res[sub_p > alpha] <- NA
  return(res)
}

#GENERAL OBJECTS (objects used in multiple scripts)
#nloci 
nloci = 10000

#number of points to sample
npts <- c(36, 81, 144, 225, 324)
#sampling strategies
sampstrats <- c("rand", "grid", "trans", "envgeo")
#landscape dimensions (square)
dim = 40

#Create dataframe with all variable combos
params <- expand.grid(K = c(2, 4), 
                      phi = c(0.1, 0.5),
                      m = c(0.25, 1.0),
                      seed = c(1, 2, 3),
                      H = c(0.05 , 0.5),
                      r = c(0.3, 0.6),
                      it = 1:10)

#TESTING PARAMS (REMOVE LATER)
#params <- expand.grid(K = c(2), 
                     # phi = c(0.5),
                     # m = c(0.25),
                     # seed = c(1),
                     # H = c(0.50),
                     # r = c(0.30),
                     # it = 1)


