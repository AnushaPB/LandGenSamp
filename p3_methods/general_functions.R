
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
                     "_r",params[i,"r"]*100,
                     "_it--",params[i,"it"])
  
  #different file patterns for different data types
  if(type == "gen"){filepath <- paste0(datadir, "mod-", paramset, "_t-500_spp-spp_0.vcf")}
  if(type == "gsd"){filepath <- paste0(datadir, "mod-", paramset, "_t-500_spp-spp_0.csv")}
  if(type == "loci"){filepath <- paste0(datadir, "nnloci_", paramset, ".csv")}
  
  return(filepath)
}


#Get gen data
get_gen <- function(filepath){
  vcf <- read.vcfR(filepath)
  x <- vcfR2genlight(vcf) #CHECK THIS
  gen <- as.matrix(x)
  return(gen)
}

#Get geospatial data
get_gsd <- function(filepath){
  gsd_df <- read.csv(filepath)
  gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
  gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
  return(gsd_df)
}

#Create dataframe with all variable combos
params <- expand.grid(K = c(2,5), #CHANGE TO 1.5 and 3
                      phi = c(0.1,0.5),
                      m = c(0.25,1.0),
                      seed = c(1,2,3),
                      H = c(0.05,0.5),
                      r = c(0.3, 0.6),
                      it = 1:10)

#TESTING PARAMS (REMOVE LATER)
params <- expand.grid(K = c(2), #CHANGE TO 1.5 and 3
                      phi = c(0.5),
                      m = c(0.5),
                      seed = c(1),
                      H = c(0.5),
                      r = c(0.6),
                      it = 1)

#define nloci 
nloci = 10000
