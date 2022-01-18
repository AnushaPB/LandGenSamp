#library to create paths
library("here")


#####################
# GENERAL FUNCTIONS #
#####################

#create filepath based on params index and data type (e.g. genetic data = gen, geospatial data = gsd, and adaptive loci = loci)
#FOR FILES NOT NESTED IN SUBFOLDERS
create_filepath <- function(i, params, type, datadir = here(dirname(getwd()), "p1_gnxsims", "parallel", "LGS_data")){
  
  #set of parameter names in filepath form
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100)
  
  #different file patterns for different data types
  if(type == "gen"){filepath <- paste0(datadir, "/mod-", paramset,
                                       "_it-", params[i,"it"], "_t-1000_spp-spp_0.vcf")}
  if(type == "gsd"){filepath <- paste0(datadir, "/mod-", paramset,
                                       "_it-",params[i,"it"], "_t-1000_spp-spp_0.csv")}
  if(type == "loci"){filepath <- paste0(datadir, "/nnloci_", paramset, ".csv")}
  
  print(filepath)
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
get_data <- function(i, params, type){
  #different file patterns for different data types
  if(type == "gen"){
    filepath <- create_filepath(i, params, type)
    print(filepath)
    df <- get_gen(filepath)
  }
  
  if(type == "gsd"){
    filepath <- create_filepath(i, params, type)
    print(filepath)
    df <- get_gsd(filepath)
  }
  
  if(type == "loci"){
    filepath <- create_filepath(i, params, type)
    print(filepath)
    df <- read.csv(filepath)
  }
  
  return(df)
}

#get list of sampling IDs that correspond with parameter set, sampling strategy, and number of samples
get_samples <- function(param_set, params = params, sampstrat, nsamp, outdir = here(dirname(getwd()), "p2_sampling", "outputs")){
  #param_set - vector of one set of parameters (e.g. params[i,])
  #sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  #nsamp - number of samples
  
  #Check if files for parameter exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  print(gen_filepath)
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  print(gsd_filepath)
  loci_filepath <- create_filepath(i, params = params, "loci")
  print(loci_filepath)
  file_exists <- TRUE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){file_exists <- FALSE}
  if(!file_exists) { 
    print("File does not exist:")
    print(params[i,]) 
  } 
  stopifnot(file_exists)
  
  subIDs <- read.csv(paste0(outdir, "/samples_", sampstrat, nsamp, ".csv"))
  
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


# quickly choose an elbow for a PC. 
# at variance below 5% per component, choose the largest % drop
# designed for variance percentages, but will also work given a full set of Evalues
#' Quickly estimate the 'elbow' of a scree plot (PCA)
#' 
#' This function uses a rough algorithm to estimate a sensible 'elbow' to
#' choose for a PCA scree plot of eigenvalues. The function looks at an initial arbitrarily 'low'
#' level of variance and looks for the first eigenvalue lower than this. If the very first eigenvalue
#' is actually lower than this (i.e, when the PCs are not very explanatory) then this 'low' value is
#' iteratively halved until this is no longer the case. After starting below this arbitrary threshold
#' the drop in variance explained by each pair of consecutive PCs is standardized by dividing over the 
#' larger of the pair. The largest percentage drop in the series below 'low' % is selected as the 'elbow'.
#' @param varpc numeric, vector of eigenvalues, or 'percentage of variance' explained datapoints for
#'  each principle component. If only using a partial set of components, should first pass to 
#'  estimate.eig.vpcs() to estimate any missing eigenvalues.
#' @param low numeric, between zero and one, the threshold to define that a principle component
#'  does not explain much 'of the variance'.
#' @param max.pc maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n')
#' @return The number of last principle component to keep, prior to the determined elbow cutoff
#' @export
#' @seealso \code{\link{estimate.eig.vpcs}}
#' @author Nicholas Cooper 
#' @examples
#' # correlated data
#' mat <- sim.cor(100,50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick.elbow(eig)
#' pca.scree.plot(eig,elbow=elb.a,M=mat) 
#' elb.b <- quick.elbow(eig,low=.05) # decrease 'low' to select more components
#' pca.scree.plot(eig,elbow=elb.b,M=mat) 
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5,3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick.elbow(result2$sdev^2)
#' pca.scree.plot(eig2,elbow=elb2,M=mat2)
quick.elbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else { 
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}


######################################################
# GENERAL OBJECTS (objects used in multiple scripts) #
######################################################

#nloci 
nloci = 10000
#landscape dimensions (square)
ldim = 100
#number of points to sample
npts <- c(36, 81, 144, 225)
#sampling strategies
sampstrats <- c("rand", "grid", "trans", "envgeo")
#Create dataframe with all variable combos
params <- expand.grid(K = c(1,2), 
                      phi = c(0.1, 0.5),
                      m = c(0.25, 1.0),
                      seed = c(1, 2, 3),
                      H = c(0.05 , 0.5),
                      r = c(0.3, 0.6),
                      it = 0)
#TEMPORARY it = 0
