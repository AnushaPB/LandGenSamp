

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
#number of sites
nsites <- c(9, 16, 25)
#sampling strategies
sampstrats <- c("rand", "equi", "envgeo")
#landscape dimensions (square)
ldim = 40

#Create dataframe with all variable combos
params <- expand.grid(K = c(1,2), 
                      phi = c(0.1, 0.5),
                      m = c(0.25, 1.0),
                      seed = c(1, 2, 3),
                      H = c(0.05 , 0.5),
                      r = c(0.3, 0.6),
                      it = 0:9)

#TESTING PARAMS (REMOVE LATER)
#params <- expand.grid(K = c(2, 4), 
                      #phi = c(0.1, 0.5),
                     # m = c(0.25, 1),
                      #seed = c(1, 2, 3),
                      #H = c(0.05, 0.5),
                     # r = c(0.30, 0.60),
                     # it = 1)

