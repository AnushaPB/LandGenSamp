#####################
# GENERAL FUNCTIONS #
#####################
# create filepath based on params index and data type (e.g. genetic data = gen/dos, geospatial data = gsd, and adaptive loci = loci)
create_filepath <- function(i, params, type) {
  # set of parameter names in filepath form
  paramset <- paste0(
    "K", params[i, "K"],
    "_phi", params[i, "phi"] * 100,
    "_m", params[i, "m"] * 100,
    "_seed", params[i, "seed"],
    "_H", params[i, "H"] * 100,
    "_r", params[i, "r"] * 100
  )

  # Create data directory path
  datadir = here::here("p1_gnxsims", "gnx", "LGS_data")

  # Create model directory path
  moddir <- here::here(datadir, paste0("GNX_mod-", paramset), paste0("it-", params[i, "it"]), "spp-spp_0")

  # Different file patterns for different data types
  if (type == "gen") {
    filepath <- paste0("/mod-", paramset, "_it-", params[i, "it"], "_t-1000_spp-spp_0.vcf")
  }
  if (type == "gsd") {
    filepath <- paste0("/mod-", paramset, "_it-", params[i, "it"], "_t-1000_spp-spp_0.csv")
  }
  if (type == "dos") {
    filepath <- paste0("/dos-", paramset, "_it-", params[i, "it"], "_t-1000_spp-spp_0.csv")
  }

  # return full filepath
  filepath <- here::here(moddir, filepath)
  print(filepath)

  return(filepath)
}

# Get gen data
get_gen <- function(filepath) {
  # read vcf
  vcf <- vcfR::read.vcfR(filepath)
  # convert to genlight from vcf
  genmat <- vcf_to_dosage(vcf)
  return(genmat)
}

# Convert from vcf to dosage
vcf_to_dosage <- function(vcf) {
  # convert to genlight from vcf
  genlight <- vcfR::vcfR2genlight(vcf)
  # convert to matrix
  genmat <- as.matrix(genlight)
  # assign IDs from genlight to matrix rownames
  rownames(genmat) <- genlight@ind.names
  return(genmat)
}

#  Get geospatial data
get_gsd <- function(filepath) {
  gsd_df <- read.csv(filepath)
  # assign IDs to rownames
  rownames(gsd_df) <- gsd_df$idx
  # extract env and z values
  gsd_df <- gsd_df %>%
    tidyr::extract(z, into = c("z1", "z2"), regex = "\\[(.*), (.*)\\]", remove = FALSE) %>%
    tidyr::extract(e, into = c("env0", "env1", "env2"), regex = "\\[(.*), (.*), (.*)\\]", remove = FALSE) %>%
    dplyr::mutate(across(c(z1, z2, env0, env1, env2), as.numeric))
  
  # if there is only one trait value get the number from the z column
  # (not relevent for our main simulations)
  if (all(is.na(gsd_df$z1))) gsd_df <- gsd_df %>% mutate(z1 = parse_number(z))

  # correct coordinates
  gsd_df$y <- -gsd_df$y

  return(gsd_df)
}

get_dos <- function(filepath){
  read.csv(filepath, row.names = 1)
}

# general function to get data
get_data <- function(i, params, type) {
  # Create file path
  filepath <- create_filepath(i, params, type)

  # different file patterns for different data types
  if (type == "gen") return(get_gen(filepath))

  if (type == "gsd") return(get_gsd(filepath))
  
  if (type == "dos") return(get_dos(filepath))
  
}

# get list of sampling IDs that correspond with parameter set, sampling strategy, and number of samples
get_samples <- function(param_set, sampstrat, nsamp, outdir = here("p2_sampling", "outputs"), site = FALSE) {
  # param_set - vector of one set of parameters (e.g. params[i,])
  # sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  # nsamp - number of samples

  if (site) {
    subIDs <- read.csv(paste0(outdir, "/site_samples_", sampstrat, nsamp, ".csv"))
  } else {
    subIDs <- read.csv(paste0(outdir, "/samples_", sampstrat, nsamp, ".csv"))
  }

  subIDs <- subIDs[subIDs$K == param_set$K &
    subIDs$phi == param_set$phi &
    subIDs$m == param_set$m &
    subIDs$seed == param_set$seed &
    subIDs$H == param_set$H &
    subIDs$r == param_set$r &
    subIDs$it == param_set$it, ]

  # confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)

  # remove parameter columnds and convert to vector of IDs
  subIDs <- subIDs[, !names(subIDs) %in% names(param_set)]
  subIDs <- unlist(subIDs)

  # confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))

  return(as.character(subIDs))
}

# check if files for params exist and skip otherwise
skip_check <- function(i, params) {
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  skip_to_next <- FALSE
  if (file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE) {
    skip_to_next <- TRUE
  }
  if (skip_to_next) {
    print("File does not exist:")
    print(params[i, ])
  }
  return(skip_to_next)
}

# function to calculate RMSE
err_coeff <- function(full_coeff, sub_coeff) {
  err <- (sub_coeff - full_coeff)
  return(err)
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
#' mat <- sim.cor(100, 50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick.elbow(eig)
#' pca.scree.plot(eig, elbow = elb.a, M = mat)
#' elb.b <- quick.elbow(eig, low = .05) # decrease 'low' to select more components
#' pca.scree.plot(eig, elbow = elb.b, M = mat)
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5, 3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick.elbow(result2$sdev^2)
#' pca.scree.plot(eig2, elbow = elb2, M = mat2)
quick.elbow <- function(varpc, low = .08, max.pc = .9) {
  ee <- varpc / sum(varpc) # ensure sums to 1
  # print(round(log(ee),3))
  while (low >= max(ee)) {
    low <- low / 2
  } # when no big components, then adjust 'low'
  lowie <- (ee < low)
  highie <- ee > low / 8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if (length(low.ones) > 0) {
    if (length(low.ones) == 1) {
      elbow <- low.ones
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set)) / (set[1:(length(set) - 1)])
      infz <- is.infinite(pc.drops)
      # print(pc.drops)
      elbow <- which(pc.drops == max(pc.drops[!infz], na.rm = T))[1] + others
    }
  } else {
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee)
  }
  if (tail(cumsum(ee[1:elbow]), 1) > max.pc) {
    elbow <- which(cumsum(ee) > max.pc)[1] - 1
  }
  if (elbow < 1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}


# get list of site IDs that correspond with parameter set, sampling strategy, and number of samples (and sample IDs)
get_sites <- function(param_set, sampstrat, nsamp, dir = here("p2_sampling", "outputs")) {
  # param_set - vector of one set of parameters (e.g. params[i,])
  # params - full set of parameters
  # sampstrat - sampling strategy (e.g. "rand", "grid", "trans", "envgeo")
  # nsamp - number of samples

  subIDs <- read.csv(paste0(dir, "/site_ids_", sampstrat, nsamp, ".csv"))

  subIDs <- subIDs[subIDs$K == param_set$K &
    subIDs$phi == param_set$phi &
    subIDs$m == param_set$m &
    subIDs$seed == param_set$seed &
    subIDs$H == param_set$H &
    subIDs$r == param_set$r &
    subIDs$it == param_set$it, ]

  # confirm there is only one set of IDs being used
  stopifnot(nrow(subIDs) == 1)

  # remove parameter columns and convert to vector of IDs
  subIDs <- subIDs[, !names(subIDs) %in% names(param_set)]
  subIDs <- unlist(subIDs)

  # confirm that final set of IDs is a vector
  stopifnot(is.vector(subIDs))

  return(subIDs)
}

# make loci_df and convert from python to R indexing by adding 1
get_loci <- function() data.frame(trait1 = 0:3, trait2 = 4:7) + 1

# utility function to return 0 instead of integer(0) for which
which0 <- function(x) {
  y <- which(x)
  if (length(y) == 0) {
    return(0)
  } else {
    return(y)
  }
}

# get list of packages to provide when parallelizing
get_packages <- function() {
  c("here", "vcfR", "adegenet", "stringr", "dplyr", "tidyr", "purrr", "lfmm", "AssocTests", "gdm", "vegan", "robust", "qvalue", "raster", "hierfstat")
}

######################################################
# GENERAL OBJECTS (objects used in multiple scripts) #
######################################################
# landscape dimensions (square)
ldim <- 100
# sampling strategies
sampstrats <- c("rand", "grid", "trans", "envgeo")
sitestrats <- c("rand", "equi", "envgeo")
nsites <- c(9, 16, 25)
nsamps <- c(36, 81, 144, 225)
# Create dataframe with all variable combos
params <- expand.grid(
  K = c(1, 2),
  phi = c(0.5, 1.0),
  m = c(0.25, 1.0),
  seed = c(1, 2, 3),
  H = c(0.05, 0.5),
  r = c(0.3, 0.6),
  it = 0:9
)
