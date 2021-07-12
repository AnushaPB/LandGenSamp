library("here") #paths
library("gdm") #GDM

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")

#set seed
set.seed(42)

###########
#   GDM   #
###########
# Sum coefficients for each predictor (each has 3 splines)
coeffs <- function(gdm.model){
  coefSums <- c()
  for (i in 1:length(gdm.model$predictors)){
    j <- (i * 3) - 2
    coefSums[i] <- sum(gdm.model$coefficients[j:(j+2)])
  }
  
  # Add those values to a simple data frame
  coeffs <- data.frame(predictor = gdm.model$predictors, coefficient = coefSums)
  return(coeffs)
}


run_gdm <- function(gen_filepath, gsd_filepath){
  
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #CALCULATING PC BASED GENETIC DISTANCE
  #convert gen data to matrix
  gen <- as.matrix(gen)
  #perform PCA
  pc <- prcomp(gen)
  #Calculate PC distance based on first three PCs (?MODIFY?)
  pc_dist <- as.matrix(dist(pc$x[,1:100], diag = TRUE, upper = TRUE))
  
  #Format gdm dataframe
  site <- 1:nrow(pc_dist) #vector of sites
  gdmGen <- cbind(site, pc_dist) #bind vector of sites with gen distances
  gdmPred <- data.frame(site = site, Longitude = gea_df$x, Latitude = gea_df$y, env1 = gea_df$env1, env2 = gea_df$env2)
  gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
  #SCALE DISTANCE FROM 0 to 1 if max(distance) >1 (gdm only works for 0<vals<1) (?MODIFY?)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  gdmData$distance <- range01(gdmData$distance) 
  
  #run GDM
  gdm.model <- gdm(gdmData, geo = TRUE)
  
  #check var importance/significance (ASK IAN IF WE WANT TO DO THIS OR JUST COMPARE THE COEFFICIENTS FROM A FULL MODEL (PROS: FASTER/EASIER))
  #vars <- gdm.varImp(gdmData, geo = TRUE, splines = NULL, nPerm=100)
  

  predictors <- coeffs(gdm.model)
  predictors

  #turn results into dataframe
  results <- data.frame(env1_coeff = predictors[predictors$predictor == "env1", "coefficient"],
                        env2_coeff = predictors[predictors$predictor == "env2", "coefficient"],
                        geo_coeff = predictors[predictors$predictor == "Geographic", "coefficient"]
                        )
  #remove rownames
  rownames(results) <- NULL
  
  return(results)
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


res_mmrr <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library(here)
  
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_gdm(gen_filepath = gen_filepath,
                      gsd_filepath = gsd_filepath)
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_mmrr)
write.csv(stats_out, "gdm_results.csv")



