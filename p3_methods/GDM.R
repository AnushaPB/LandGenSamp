library("here") #paths
library("gdm") #GDM
library("vcfR")
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


run_gdm <- function(gen, gsd_df){
  
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
  gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, env1 = gsd_df$env1, env2 = gsd_df$env2)
  gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
  #SCALE DISTANCE FROM 0 to 1 if max(distance) >1 (gdm only works for 0<vals<1) (?MODIFY?)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  gdmData$distance <- range01(gdmData$distance) 
  
  #run GDM
  gdm.model <- gdm(gdmData, geo = TRUE)
  
  #check var importance/significance (ASK IAN IF WE WANT TO DO THIS OR JUST COMPARE THE COEFFICIENTS FROM A FULL MODEL (PROS: FASTER/EASIER))
  system.time(vars <- gdm.varImp(gdmData, geo = TRUE, splines = NULL, nPerm=50))
  

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


res_gdm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  
  skip_to_next <- FALSE
  if(file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run gdm
  if(skip_to_next == FALSE){
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    
    #subsample full data randomly
    s <- sample(2000, nrow(gsd_df), replace = FALSE)
    gen <- gen[s,]
    gsd_df <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_gdm(gen, gsd_df, loci_df)
    result <- data.frame(sampstrat = "full", nsamp = nrow(gsd_df), full_result, env1_rmse = NA, env2_rmse = NA, geo_rmse = NA)
    
    #write full datafile (temp)
    csv_file <- paste0("gdm_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_gdm(subgen, subgsd_df)
        
        #calculate RMSE
        env1_rmse <- rmse_coeff(full_result$env1_coeff, sub_result$env1_coeff)
        env2_rmse <- rmse_coeff(full_result$env2_coeff, sub_result$env2_coeff)
        geo_rmse <- rmse_coeff(full_result$geo_coeff, sub_result$geo_coeff)
        
        #save and format new result
        sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, sub_result, 
                                env1_rmse = env1_rmse, env2_rmse = env2_rmse, geo_rmse = geo_rmse)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  return(result)
  
  gc()
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_mmrr)
write.csv(stats_out, "gdm_results.csv")



