
run_gdm_params <- function(i, params, path, mode = "ind"){
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
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_gdm(gen_2k, gsd_df_2k, distmeasure = "euc")
    result <- data.frame(params[i,], sampstrat = "full", nsamp = 2000, full_result, env1_err = NA, env2_err = NA, geo_err = NA)
    
    #write full datafile (temp)
    csv_file <- paste0(path, paramset,".csv")
    write.csv(result, csv_file, row.names = FALSE)
    
    #run subset data
    sampcombos <- expand.grid(sampstrats, npts)
    sub_result <- map_dfr(sampcombos, run_sub_gdm, i, params, gen, gsd_df, mode)
    
    #bind results
    result <- rbind.data.frame(result, sub_result)
    
    #export data to csv (temp)
    csv_df <- read.csv(csv_file)
    csv_df <- rbind(csv_df, sub_result)
    write.csv(csv_df, csv_file, row.names = FALSE)
    
  }
  
  return(result)
}

run_gdm <- function(gen, gsd_df, distmeasure = "euc"){
  
  if(distmeasure == "bray"){
    K <- nrow(gen)
    nloc <- ncol(gen)
    ret <- matrix(0,K,K)
    rownames(ret) <- colnames(ret) <- rownames(gen)
    for( i in 1:K){
      for( j in 1:i){
        if( i != j){
          ret[i,j] <- ret[j,i] <- sum(apply(gen[ c(i,j), ], 2, min)) / nloc
        }
      }
    }
    gendist <- ret
  } else if(distmeasure == "dps"){
    #DPS GENETIC DISTANCE
    gen[gen == 0] <- "11"
    gen[gen == 1] <- "12"
    gen[gen == 2] <- "22"
    
    genindobj <- df2genind(gen, ploidy=2, ncode=1)
    psh <- propShared(genindobj)
    dps <- 1 - psh
    gendist <- dps
  } else if(distmeasure == "pca"){
    #perform PCA
    pc <- prcomp(gen)
    #Calculate PC distance based on  PCs (?MODIFY?)
    #use npcs based on sample size
    npcs <- round(nrow(gen)*0.5,0)
    
    pc_dist <- as.matrix(dist(pc$x[,1:npcs], diag = TRUE, upper = TRUE))
    gendist <- pc_dist
  } else if(distmeasure == "euc"){
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  } else {
    print("appropriate gen dist measure not specified")
  }
  
  #Format gdm dataframe
  site <- 1:nrow(gendist) #vector of sites
  gdmGen <- cbind(site, gendist) #bind vector of sites with gen distances
  gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, env1 = gsd_df$env1, env2 = gsd_df$env2)
  gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
  #SCALE DISTANCE FROM 0 to 1 if max(distance) >1 (gdm only works for 0<vals<1) (?MODIFY?)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  gdmData$distance <- range01(gdmData$distance) 
  
  #run GDM
  gdm.model <- gdm(gdmData, geo = TRUE)
  
  #check var importance/significance (ASK IAN IF WE WANT TO DO THIS OR JUST COMPARE THE COEFFICIENTS FROM A FULL MODEL (PROS: FASTER/EASIER))
  #system.time(vars <- gdm.varImp(gdmData, geo = TRUE, splines = NULL, nPerm=50))
  
  if(is.null(gdm.model)){
    #turn results into dataframe
    results <- data.frame(env1_coeff = "NULL",
                          env2_coeff = "NULL",
                          geo_coeff = "NULL"
    )
  } else {
    predictors <- coeffs(gdm.model)
    predictors
    
    #turn results into dataframe
    results <- data.frame(env1_coeff = predictors[predictors$predictor == "env1", "coefficient"],
                          env2_coeff = predictors[predictors$predictor == "env2", "coefficient"],
                          geo_coeff = predictors[predictors$predictor == "Geographic", "coefficient"]
    )
  }
  
  
  #remove rownames
  rownames(results) <- NULL
  
  return(results)
}

run_sub_gdm <- function(sampcombo, i, params, gen, gsd_df, full_result, mode = "ind"){
  
  # get nsamp and sampstrat
  nsamp <- sampcombos[1]
  sampstrat <- sampcombos[2]
  
  #subsample from data based on sampling strategy and number of samples
  subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
  subgen <- gen[subIDs,]
  subgsd_df <- gsd_df[subIDs,]
  
  #run analysis using subsample
  sub_result <- run_gdm(subgen, subgsd_df, distmeasure = "euc")
  
  #calculate err
  env1_err <- err_coeff(full_result$env1_coeff, sub_result$env1_coeff)
  env2_err <- err_coeff(full_result$env2_coeff, sub_result$env2_coeff)
  geo_err <- err_coeff(full_result$geo_coeff, sub_result$geo_coeff)
  
  #save and format new result
  sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsamp, sub_result, 
                           env1_err = env1_err, env2_err = env2_err, geo_err = geo_err)
  
  #bind results
  return(sub_result)
}

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
