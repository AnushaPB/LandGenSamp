
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

#for scaling genetic distances from 0 to 1 for GDM
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

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
    #Calculate PC distance based on  PCs (MODIFY to make based on % var explained or tracy widom test)
    #use npcs based on sample size
    #npcs <- round(nrow(gen)*0.5,0)
    pc_dist <- as.matrix(dist(pc$x[,1:npcs], diag = TRUE, upper = TRUE))
    gendist <- range01(pc_dist)
  } else if(distmeasure == "euc"){
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  } else {
    print("appropriate gen dist measure not specified, defaulting to euclidean")
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  }
  
  #Format gdm dataframe
  site <- 1:nrow(gendist) #vector of sites
  gdmGen <- cbind(site, gendist) #bind vector of sites with gen distances
  gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, env1 = gsd_df$env1, env2 = gsd_df$env2)
  gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
  
  #scale distance from 01
  #!THINK THIS THROUGH!
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