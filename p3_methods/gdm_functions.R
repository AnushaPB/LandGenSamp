
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
  #Format data for GDM  
  gendist <- calc_dist(gen, distmeasure)
  
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
  
  if(is.null(gdm.model)){
    #turn results into dataframe
    results <- data.frame(env1_coeff = "NULL",
                          env2_coeff = "NULL",
                          geo_coeff = "NULL",
                          env1_p = NA,
                          env2_p = NA,
                          geo_p = NA)
  } else {
    predictors <- coeffs(gdm.model)
    predictors
    
    # turn results into dataframe
    results <- data.frame(env1_coeff = predictors[predictors$predictor == "env1", "coefficient"],
                          env2_coeff = predictors[predictors$predictor == "env2", "coefficient"],
                          geo_coeff = predictors[predictors$predictor == "Geographic", "coefficient"])
    
    # get pvalues
    safe_gdm.varImp <- safely(gdm.varImp)
    modTest <- safe_gdm.varImp(gdmData, geo=T, nPerm=50, parallel=F, predSelect=F)
    if (is.null(modTest$error)) {
      pvals <- modTest$result$`Predictor p-values`
      pvals$var <- row.names(pvals)
      pvals <- left_join(data.frame(var = c("Geographic", "env1", "env2")), pvals)
      results <- data.frame(results,
                 env1_p = pvals[pvals$var == "env1", 2],
                 env2_p = pvals[pvals$var == "env2", 2],
                 geo_p = pvals[pvals$var == "Geographic", 2])
    } else {
      results <- data.frame(results,
                            env1_p = NA,
                            env2_p = NA,
                            geo_p = NA)
    }
    
  }
  
  
  #remove rownames
  rownames(results) <- NULL
  
  return(results)
}
