
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


# run GDM with each environmental variable treated separately
run_gdm2 <- function(gendist, gsd_df){
  #Format gdm dataframe
  site <- 1:nrow(gendist) #vector of sites
  gdmGen <- cbind(site, gendist) #bind vector of sites with gen distances
  
  # model combo
  gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, env1 = gsd_df$env1, env2 = gsd_df$env2)
  
  gdmData <-
    formatsitepair(
      gdmGen,
      bioFormat = 3,
      predData = gdmPred,
      XColumn = "Longitude",
      YColumn = "Latitude",
      siteColumn = "site"
    )
  
  #scale distance from 01
  gdmData$distance <- gdmData$distance/100
  
  #run GDM
  gdm.model <- gdm(gdmData, geo = TRUE)
  
  if (is.null(gdm.model)){
    #turn results into dataframe
    results <- data.frame(env1_coeff = NA,
                          env2_coeff = NA,
                          geo_coeff = NA,
                          ratio = NA,
                          env1_p = NA,
                          env2_p = NA,
                          geo_p = NA)
  } else {
    predictors <- coeffs(gdm.model)
    
    # turn results into dataframe
    results <- data.frame(env1_coeff = predictors[predictors$predictor == "env1", "coefficient"],
                          env2_coeff = predictors[predictors$predictor == "env2", "coefficient"],
                          geo_coeff = predictors[predictors$predictor == "Geographic", "coefficient"])
    
    results$ratio <- sum(abs(results$env1_coeff) + abs(results$env2_coeff))/abs(results$geo_coeff)
    
    
    # check if both env coefficients are 0
    # modTest cannot be run if both are 0
    do_modTest <- !(results$env1_coeff == 0 & results$env2_coeff == 0)
      
    # get pvalues
    if (do_modTest) modTest <- gdm.varImp_custom(gdmData, geo = TRUE, nPerm = 50, parallel = F, predSelect = F) else modTest <- NULL
    
    if (!is.null(modTest)) {
      pvals <- modTest$`Predictor p-values`
      pvals$var <- row.names(pvals)
      pvals <- left_join(data.frame(var = c("env1", "env2", "Geographic")), pvals, by = "var")
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

# Customized variable importance function (changed code so that if null models pop up they don't result in an error, they just aren't counted)
# 
# All changes are marked by "CHANGE"
#
# *Summary of issue*
# Problem: In the loop to create "permVarDev" if there is a NULL gdm model for a variable (k) the value assigned to permVarDev[[k]] is NULL which means that the list ends up having a length less than the number of variables such that the line:
# names(permVarDev) <- varNames.x throws the error: 'names' attribute [2] must be the same length as the vector [1]".
# Change: I think one solution to this is to have a value of NA assigned to permVarDev[[k]] if the GDM model is NULL. 
gdm.varImp_custom <- function(spTable, geo, splines=NULL, knots=NULL, predSelect=FALSE,
                       nPerm=50, pValue=0.05, parallel=FALSE, cores=2, sampleSites=1,
                       sampleSitePairs=1, outFile=NULL){
  
  ##assign k to prevent issues with cran checking
  k <- NULL
  
  ##error checking for input objects
  ##checks to see if in site-pair format from formatsitepair function
  if(!is(spTable, "gdmData")){
    warning("The spTable object is not of class 'gdmData'. See the formatsitepair function for help.")
  }
  ##checks to makes sure data is a matrix or data frame
  if(!(is(spTable, "gdmData") | is(spTable, "matrix") | is(spTable, "data.frame"))){
    stop("spTable argument needs to be of class 'gdmData', 'matrix', or 'data frame'")
  }
  
  ##sanity check on the data table
  if(ncol(spTable) < 6){
    stop("spTable object requires at least 6 columns: distance, weights, s1.xCoord, s1.yCoord, s2.xCoord, s2.yCoord")
  }
  if(nrow(spTable) < 1){
    stop("The spTable object contains zero rows of data.")
  }
  
  ##checks that geo has either TRUE or FALSE
  if(!(geo==TRUE | geo==FALSE)){
    stop("The geo argument must be either TRUE or FALSE.")
  }
  ##makes sure splines is a numeric vector
  if(is.null(splines)==FALSE & !is(splines, "numeric")){
    stop("The splines argument needs to be a numeric data type.")
  }
  ##checks knots inputs
  if(is.null(knots)==FALSE & !is(knots, "numeric")){
    stop("The knots argument needs to be a numeric data type.")
  }
  ##checks that predSelect has either TRUE or FALSE
  if(!(predSelect==TRUE | predSelect==FALSE)){
    stop("The predSelect argument must be either TRUE or FALSE.")
  }
  ##makes sure that nPerm is a positive integer
  if((is.null(nPerm)==FALSE & is.numeric(nPerm)==FALSE) | nPerm<1){
    stop("The nPerm argument needs to be a positive integer.")
  }
  ##checks that parallel has either TRUE or FALSE
  if(!(parallel==TRUE | parallel==FALSE)){
    stop("The parallel argument must be either TRUE or FALSE.")
  }
  ##makes sure that cores has a value when parallel is true
  if(parallel==TRUE & is.null(cores)==TRUE){
    stop("If parallel==TRUE, the number of cores must be specified.")
  }
  ##makes sure that cores is a positive integer
  if((is.null(cores)==FALSE & is.numeric(cores)==FALSE) | cores<1){
    stop("The cores argument needs to be a positive integer.")
  }
  ##makes sure that both sampleSites and sampleSitePairs are a number between 0 and 1,
  ##and that neither is equal to 0
  if(is.numeric(sampleSites)==FALSE | sampleSites<0 | sampleSites>1){
    stop("The sampleSites argument needs to be a positive number between 0 and 1.")
  }
  if(is.numeric(sampleSitePairs)==FALSE | sampleSitePairs<0 | sampleSitePairs>1){
    stop("The sampleSitePairs argument needs to be a positive number between 0 and 1.")
  }
  if(sampleSites==0){
    stop("A sampleSites value of 0 will remove all sites from the analysis.")
  }
  if(sampleSitePairs==0){
    stop("A sampleSitePairs value of 0 will remove all sites from the analysis.")
  }
  ##checks to see if the user has requested for an output file to be written, and if so
  ##makes sure that it is formatted correctly
  if(is.null(outFile)==FALSE){
    ##first makes sure outFile is a string
    if(is.character(outFile)==FALSE){
      stop("The outFile argument needs to be a character string of the directory and file name you wish the tables to be written to")
    }
    ##makes sure that text has ".RData" in it, if not, adds it
    outFileChar <- nchar(outFile)
    if(substr(outFile, outFileChar-5, outFileChar)!=".RData"){
      outFile <- paste(outFile, ".RData", sep="")
    }
    ##checks to see if there is a path as well as a file name
    if(length(strsplit(outFile,"/")[[1]])>1){
      splitOutFile <- strsplit(outFile,"/")[[1]][-length(strsplit(outFile,"/")[[1]])]
      dir.create(paste(splitOutFile, collapse="/"))
    }else{
      outFile <- paste("./", outFile, sep="")
    }
  }
  
  ##double makes sure these values are integers, seems to truncate if not
  nPerm <- as.integer(nPerm)
  cores <- as.integer(cores)
  
  ##removes a user specified number of sites from the site-pair table
  if(sampleSites<1){
    spTable <- subsample.sitepair(spTable, sampleSites=sampleSites)
    ##throws warning if sampleSitePairs<1 as well
    if(sampleSitePairs<1){
      warning("You have selected to randomly remove sites and/or site-pairs.")
    }
  }
  ##removes a user specified number of site-pairs from the site-pair table
  if(sampleSitePairs<1){
    ##determine which rows to remove
    numRm <- sample(1:nrow(spTable), round(nrow(spTable)*(1-sampleSitePairs)))
    spTable <- spTable[-c(numRm),]
  }
  
  ##check that the response data is [0..1]
  rtmp <- spTable[,1]
  if(length(rtmp[rtmp<0]) > 0){
    stop("The spTable contains negative distance values. Must be between 0 - 1.")
  }
  if (length(rtmp[rtmp>1]) > 0){
    stop("The spTable contains distance values greater than 1. Must be between 0 - 1.")
  }
  
  # number of variables in the site-pair table, adds 1 if geo=TRUE
  nVars <- (ncol(spTable)-6)/2
  # create vector of variable names
  varNames <- colnames(spTable[c(7:(6+nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if(geo==TRUE){
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  
  if(nVars<2){
    stop("Function requires at least two predictor variables.")
  }
  
  # run initial GDM to see if any vars have zero I-spline coeffs
  message(paste0("Fitting initial model with all ", nVars,  " predictors..."))
  Sys.sleep(0.5)
  fullGDM <- gdm(spTable, geo=geo, splines=splines, knots=knots)
  
  # check for zero coeffs
  thiscoeff <- 1
  thisquant <- 1
  sumCoeff <- NULL
  for(i in 1:length(fullGDM$predictors)){
    numsplines <- fullGDM$splines[[i]]
    holdCoeff <- NULL
    for(j in 1:numsplines){
      holdCoeff[j] <- fullGDM$coefficients[[thiscoeff]]
      thiscoeff <- thiscoeff + 1
    }
    sumCoeff[i] <- sum(holdCoeff)
  }
  
  # remove any predictors with sumCoeff=0
  zeroSum <- fullGDM$predictors[which(sumCoeff==0)]
  if(length(zeroSum)>0){
    for(p in 1:length(zeroSum)){
      message(paste0("Sum of I-spline coefficients for predictor ", zeroSum[p]," = 0"))
      Sys.sleep(0.5)}
    #message("\n")
    for(p in 1:length(zeroSum)){
      if(zeroSum[p]=="Geographic"){
        message("Setting Geo=FALSE and proceeding with permutation testing...")
      } else {
        message(paste0("Removing ", zeroSum[p], " and proceeding with permutation testing..."))
      }
      Sys.sleep(0.5)
    }
    
    if(length(grep("Geographic", zeroSum))==1){
      geo <- FALSE
      zeroSum <- zeroSum[-grep("Geographic", zeroSum)]
    }
    
    for(z in zeroSum){
      ##select variable columns to be removed from original site-pair table
      testVarCols1 <- grep(paste("^s1.", z, "$", sep=""), colnames(spTable))
      testVarCols2 <- grep(paste("^s2.", z, "$", sep=""), colnames(spTable))
      spTable <- spTable[,-c(testVarCols1, testVarCols2)]
    }
  }
  
  # recalculate the number of variables in the site-pair table, adds 1 if geo=TRUE
  nVars <- (ncol(spTable)-6)/2
  # create vector of variable names
  varNames <- colnames(spTable[c(7:(6+nVars))])
  varNames <- sapply(strsplit(varNames, "s1."), "[[", 2)
  if(geo==TRUE){
    nVars <- nVars + 1
    varNames <- c("Geographic", varNames)
  }
  
  # stop routine if fewer than two variables remain
  if(nVars<2){
    # CHANGE 1: changed stop() to warning() so that NULL is returned
    #stop("Function requires at least two predictor variables.")
    warning("Function requires at least two predictor variables. Returning NULL")
    return(NULL)
  }
  
  # reduce number of cores to nVars
  if(cores>nVars){
    cores <- nVars
  }
  
  # crate spline object
  splines <- rep(unique(splines), nVars)
  
  ##First create a spTable to determine the index of each site in the site-pair table
  sortMatX <- sapply(1:nrow(spTable), function(i, spTab){c(spTab[i,3], spTab[i,5])}, spTab=spTable)
  sortMatY <- sapply(1:nrow(spTable), function(i, spTab){c(spTab[i,4], spTab[i,6])}, spTab=spTable)
  sortMatNum <- sapply(1:nrow(spTable), function(i){c(1,2)})
  sortMatRow <- sapply(1:nrow(spTable), function(i){c(i,i)})
  ##adds a column of NA for index to be added to
  fullSortMat <- cbind(as.vector(sortMatX), as.vector(sortMatY), as.vector(sortMatNum), as.vector(sortMatRow), rep(NA, length(sortMatX)))
  ##assigns sites by unique coordinates
  siteByCoords <- as.data.frame(unique(fullSortMat[,1:2]))
  ##number of sites to expect by coordinates
  numSites <- nrow(siteByCoords)
  ##assigns site index based on coordinates
  for(i in 1:numSites){
    fullSortMat[which(fullSortMat[,1]==siteByCoords[i,1] & fullSortMat[,2]==siteByCoords[i,2]),5] <- i
  }
  
  ##create index table to know where each site is in input site-pair table
  indexTab <- matrix(NA,nrow(spTable),2)
  for(iRow in 1:nrow(fullSortMat)){
    indexTab[fullSortMat[iRow,4],fullSortMat[iRow,3]] <- fullSortMat[iRow,5]
  }
  
  ## And remove the sorting table and supporting objects to free up memory
  rm(fullSortMat)
  rm(sortMatX)
  rm(sortMatY)
  rm(sortMatNum)
  rm(sortMatRow)
  rm(siteByCoords)
  
  ##create site x predictor table, to be able to rebuild site-pair table later in function
  exBySite <- lapply(1:numSites, function(i, index, tab){
    rowSites <- which(index[,1] %in% i)
    if(length(rowSites)<1){
      rowSites <- which(index[,2] %in% i)
    }
    exSiteData <- tab[rowSites[1],]
    return(exSiteData)
  }, index=indexTab, tab=spTable)
  ##identifies the one site not in the first column of the index table
  outSite <- which(!(1:numSites %in% indexTab[,1]))
  
  ##sets up siteXvar table, uses for loop to make sure have steps correct
  for(i in 1:length(exBySite)){
    ##grabs row and identify if should take s1 or s2 by rather or not number appears in outsite
    siteRow <- exBySite[[i]]
    if(i %in% outSite){
      ##extracts the data from the site-pair table by site
      siteRow <- siteRow[grep("s2.", colnames(siteRow))]
      colnames(siteRow) <- sapply(strsplit(colnames(siteRow), "s2."), "[[", 2)
    }else{
      ##extracts the data from the site-pair table by site
      siteRow <- siteRow[grep("s1.", colnames(siteRow))]
      colnames(siteRow) <- sapply(strsplit(colnames(siteRow), "s1."), "[[", 2)
    }
    exBySite[[i]] <- siteRow
  }
  
  ##transforms data from list to table
  siteData <- do.call("rbind", exBySite)
  
  ##sets up objects to be returned by the function
  modelTestValues <- matrix(NA,4,nVars*5,dimnames = list(c("Model deviance",
                                                           "Percent deviance explained",
                                                           "Model p-value",
                                                           "Fitted permutations"),
                                                         c("All predictors",
                                                           "1-removed",
                                                           paste(seq(2,nVars*5-1), "-removed", sep=""))))
  ##deviance reduction predictor table
  varImpTable <- matrix(NA, nVars, nVars*5-1)
  rownames(varImpTable) <- varNames
  colnames(varImpTable) <- c("All predictors",
                             "1-removed",
                             paste(seq(2,nVars*5-2), "-removed", sep=""))
  ##p value predictor table
  pValues <- nModsConverge <- varImpTable
  
  ##assigns given site-pair table to new variable, to prevent changing the original input
  currSitePair <- spTable
  nullGDMFullFit <- 0  ##a variable to track if the fully fitted gdm model returned a NULL object
  
  # create the set up permuted site-pair tables to be used for all
  # downstream analyses
  message(paste0("Creating ", nPerm, " permuted site-pair tables..."))
  
  # use replicate if number of perms is small or parallel == FALSE
  if(parallel==F | nPerm <= 25){
    # CHANGE 2 (trivial): added "pbapply::" instead of importing
    permSpt <- pbapply::pbreplicate(nPerm, list(permutateSitePair(currSitePair,
                                                                  siteData,
                                                                  indexTab,
                                                                  varNames)))}
  # create permuted tables in parallel otherwise
  if(parallel==T & nPerm > 25){
    # set up parallel processing
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    permSpt <- foreach(k=1:nPerm,
                       .verbose=F,
                       .packages=c("gdm"),
                       .export=c("permutateSitePair")) %dopar%
      permutateSitePair(currSitePair, siteData, indexTab, varNames)
    stopCluster(cl)}
  
  # create new vector to track varNames
  varNames.x <- varNames
  message("Starting model assessment...")
  for(v in 1:length(varNames)){
    # ends the loop if only 1 variable remains
    if(length(varNames.x)<2){
      stop("Only one predictor remains...variable assessment stopped.")
    }
    
    if(v>1){
      message(paste0("Removing ", names(elimVar), " and proceeding with the next round of permutations."))
    }
    
    if(is.numeric(splines)){
      splines <- rep(unique(splines), length(varNames.x))
    }
    
    # runs gdm, first time on site-pair table with all variables
    # As variables are eliminated the "full" site-pair table will
    # contain fewer variables
    fullGDM <- gdm(currSitePair, geo=geo, splines=splines, knots=knots)
    message(paste0("Percent deviance explained by the full model =  ", round(fullGDM$explained,3)))
    
    if(is.null(fullGDM)==TRUE){
      warning(paste("The model did not converge when testing variable: ", varNames.x[v],
                    ". Terminating analysis and returning output completed up to this point.", sep=""))
      break
    }
    
    # fit gdm to each permuted site-pair table
    message("Fitting GDMs to the permuted site-pair tables...")
    permGDM <- lapply(permSpt, function(x){
      gdm(x, geo=geo, splines=splines, knots=knots)
    })
    
    ##extracts deviance of permuted gdms
    permModelDev <- sapply(permGDM, function(mod){mod$gdmdeviance})
    ##if needed, removes nulls from output permModelDev
    modPerms <- length(which(sapply(permModelDev,is.null)==TRUE))
    if(modPerms>0){
      permModelDev <- unlist(permModelDev[-(which(sapply(permModelDev,is.null)==T))])
    }
    
    if(v>1){
      colnames(modelTestValues)[v] <- colnames(pValues)[v] <- colnames(varImpTable)[v] <- colnames(nModsConverge)[v] <- paste0(names(elimVar),"-removed")
    }
    
    ##begins to fill in the output table with data from fully fitted model
    modelTestValues[1,v] <- round(fullGDM$gdmdeviance,3)
    modelTestValues[2,v] <- round(fullGDM$explained,3)
    #p-value
    modelTestValues[3,v] <- round(sum(permModelDev<=fullGDM$gdmdeviance)/(nPerm-modPerms),3)
    ##fitted permutations
    modelTestValues[4,v] <- round(nPerm-modPerms,0)
    
    # now permutate each variable individually and fit models
    if(parallel==TRUE){
      if(length(varNames.x)<cores){
        cores <- length(varNames.x)
      }
      # set up parallel processing
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      
      # foreach function to create site-pair tables with each variable permuted,
      # fit gdms and extract deviance.
      permVarDev <- foreach(k=1:length(varNames.x),
                            .verbose=F,
                            .packages=c("gdm"),
                            .export = c("currSitePair")) %dopar%{
                              if(varNames.x[k]!="Geographic"){
                                # permute a single variable
                                lll <- lapply(permSpt, function(x, spt=currSitePair){
                                  idx1 <- grep(paste("^s1.", varNames.x[k], "$", sep=""), colnames(x))
                                  idx2 <- grep(paste("^s2.", varNames.x[k], "$", sep=""), colnames(x))
                                  spt[,c(idx1, idx2)] <- x[,c(idx1, idx2)]
                                  return(spt)
                                })}
                              
                              if(varNames.x[k]=="Geographic"){
                                # permute a single variable
                                lll <- lapply(permSpt, function(x, spt=currSitePair){
                                  s1 <- sample(1:nrow(spt), nrow(spt))
                                  s2 <- sample(1:nrow(spt), nrow(spt))
                                  s3 <- sample(1:nrow(spt), nrow(spt))
                                  s4 <- sample(1:nrow(spt), nrow(spt))
                                  spt[,3] <- spt[s1,3]
                                  spt[,4] <- spt[s2,4]
                                  spt[,5] <- spt[s3,5]
                                  spt[,6] <- spt[s4,6]
                                  return(spt)
                                })}
                              
                              gdmPermVar <- lapply(lll, function(x){
                                try(gdm(x, geo=geo, splines=splines, knots=knots))
                              })
                              
                              ##extracts deviance of permuted gdms
                              permModelDev <- sapply(gdmPermVar, function(mod){mod$gdmdeviance})
                              return(permModelDev)
                            }
      
      ##closes cores
      #close(pb)
      stopCluster(cl)
    }
    
    if(parallel==FALSE){
      # for-loop to create site-pair tables with each variable permuted,
      # fit gdms and extract deviance.
      permVarDev <- list()
      for(k in 1:length(varNames.x)){
        if(varNames.x[k]!="Geographic"){
          message(paste0("Assessing importance of ", varNames.x[k], "..."))
          # permute a single variable
          lll <- lapply(permSpt, function(x, spt=currSitePair){
            idx1 <- grep(paste("^s1.", varNames.x[k], "$", sep=""), colnames(x))
            idx2 <- grep(paste("^s2.", varNames.x[k], "$", sep=""), colnames(x))
            spt[,c(idx1, idx2)] <- x[,c(idx1, idx2)]
            return(spt)
          })}
        
        if(varNames.x[k]=="Geographic"){
          message("Assessing importance of geographic distance...")
          # permute a single variable
          lll <- lapply(permSpt, function(x, spt=currSitePair){
            s1 <- sample(1:nrow(spt), nrow(spt))
            s2 <- sample(1:nrow(spt), nrow(spt))
            s3 <- sample(1:nrow(spt), nrow(spt))
            s4 <- sample(1:nrow(spt), nrow(spt))
            spt[,3] <- spt[s1,3]
            spt[,4] <- spt[s2,4]
            spt[,5] <- spt[s3,5]
            spt[,6] <- spt[s4,6]
            return(spt)
          })}
        
        gdmPermVar <- lapply(lll, function(x){
          try(gdm(x, geo=geo, splines=splines, knots=knots))
        })
        
        # CHANGE 3 (major): if the result is NULL replace it with a vector of NA values the same length as nPerm
        ##extracts deviance of permuted gdms
        result <-  unlist(sapply(gdmPermVar, function(mod){mod$gdmdeviance}))
        if (is.null(result)) result <- rep(NA, nPerm)
        permVarDev[[k]] <- result
      }
    }
    
    names(permVarDev) <- varNames.x
    nullDev <- fullGDM$nulldeviance
    
    for(var in varNames.x){
      grepper <- grep(paste0("^",var,"$"), names(permVarDev))
      varDevTab <- permVarDev[[grepper]]
      
      # number of perms for which GDM converged
      # CHANGE 4: only count cases where varDevTab is not NA 
      #nConv <- length(varDevTab)
      nConv <- length(na.omit(varDevTab))
      nModsConverge[which(rownames(varImpTable) == var),v] <- nConv
      
      # remove NULLs (GDMs that did not converge)
      #if(nConv>0){
      #varDevTab <- unlist(varDevTab[-(which(sapply(is.null(varDevTab))))])
      #}
      
      # calculate variable importance
      varDevExplained <- 100*(nullDev-varDevTab)/nullDev
      varImpTable[which(rownames(varImpTable) == var),v] <- median(100 * abs((varDevExplained - fullGDM$explained)/fullGDM$explained))
      
      if(var!="Geographic"){
        ##select variable columns to be removed from original site-pair table
        testVarCols1 <- grep(paste("^s1.", var, "$", sep=""), colnames(currSitePair))
        testVarCols2 <- grep(paste("^s2.", var, "$", sep=""), colnames(currSitePair))
        testSitePair <- currSitePair[,-c(testVarCols1, testVarCols2)]
        ##run gdm for the missing variable
        noVarGDM <- gdm(testSitePair, geo=geo, splines=splines[-1], knots=knots)
      } else {
        noVarGDM <- gdm(currSitePair, geo=F, splines=splines[-1], knots=knots)
      }
      
      # calculate p-Value
      # CHANGE 5: if all varDevTab is NA, p-values are NA
      #permDevReduct <- noVarGDM$gdmdeviance - varDevTab	
      #pValues[which(rownames(pValues) == var),v] <- sum(permDevReduct>=(varDevTab - fullGDM$gdmdeviance))/(nConv)
      if (is.null(noVarGDM$gdmdeviance) | all(is.na(varDevTab))){
        pValues[which(rownames(pValues) == var),v] <- NA
      } else {
        permDevReduct <- noVarGDM$gdmdeviance - varDevTab
        pValues[which(rownames(pValues) == var),v] <- sum(permDevReduct>=(varDevTab - fullGDM$gdmdeviance))/(nConv)
      }
    }
    
    if(max(na.omit(pValues[,v]))<pValue){
      message("All remaining predictors are significant, ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", round(fullGDM$explained,3)))
      message("Final set of predictors returned: ")
      for(vvv in 1:length(fullGDM$predictors)){
        message(fullGDM$predictors[vvv])
      }
      break
    }
    
    if(predSelect==T){
      # eliminate least important variable
      #elimVar <- which.min(varImpTable[,v])
      # eliminate variable with highest p-value
      elimVar <- which.max(pValues[,v])
      
      if(names(elimVar)!="Geographic"){
        ##select variable columns to be removed from original site-pair table
        remVarCols1 <- grep(paste("^s1.", names(elimVar), "$", sep=""), colnames(currSitePair))
        remVarCols2 <- grep(paste("^s2.", names(elimVar), "$", sep=""), colnames(currSitePair))
        
        # remove columns from site-pair table
        currSitePair <- currSitePair[,-c(remVarCols1, remVarCols2)]
        # remove columns from permuted site-pair tables
        permSpt <- lapply(permSpt, function(x){
          x[,-c(remVarCols1, remVarCols2)]
        })
      } else {
        geo <- F
      }
      
      varNames.x <- varNames.x[-which(varNames.x==names(elimVar))]
    }
    
    if(v==1 & predSelect==F){
      message("Backwards elimination not selected by user (predSelect=F). Ceasing assessment.")
      message(paste0("Percent deviance explained by final model = ", round(fullGDM$explained,3)))
      message("Final set of predictors returned: ")
      for(vvv in 1:length(fullGDM$predictors)){
        message(fullGDM$predictors[vvv])
      }
      break
    }
  }
  
  
  if(v==1 & predSelect==F){
    # Model assessment
    modelTestVals <- data.frame(matrix(round(modelTestValues[,1], 3), ncol=1))
    rownames(modelTestVals) <- rownames(modelTestValues)
    colnames(modelTestVals) <- "All predictors"
    #Variable importance
    varImpTab <- data.frame(matrix(round(varImpTable[,1], 3), ncol=1))
    rownames(varImpTab) <- rownames(varImpTable)
    colnames(varImpTab) <- "All predictors"
    #Variable selection p-values
    pVals <- varImpTab
    pVals[,1] <- round(pValues[,1], 3)
    #Variable selection model convergence
    nModsConv <- varImpTab
    nModsConv[,1] <- round(nModsConverge[,1], 3)
    outObject <- list(modelTestVals, varImpTab, pVals, nModsConv)
    names(outObject) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
  } else{
    outObject <- list(round(modelTestValues[,1:v], 3), round(varImpTable[,1:v],3), round(pValues[,1:v],3), nModsConverge[,1:v])
    names(outObject) <- c("Model assessment", "Predictor Importance", "Predictor p-values", "Model Convergence")
  }
  
  ##if given, writes out files to space on disk
  if(is.null(outFile)==FALSE){
    save(outObject, file=outFile)
  }
  return(outObject)
}

# GDM function needed by gdm.varImp_custom
permutateSitePair <- function(spTab, siteVarTab, indexTab, vNames){
  
  #################
  #spTab <- currSitePair    ##site-pair table
  #siteVarTab <- siteData   ##siteXvar table
  #indexTab <- indexTab     ##table of the index of sites
  #vNames <- varNames       ##variables names
  #vNames <- c("awcA", "phTotal", "shcA", "solumDepth", "bio5", "bio19")
  #################
  
  ##randomizes the row order of the given siteXvar table
  randVarTab <- siteVarTab[sample(nrow(siteVarTab), nrow(siteVarTab)), ]
  
  #site1x <- siteVarTab$xCoord[1]
  #site1y <- siteVarTab$yCoord[1]
  #checkingIn <- siteVarTab[siteVarTab$xCoord==site1x & siteVarTab$yCoord==site1y,]
  #checkX <- siteVarTab[siteVarTab$xCoord==site1x,]
  #checkingRand <- randVarTab[randVarTab$xCoord==site1x & randVarTab$yCoord==site1y,]
  
  ##sets up the coordinate values for the randomized site-pair table
  s1xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],1]})
  s1yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,1],2]})
  s2xCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],1]})
  s2yCoord <- sapply(1:nrow(spTab), function(i){randVarTab[indexTab[i,2],2]})
  
  #print(vNames)
  ##extracts values of other variables
  varLists <- lapply(vNames, function(vn, rvTab, spt, inT){if(vn!="Geographic"){
    ###################
    #vn <- vNames[[2]]
    #rvTab=randVarTab
    #spt=spTab
    #inT=indexTab
    ###################
    ##identifies variable columns in randVarTab
    randCols <- grep(paste("^", vn, "$", sep=""), colnames(rvTab))
    #print(randCols)
    ##identifies variable columns in site-pair table
    spCols <- grep(vn, colnames(spt))
    
    s1var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,1],randCols]})
    s2var <- sapply(1:nrow(spt), function(i){rvTab[inT[i,2],randCols]})
    
    return(list(s1var, s2var))
  }
  }, rvTab=randVarTab, spt=spTab, inT=indexTab)
  
  # unravels the varList into a data.frame of the variable portion of a site-pair table
  bySite <- lapply(1:2, function(i,vlist){sapply(vlist, function(vl,k){vl[[k]]},k=i)}, vlist=varLists)
  
  if(is(bySite[[1]], "list")){
    site1Vars <- do.call("cbind", bySite[[1]])
    site2Vars <- do.call("cbind", bySite[[2]])
  }else{
    site1Vars <- bySite[[1]]
    site2Vars <- bySite[[2]]
  }
  
  ##sets up new site-pair table
  newSP <- as.data.frame(cbind(spTab$distance, spTab$weights, s1xCoord, s1yCoord, s2xCoord, s2yCoord, site1Vars, site2Vars))
  colnames(newSP) <- colnames(spTab)
  class(newSP) <- c(class(spTab))
  
  #getCoords1 <- newSP[newSP$s1.xCoord==site1x,]
  #getCoords2 <- newSP[newSP$s2.xCoord==site1x,]
  
  return(newSP)
}



############
#   MMRR   #
############

# MMRR FUNCTIONS:
# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)
mmrr <- function(Y,X,nperm=50){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR
unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  x<-scale(x, center=TRUE, scale=TRUE)  # Comment this line out if you wish to perform the analysis without standardizing the distance matrices! 
  return(x)
}


# Run MMRR with one enviornmental matrix of combined distances (not used in final analysis)
run_mmrr <- function(gendist, gsd_df){
  
  ##get env vars and coords
  env_dist <- as.matrix(dist(gsd_df[,c("env1", "env2")], diag = TRUE, upper = TRUE))
  geo_dist <- as.matrix(dist(gsd_df[,c("x", "y")], diag = TRUE, upper = TRUE))

  #Run  MMRR
  mmrr_res <- mmrr(gendist, list(geo = geo_dist, env = env_dist), nperm = 50)
  
  #turn results into dataframe
  results <- mmrr_results_df(mmrr_res)
  
  return(results)
}

#create data frame of results
mmrr_results_df <- function(x, name = NULL){
  df <- 
    data.frame(coeff = x$coefficients, p = x$tpvalue, var = names(x$coefficients)) %>%
    filter(var != "Intercept") %>%
    pivot_wider(names_from = var, values_from = c(coeff, p), names_glue = "{var}_{.value}")
  
  env_cols <- grepl("env", names(df)) & grepl("coeff", names(df))
  geo_cols <- grepl("geo", names(df)) & grepl("coeff", names(df))
  df$ratio <- sum(abs(df[,env_cols]), na.rm = TRUE)/abs(as.numeric(df[,geo_cols]))
  
  if (!is.null(name)) colnames(df) <- paste0(colnames(df), "_", name)
  
  return(df)
}

###############
#   GENERAL   #
###############

# calculate different genetic distance measures
calc_dist <- function(gen, distmeasure = "euc"){
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
    
    #SCALE DISTANCE FROM 0 to 1 if max(distance) >1 (gdm only works for 0<vals<1) (?MODIFY?)
    gendist <- range01(pc_dist)
  } else if(distmeasure == "euc"){
    gendist <- as.matrix(dist(gen, diag = TRUE, upper = TRUE))
  } else {
    print("appropriate gen dist measure not specified")
  }
  
  return(gendist)
}

# calculate statistics for MMRR and GDM
stat_ibdibe <- function(sub, full, sig = 0.05){
  err_cols <- colnames(sub)[grepl("coeff", colnames(full)) | grepl("ratio", colnames(full))]
  err <- map(err_cols, ~err_coeff(full[.x], sub[.x])) %>% bind_cols()
  ae <- abs(err)
  colnames(err) <- paste0(colnames(err), "_", "err")
  colnames(ae) <- paste0(colnames(err), "_", "ae")
  
  p_cols <- colnames(full)[grepl("_p", colnames(full))]
  
  # replace NA pvalue with 1 for calculations because NA values means the coefficient was zero so the significance test  should treat it as a negative
  sub[,p_cols] <- purrr::map_dbl(sub[,p_cols], ~ifelse(is.na(.x), 1, .x))
  
  # True positive rate
  TPR <- ((sub[,p_cols] < sig) & (full[,p_cols] < sig))/(full[,p_cols] < sig)
  colnames(TPR) <- paste0(colnames(TPR), "_", "TPR")
  
  # False discovery rate
  FDR <- ((sub[,p_cols] < sig) & !(full[,p_cols] < sig))/(sub[,p_cols] < sig)
  # replace NA with 0 because NA occurs when denominator is 0 (in which case numerator would also be 0 for this calc)
  FDR <- ifelse(is.na(FDR), 0, FDR)
  colnames(FDR) <- paste0(colnames(FDR), "_", "FDR")
  
  df <- data.frame(err, ae, TPR, FDR)
  return(df)
}
