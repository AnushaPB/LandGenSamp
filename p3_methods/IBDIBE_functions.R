
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
  
  # model combo
  distPreds <- cbind(site, as.matrix(dist(gsd_df[,c("env1", "env2")], method = "euclidean", diag = TRUE, upper = TRUE)))
  gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, REMOVE = rep(1, nrow(gsd_df)))
   
  gdmData <-
    formatsitepair(
      gdmGen,
      bioFormat = 3,
      predData = gdmPred,
      XColumn = "Longitude",
      YColumn = "Latitude",
      siteColumn = "site",
      distPreds = list(env = distPreds)
    )
  
  #remove placeholder column
  gdmData <- gdmData[,!grepl("*REMOVE*", colnames(gdmData))]
  
  #scale distance from 01
  gdmData$distance <- range01(gdmData$distance) 
  
  #run GDM
  gdm.model <- gdm(gdmData, geo = TRUE)
  
  if(is.null(gdm.model)){
    #turn results into dataframe
    results <- data.frame(env_coeff = NA,
                          geo_coeff = NA,
                          ratio = NA,
                          env_p = NA,
                          geo_p = NA)
  } else {
    predictors <- coeffs(gdm.model)
    
    # turn results into dataframe
    results <- data.frame(env_coeff = predictors[predictors$predictor == "matrix_1", "coefficient"],
                          geo_coeff = predictors[predictors$predictor == "Geographic", "coefficient"])
    results$ratio <- abs(results$env_coeff)/abs(results$geo_coeff)
    
    # get pvalues
    safe_gdm.varImp <- safely(gdm.varImp)
    modTest <- safe_gdm.varImp(gdmData, geo = T, nPerm = 50, parallel = F, predSelect = F)
    if (is.null(modTest$error)) {
      pvals <- modTest$result$`Predictor p-values`
      pvals$var <- row.names(pvals)
      pvals <- left_join(data.frame(var = c("Geographic", "matrix_1")), pvals)
      results <- data.frame(results,
                 env_p = pvals[pvals$var == "matrix_1", 2],
                 geo_p = pvals[pvals$var == "Geographic", 2])
    } else {
      results <- data.frame(results,
                            env_p = NA,
                            geo_p = NA)
    }
    
  }
  
  
  #remove rownames
  rownames(results) <- NULL
  
  return(results)
}

############
#   MMRR   #
############

# MMRR FUNCTIONS:
# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

# NOTE: adjusted for sims
MMRR<-function(Y,X,nperm=50){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X))) names(X) <- paste("X",1:length(X),sep="")
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
  
  stat_df <- data.frame(coeffs, vars = names(coeffs))
  stat_df <- merge(stat_df, data.frame(tstat, vars = names(tstat)), all = TRUE)
  stat_df <- merge(stat_df, data.frame(tp, vars = names(tstat)), all = TRUE)
  stat_df$vars <- c("Intercept",names(X))
  tstat <- stat_df$tstat
  tp <- stat_df$tp
  coeffs <- stat_df$coeffs 
  names(coeffs) <- names(tstat) <- names(tp) <- stat_df$vars
  
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


# Run MMRR
run_mmrr <- function(gen, gsd_df, distmeasure= "euc"){
  #Format data for MMRR  
  gendist <- calc_dist(gen, distmeasure)
  
  ##get env vars and coords
  env_dist <- as.matrix(dist(gsd_df[,c("env1", "env2")], diag = TRUE, upper = TRUE))
  geo_dist <- as.matrix(dist(gsd_df[,c("x", "y")], diag = TRUE, upper = TRUE))

  #Run  MMRR
  mmrr_res <- MMRR(gendist, list(geo = geo_dist, env = env_dist), nperm = 50)
  
  #turn results into dataframe
  results <- mmrr_results_df(mmrr_res)
  
  return(results)
}


mmrr_results_df <- function(x, name = NULL){
  #create data frame of results
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


stat_ibdibe <- function(sub, full, sig = 0.05){
  err_cols <- colnames(sub)[grepl("coeff", colnames(full)) | grepl("ratio", colnames(full))]
  err <- map(err_cols, ~err_coeff(full[.x], sub[.x])) %>% bind_cols()
  ae <- abs(err)
  colnames(err) <- paste0(colnames(err), "_", "err")
  colnames(ae) <- paste0(colnames(err), "_", "ae")
  
  p_cols <- colnames(full)[grepl("_p", colnames(full))]
  TPR <- ((sub[,p_cols] < sig) & (full[,p_cols] < sig))/(full[,p_cols] < sig)
  colnames(TPR) <- paste0(colnames(TPR), "_", "TPR")
  
  df <- data.frame(err, ae, TPR)
  return(df)
}
