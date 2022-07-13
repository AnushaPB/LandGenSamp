############
#   MMRR   #
############

# MMRR FUNCTIONS:
# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=50){
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


# Run MMRR
run_mmrr <- function(gen, gsd_df, distmeasure= "euc"){
  #Format data for MMRR  
  gendist <- calc_dist(gen, distmeasure)
  
  ##get env vars and coords
  env_dist1 <- as.matrix(dist(gsd_df$env1, diag = TRUE, upper = TRUE))
  env_dist2 <- as.matrix(dist(gsd_df$env2, diag = TRUE, upper = TRUE))
  geo_dist <- as.matrix(dist(gsd_df[,c("x", "y")], diag = TRUE, upper = TRUE))
  
  #format X matrices
  Xmats <- list(env1 = env_dist1, env2 = env_dist2, geography = geo_dist)
  
  #Run  MMRR
  mmrr_res <- MMRR(gendist, Xmats, nperm = 50)
  
  #create data frame of results
  mmrr_df <- cbind.data.frame(mmrr_res$coefficients, mmrr_res$tpvalue)
  
  #turn results into dataframe
  results <- data.frame(env1_coeff = mmrr_res$coefficients["env1"],
                        env2_coeff = mmrr_res$coefficients["env2"],
                        geo_coeff = mmrr_res$coefficients["geography"],
                        env1_p = mmrr_res$tpvalue["env1(p)"],
                        env2_p = mmrr_res$tpvalue["env2(p)"],
                        geo_p = mmrr_res$tpvalue["geography(p)"]
  )
  #remove rownames
  rownames(results) <- NULL
  
  return(results)
}

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
