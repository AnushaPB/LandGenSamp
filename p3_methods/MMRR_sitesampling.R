set.seed(42)

library(here)

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")

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



run_mmrr <- function(gen, gsd_df, distmeasure= "euc"){
  #Format data for MMRR  
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


#register cores
cores <- 5
cl <- makeCluster(cores)
registerDoParallel(cl)


res_mmrr <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("adegenet")
  library("stringr")

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
  
  #run MMRR
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_mmrr(gen_2k, gsd_df_2k)
    fullratio <- (full_result$env1_coeff + full_result$env2_coeff)/full_result$geo_coeff
    result <- data.frame(params[i,], sampstrat = "full", nsamp = nrow(gsd_df_2k), full_result, 
                         ratio = fullratio, 
                         env1_err = NA, env2_err = NA, geo_err = NA, ratio_err = NA)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/MMRR/MMRR_sitesampling_results_",paramset,".csv")
    write.csv(result, csv_file, row.names = FALSE)
    
    for(nsite in nsites){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsite)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #get sites
        siteIDs <- get_sites(params[i,], params, sampstrat, nsite)
        #confirm that number of sites matches number of sample IDs
        stopifnot(length(subIDs) == length(siteIDs))
        #calculate allele frequency by site (average)
        sitegen <- data.frame(aggregate(subgen, list(siteIDs), FUN=mean)[,-1])
        #calculate env values by site
        sitegsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN=mean)[,-1]) 
        
        #run analysis using subsample
        sub_result <- run_mmrr(subgen, subgsd_df)
        
        #calculate err
        env1_err <- err_coeff(full_result$env1_coeff, sub_result$env1_coeff)
        env2_err <- err_coeff(full_result$env2_coeff, sub_result$env2_coeff)
        geo_err <- err_coeff(full_result$geo_coeff, sub_result$geo_coeff)
        
        subratio <- (sub_result$env1_coeff + sub_result$env2_coeff)/sub_result$geo_coeff
        ratio_err <- err_coeff(fullratio, subratio)
        
        #save and format new result
        sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsite, sub_result, 
                                 ratio = subratio,
                                 env1_err = env1_err, env2_err = env2_err, geo_err = geo_err,
                                 ratio_err = ratio_err)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, sub_result)
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
        
      }
    }
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

write.csv(res_mmrr, "outputs/mmrr_sitesampling_results.csv", row.names = FALSE)
