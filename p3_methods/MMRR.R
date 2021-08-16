set.seed(42)

library(here) #paths

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")


############
#   MMRR   #
############

# MMRR FUNCTIONS:
# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=999){
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



run_mmrr <- function(gen, gsd_df, npcs = 20){
  #Format data for MMRR  
  ##calculate genetic distance based on pca
  Y <- as.matrix(gen)
  pc <- prcomp(Y)
  pc_dist <- as.matrix(dist(pc$x[,1:npcs], diag = TRUE, upper = TRUE)) #CHANGE NUMBER OF PCS? (see Shirk et al. 2016:  10.1111/1755-0998.12684)
  
  ##get env vars and coords
  env_dist1 <- as.matrix(dist(gsd_df$env1, diag = TRUE, upper = TRUE))
  env_dist2 <- as.matrix(dist(gsd_df$env2, diag = TRUE, upper = TRUE))
  geo_dist <- as.matrix(dist(gsd_df[,c("x", "y")], diag = TRUE, upper = TRUE))
  
  #format X matrices
  Xmats <- list(env1 = env_dist1, env2 = env_dist2, geography = geo_dist)
  
  #Run  MMRR
  mmrr_res <- MMRR(pc_dist, Xmats, nperm = 99)
  
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
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


res_mmrr <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100)
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  skip_to_next <- FALSE
  if(file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run MMRR
  if(skip_to_next == FALSE){
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen <- gen[s,]
    gsd_df <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_mmrr(gen, gsd_df)
    result <- data.frame(sampstrat = "full", nsamp = nrow(gsd_df), full_result, env1_rmse = NA, env2_rmse = NA, geo_rmse = NA)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/MMRR/MMRR_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_mmrr(subgen, subgsd_df)
        
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
  
}


stats_out <- cbind.data.frame(params, res_mmrr)
write.csv(stats_out, "outputs/MMRR/MMRR_results_coeffs.csv")


#stop cluster
stopCluster(cl)