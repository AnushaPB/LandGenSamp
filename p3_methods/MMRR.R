library(here) #paths

#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")



############
#   TEST   #
############
#define nloci 
nloci = 10000

#read in geospatial data
file_path = here("data","mod-10k_K2_phi50_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.csv")
gsd_df <- read.csv(file_path)
#Extract env values from lists
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
#Extract phenotype values from lists
gsd_df$z1 <- as.numeric(stringr::str_extract(gsd_df$z, '(?<=\\[)[^,]+(?=,)'))
gsd_df$z2 <- as.numeric(stringr::str_extract(gsd_df$z, '(?<=, )[^,]+(?=\\])')) 
head(gsd_df)

#read in genetic data
file_path = here("data","mod-10k_K2_phi50_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.vcf")
vcf <- read.vcfR(file_path)
#convert to matrix
x <- vcfR2genlight(vcf) 
gen <- as.matrix(x)

#FOR TESTING USE SUBSAMPLE(!!!COMMENT OUT!!!)
set.seed(42)
s <- sample(1:nrow(gea_df),200)
gsd_df <- gsd_df[s,]
gen <- gen[s,]

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



run_mmrr <- function(gen_filepath, gsd_filepath){
  
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #Format data for MMRR  
  ##calculate genetic distance based on pca
  Y <- as.matrix(gen)
  pc <- prcomp(Y)
  pc_dist <- as.matrix(dist(pc$x[,1:100], diag = TRUE, upper = TRUE)) #CHANGE NUMBER OF PCS? (see Shirk et al. 2016:  10.1111/1755-0998.12684)
  
  ##get env vars and coords
  env_dist1 <- as.matrix(dist(gsd_df$env1, diag = TRUE, upper = TRUE))
  env_dist2 <- as.matrix(dist(gsd_df$env2, diag = TRUE, upper = TRUE))
  geo_dist <- as.matrix(dist(gsd_df[,c("x", "y")], diag = TRUE, upper = TRUE))
  
  #format X matrices
  Xmats <- list(env1 = env_dist1, env2 = env_dist2, geography = geo_dist)
  
  #Run  MMRR
  mmrr_res <- MMRR(pc_dist, Xmats, nperm = 999)
  
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
  library(here)
  
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(exists(loci_filepath) == FALSE | exists(gen_filepath) == FALSE | exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_mmrr(gen_filepath = gen_filepath,
                       gsd_filepath = gsd_filepath)
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_mmrr)
write.csv(stats_out, "MMRR_results.csv")
