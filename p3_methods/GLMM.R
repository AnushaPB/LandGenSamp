
## GLMM ##
library(spaMM)
library(statgenGWAS)
#parallel
library(foreach)
library(doParallel)
#for reading in files
library(here)
library(vcfR)
#read in general functions and objects
source("general_functions.R")


# create class which holds multiple results for each parallel loop iteration.
# Each loop iteration populates three properties: $result1 and $result2 and $result3
# For a great tutorial on S3 classes, see: 
# http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
multiResultClass <- function(pVal=NULL,coeffs=NULL)
{
  me <- list(
    pVal = pVal,
    coeffs = coeffs
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}


run_glmm <- function(gen_filepath, gsd_filepath, loci_filepath){
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #get adaptive loci
  loci_df <- read.csv(loci_filepath)
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
    
  #make kinship matrix (method = c("astle", "IBS", "vanRaden"))
  Kmatrix <- kinship(gen, method = "astle")
  #Uncomment code below to plot the kinship matrix (with a subsample of 50 individuals)
  #s <- sample(1:nrow(gea_df),50)
  #ks <- Kmatrix[s,s]
  #diag(ks) <- NA
  #heatmap(ks, keep.dendro=FALSE)
  
  
  #WRITTEN AS A PARALLEL LOOP but this won't work if the outerloop is parallelized (in which case it will just run sequentially, so I kept the code the same)
  #register cores
  #cores=detectCores()
  #cl <- makeCluster(cores[1]-2) #not to overload your computer
  #registerDoParallel(cl)
  
  
  modres <- foreach(loci=1:ncol(gen)) %dopar% {
    library(spaMM)
    result <- multiResultClass()
    
    # Data frame for GLMM analysis of each SMV
    df <- data.frame(pop=1:nrow(gen), env1 = gsd_df$env1, env2 = gsd_df$env2, SNP=gen[,loci])
    
    # Fit GLMM using spaMM
    models <- list()
    models[[1]] <- corrHLfit(SNP ~ env1 + corrMatrix(1|pop), data = df, corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
    
    models[[2]] <- corrHLfit(SNP ~ env2 + corrMatrix(1|pop), data = df, corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
    
    models[[3]] <- corrHLfit(SNP ~ env1 + env2 + corrMatrix(1|pop), data = df, corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
    
    null.fit <- corrHLfit(SNP ~ 1 + corrMatrix(1|pop), data = df,
                          corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
    
    LLs <- sapply(models, logLik)
    crit.scores <- sapply(models, AIC)
    cond.AICs <- crit.scores[2,]
    best <- which(cond.AICs == min(cond.AICs))
    result$pVal <- 1-pchisq(2*(LLs[best]-logLik(null.fit)), df=1)
    result$coeffs <- models[[best]]$fixef
    return(result)
  }


  results.tab <- matrix(nrow=ncol(gen), ncol=3)
  for(z in 1:ncol(gen)){
    results.tab[z, 1] <- modres[[z]]$pVal
    results.tab[z, 2] <- modres[[z]]$coeffs["env1"]
    results.tab[z, 3] <- modres[[z]]$coeffs["env2"]
  }
  
  pVal <- c()
  for(i in 1:length(modres)){pVal[i] <- modres[[i]]$pVal}
  
  results.tab <- cbind(results.tab, p.adjust(pVal, method="fdr"))
  colnames(results.tab) = c("p", "env1", "env2", "p.adj")
  results.tab
  
  #env1 candidate loci
  #Identify glmm cand loci (based on pvalue < 0.05 and a coeff that is not NA for the env variable of interest)
  glmm_loci1 <- which(results.tab[,"p.adj"] < 0.05 & !is.na(results.tab[,"env1"]))
  #calc True Positive Rate
  TP <- sum(glmm_loci1 %in% loci_trait1)
  TPR1 <- TP/length(loci_trait1)
  #calc False Discovery Rate 
  FD <- sum(glmm_loci1 %in% neutral_loci) + sum(glmm_loci1 %in% loci_trait2)
  FDR1 <- FD/length(glmm_loci1)
  
  #env2 candidate loci
  #Identify glmm cand loci
  glmm_loci2 <- which(results.tab[,"p.adj"] < 0.05 & !is.na(results.tab[,"env2"]))
  #calc True Positive Rate
  TP <- sum(glmm_loci2 %in% loci_trait2)
  TPR2 <- TP/length(loci_trait2)
  #calc False Discovery Rate 
  FD <- sum(glmm_loci2 %in% neutral_loci) + sum(glmm_loci2 %in% loci_trait1)
  FDR2 <- FD/length(glmm_loci2)
  
  #stats for all loci 
  glmm_loci <- c(glmm_loci1, glmm_loci2)
  #calc True Positive Rate
  TP <- sum(glmm_loci %in% adaptive_loci)
  TPR <- TP/length(adaptive_loci)
  #calc False Discovery Rate 
  FD <- sum(glmm_loci %in% neutral_loci) + sum(glmm_loci %in% adaptive_loci)
  FDR <- FD/length(adaptive_loci)
  
  return(data.frame(TPR = TPR, FDR = FDR, TPR1 = TPR1, FDR1 = FDR1, TPR2 = TPR2, FDR2 = FDR2))
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library(here)
  
  gen_filepath <- create_filepath(i, "gen")
  
  gsd_filepath <- create_filepath(i, "gsd")
  
  loci_filepath <- create_filepath(i, "loci")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(exists(loci_filepath) == FALSE | exists(gen_filepath) == FALSE | exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_lfmm(gen_filepath = gen_filepath,
                       gsd_filepath = gsd_filepath,
                       loci_filepath = loci_filepath)
  }
  
  return(result)
  
}
#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library(here)
  
  gen_filepath <- create_filepath(i, "gen")
  
  gsd_filepath <- create_filepath(i, "gsd")
  
  loci_filepath <- create_filepath(i, "loci")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(exists(loci_filepath) == FALSE | exists(gen_filepath) == FALSE | exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_lfmm(gen_filepath = gen_filepath,
                       gsd_filepath = gsd_filepath,
                       loci_filepath = loci_filepath)
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_glmm)
write.csv(stats_out, "GLMM_results.csv")




