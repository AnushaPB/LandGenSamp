library("here") #paths
library("vegan") #RDA
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")

##########
#  LFMM  #
##########

run_lfmm <- function(gen_filepath, gsd_filepath, loci_filepath){
  
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #get adaptive loci
  loci_df <- read.csv(loci_filepath)
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #PCA to determine number of latent factors
  pc <- prcomp(gen)
  par(pty="s",mfrow=c(1,1))
  plot(pc$sdev[1:100]^2, xlab = 'PC', ylab = "Variance explained")
  K <- 10 #NUMBER OF LATENT FACTORS (NEED TO MODIFY TO MAKE AUTO)
  
  #gen matrix
  genmat = as.matrix(gen)
  #env matrix
  env1mat = as.matrix(gsd_df$env1)
  env2mat = as.matrix(gsd_df$env2)
  envmat = cbind(env1mat, env2mat)
  
  #ENV1
  #run model
  lfmm_mod <- lfmm_ridge(genmat, env1mat, K = K)
  #performs association testing using the fitted model:
  pv <- lfmm_test(Y = genmat, 
                  X = env1mat, 
                  lfmm = lfmm_mod, 
                  calibrate = "gif")
  #adjust pvalues
  pvalues <- data.frame(env1=p.adjust(pv$calibrated.pvalue[,1], method="fdr"))
  #env1 candidate loci
  #Identify LFMM cand loci
  lfmm_loci1 <- which(pvalues[,1] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci1 %in% loci_trait1)
  TPR1 <- TP/length(loci_trait1)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci1 %in% neutral_loci) + sum(lfmm_loci1 %in% loci_trait2)
  FDR1 <- FD/length(lfmm_loci1)
  
  #ENV2
  #run model
  lfmm_mod <- lfmm_ridge(genmat, env2mat, K = K)
  #performs association testing using the fitted model:
  pv <- lfmm_test(Y = genmat, 
                  X = env2mat, 
                  lfmm = lfmm_mod, 
                  calibrate = "gif")
  #adjust pvalues
  pvalues <- data.frame(env2=p.adjust(pv$calibrated.pvalue[,1], method="fdr"))
  #env1 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,1] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci2 %in% loci_trait2)
  TPR2 <- TP/length(loci_trait2)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci2 %in% neutral_loci) + sum(lfmm_loci2 %in% loci_trait1)
  FDR2 <- FD/length(lfmm_loci2)
  
  #BOTH ENV
  #run model
  lfmm_mod <- lfmm_ridge(genmat, envmat, K = K)
  #performs association testing using the fitted model:
  pv <- lfmm_test(Y = genmat, 
                  X = envmat, 
                  lfmm = lfmm_mod, 
                  calibrate = "gif")
  #adjust pvalues
  pvalues <- data.frame(env1=p.adjust(pv$calibrated.pvalue[,1], method="fdr"),
                        env2=p.adjust(pv$calibrated.pvalue[,2], method="fdr"))
  #env1 candidate loci
  #Identify LFMM cand loci
  lfmm_loci1 <- which(pvalues[,1] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci1 %in% loci_trait1)
  TPR1COMBO <- TP/length(loci_trait1)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci1 %in% neutral_loci) + sum(lfmm_loci1 %in% loci_trait2)
  FDR1COMBO <- FD/length(lfmm_loci1)
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,2] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci2 %in% loci_trait2)
  TPR2COMBO <- TP/length(loci_trait2)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci2 %in% neutral_loci) + sum(lfmm_loci2 %in% loci_trait1)
  FDR2COMBO <- FD/length(lfmm_loci2)
  
  #stats for all loci 
  lfmm_loci <- c(lfmm_loci1, lfmm_loci2)
  #calc True Positive Rate
  TP <- sum(lfmm_loci %in% adaptive_loci)
  TPRCOMBO <- TP/length(adaptive_loci)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci %in% neutral_loci) + sum(lfmm_loci %in% adaptive_loci)
  FDRCOMBO <- FD/length(adaptive_loci)

  #PLOT TO CHECK RESULTS
  
  par(mfrow=c(1,2))
  plot(-log10(pvalues[,1]), 
       pch = 19, 
       cex = .2, 
       xlab = "SNP", ylab = "-Log P",
       col = "grey",
       main = "env1")
  points(loci_trait1, 
         -log10(pvalues[,1])[loci_trait1], 
         col = "red", 
         cex = 1.5)
  abline(h = -log10(0.05), col="red", lty=2)
  
  plot(-log10(pvalues[,2]), 
       pch = 19, 
       cex = .2, 
       xlab = "SNP", ylab = "-Log P",
       col = "grey",
       main = "env2")
  points(loci_trait2, 
         -log10(pvalues[,2])[loci_trait2], 
         col = "red", 
         cex = 1.5)
  abline(h = -log10(0.05), col="red", lty=2)
  
  return(data.frame(TPRCOMBO = TPRCOMBO, FDRCOMBO = FDRCOMBO, 
                    TPR1COMBO = TPR1COMBO, FDR1COMBO = FDR1COMBO, 
                    TPR2 = TPR2COMBO, FDR2 = FDR2COMBO,
                    TPR1 = TPR1, FDR1 = FDR1, 
                    TPR2 = TPR2, FDR2 = FDR2))
}



#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  
  gen_filepath <- create_filepath(i, "gen")
  
  gsd_filepath <- create_filepath(i, "gsd")
  
  loci_filepath <- create_filepath(i, "loci")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
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
  
  gc()
}
  
#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)
write.csv(stats_out, "LFMM_results.csv")

