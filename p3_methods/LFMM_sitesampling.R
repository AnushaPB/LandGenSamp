set.seed(42)

library("here") #paths
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")
source("sitesampling_functions.R")

##########
#  LFMM  #
##########


run_lfmm_full <- function(gen, gsd_df, loci_df){
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #PCA to determine number of latent factors
  pc <- prcomp(gen)
  par(pty="s",mfrow=c(1,1))
  eig <- pc$sdev[1:100]^2
  #estimate number of latent factors using quick.elbow (see general functions for description of how this function works)
  #this is a crude way to determine the number of latent factors that is based on an arbitrary "low" value 
  #(low defaults to 0.08, but this was too high imo so I changed it t0 0.05)
  K <- quick.elbow(eig, low = 0.05, max.pc = 0.9)
  plot(eig, xlab = 'PC', ylab = "Variance explained")
  abline(v = K, col= "red", lty="dashed")
  

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
  FDR1 <- FD/(FD + TP)
  
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
  FDR2 <- FD/(FD + TP)
  
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
  FDR1COMBO <- FD/(FD + TP)
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,2] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci2 %in% loci_trait2)
  TPR2COMBO <- TP/length(loci_trait2)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci2 %in% neutral_loci) + sum(lfmm_loci2 %in% loci_trait1)
  FDR2COMBO <- FD/(FD + TP)
  
  #stats for all loci 
  lfmm_loci <- c(lfmm_loci1, lfmm_loci2)
  #calc True Positive Rate
  TP <- sum(lfmm_loci1 %in% loci_trait1) + sum(lfmm_loci2 %in% loci_trait2)
  TPRCOMBO <- TP/length(adaptive_loci)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci %in% neutral_loci) + sum(lfmm_loci1 %in% loci_trait2) + sum(lfmm_loci2 %in% loci_trait1)
  FDRCOMBO <- FD/(FD + TP)
  
  #PLOT TO CHECK RESULTS (for debugging, remove later)
  
  #par(mfrow=c(1,2))
  #plot(-log10(pvalues[,1]), 
  #     pch = 19, 
  #     cex = .2, 
  #     xlab = "SNP", ylab = "-Log P",
  #     col = "grey",
  #     main = "env1")
  #points(loci_trait1, 
  #       -log10(pvalues[,1])[loci_trait1], 
  #       col = "red", 
  #       cex = 1.5)
  #abline(h = -log10(0.05), col="red", lty=2)
  #
  #plot(-log10(pvalues[,2]), 
  #     pch = 19, 
  #     cex = .2, 
  #     xlab = "SNP", ylab = "-Log P",
  #     col = "grey",
  #     main = "env2")
  #points(loci_trait2, 
  #       -log10(pvalues[,2])[loci_trait2], 
  #       col = "red", 
  #       cex = 1.5)
  #abline(h = -log10(0.05), col="red", lty=2)
  
  return(data.frame(K = K,
                    TPRCOMBO = TPRCOMBO, FDRCOMBO = FDRCOMBO, 
                    TPR1COMBO = TPR1COMBO, FDR1COMBO = FDR1COMBO, 
                    TPR2COMBO = TPR2COMBO, FDR2COMBO = FDR2COMBO,
                    TPR1 = TPR1, FDR1 = FDR1, 
                    TPR2 = TPR2, FDR2 = FDR2,
                    TOTALN = length(lfmm_loci), TOTALT = TP, TOTALF = FD))
}




run_lfmm <- function(gen, gsd_df, loci_df, K){
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #PCA to determine number of latent factors
  pc <- prcomp(gen)
  
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
  FDR1 <-FD/(FD + TP)
  
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
  FDR2 <- FD/(FD + TP)
  
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
  TP1 <- sum(lfmm_loci1 %in% loci_trait1)
  TPR1COMBO <- TP1/length(loci_trait1)
  #calc False Discovery Rate 
  FD1 <- sum(lfmm_loci1 %in% neutral_loci) + sum(lfmm_loci1 %in% loci_trait2)
  FDR1COMBO <- FD1/(FD1 + TP1)
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,2] < 0.05) 
  #calc True Positive Rate
  TP2 <- sum(lfmm_loci2 %in% loci_trait2)
  TPR2COMBO <- TP2/length(loci_trait2)
  #calc False Discovery Rate 
  FD2 <- sum(lfmm_loci2 %in% neutral_loci) + sum(lfmm_loci2 %in% loci_trait1)
  FDR2COMBO <- FD2/(FD2 + TP2)
  
  #stats for all loci 
  lfmm_loci <- c(lfmm_loci1, lfmm_loci2)
  
  #calc True Positive Rate
  TP <- sum(lfmm_loci %in% adaptive_loci)
  TPRCOMBO <- TP/length(adaptive_loci)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci %in% neutral_loci) + sum(lfmm_loci %in% adaptive_loci)
  FDRCOMBO <- FD/(FD + TP)

  
  return(data.frame(K = K,
                    TPRCOMBO = TPRCOMBO, FDRCOMBO = FDRCOMBO, 
                    TPR1COMBO = TPR1COMBO, FDR1COMBO = FDR1COMBO, 
                    TPR2COMBO = TPR2COMBO, FDR2COMBO = FDR2COMBO,
                    TPR1 = TPR1, FDR1 = FDR1, 
                    TPR2 = TPR2, FDR2 = FDR2,
                    TOTALN = length(lfmm_loci), TOTALT = TP, TOTALF = FD))
}


nsites <- c(9, 16, 25)

#register cores
cores <- 10
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

system.time(
res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #create pdf to store plots
  #pdf(paste0("outputs/LFMM/plots/lfmm_plots_",paramset,".pdf"))
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  loci_filepath <- create_filepath(i, params = params, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
                      print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_lfmm_full(gen_2k, gsd_df_2k, loci_df)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = 2000, full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LFMM/LFMM_sitesampling_results_",paramset,".csv")
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
        #sub_result <- run_lfmm(subgen, subgsd_df, loci_df, K = full_result$K)
        sub_result <- run_lfmm(sitegen, sitegsd_df, loci_df, K = full_result$K)
        
        #save and format new result
        sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsite, sub_result)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, sub_result)
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  #end pdf()
  #dev.off()
  
  return(result)
  
  gc()
}
)

#stop cluster
stopCluster(cl)

write.csv(res_lfmm, "outputs/LFMM/lfmm_sitesampling_results.csv", row.names = FALSE)

