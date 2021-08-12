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

##########
#  LFMM  #
##########
# quickly choose an elbow for a PC. 
# at variance below 5% per component, choose the largest % drop
# designed for variance percentages, but will also work given a full set of Evalues
#' Quickly estimate the 'elbow' of a scree plot (PCA)
#' 
#' This function uses a rough algorithm to estimate a sensible 'elbow' to
#' choose for a PCA scree plot of eigenvalues. The function looks at an initial arbitrarily 'low'
#' level of variance and looks for the first eigenvalue lower than this. If the very first eigenvalue
#' is actually lower than this (i.e, when the PCs are not very explanatory) then this 'low' value is
#' iteratively halved until this is no longer the case. After starting below this arbitrary threshold
#' the drop in variance explained by each pair of consecutive PCs is standardized by dividing over the 
#' larger of the pair. The largest percentage drop in the series below 'low' % is selected as the 'elbow'.
#' @param varpc numeric, vector of eigenvalues, or 'percentage of variance' explained datapoints for
#'  each principle component. If only using a partial set of components, should first pass to 
#'  estimate.eig.vpcs() to estimate any missing eigenvalues.
#' @param low numeric, between zero and one, the threshold to define that a principle component
#'  does not explain much 'of the variance'.
#' @param max.pc maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n')
#' @return The number of last principle component to keep, prior to the determined elbow cutoff
#' @export
#' @seealso \code{\link{estimate.eig.vpcs}}
#' @author Nicholas Cooper 
#' @examples
#' # correlated data
#' mat <- sim.cor(100,50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick.elbow(eig)
#' pca.scree.plot(eig,elbow=elb.a,M=mat) 
#' elb.b <- quick.elbow(eig,low=.05) # decrease 'low' to select more components
#' pca.scree.plot(eig,elbow=elb.b,M=mat) 
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5,3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick.elbow(result2$sdev^2)
#' pca.scree.plot(eig2,elbow=elb2,M=mat2)
quick.elbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones 
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else { 
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee) 
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}


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
  
  #PLOT TO CHECK RESULTS (for debugging, remove later)
  
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
  
  return(data.frame(K = K,
                    TPRCOMBO = TPRCOMBO, FDRCOMBO = FDRCOMBO, 
                    TPR1COMBO = TPR1COMBO, FDR1COMBO = FDR1COMBO, 
                    TPR2 = TPR2COMBO, FDR2 = FDR2COMBO,
                    TPR1 = TPR1, FDR1 = FDR1, 
                    TPR2 = TPR2, FDR2 = FDR2))
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
  
  return(data.frame(K = K,
                    TPRCOMBO = TPRCOMBO, FDRCOMBO = FDRCOMBO, 
                    TPR1COMBO = TPR1COMBO, FDR1COMBO = FDR1COMBO, 
                    TPR2 = TPR2COMBO, FDR2 = FDR2COMBO,
                    TPR1 = TPR1, FDR1 = FDR1, 
                    TPR2 = TPR2, FDR2 = FDR2))
}



#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
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
                     "_r",params[i,"r"]*100)
  
  #create pdf to store plots
  pdf(paste0("lfmm_plots_",paramset,".pdf"))
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  loci_filepath <- create_filepath(i, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
                      print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    loci_df <- get_data(i, "loci")
    
    #subsample full data randomly
    s <- sample(2000, nrow(gsd_df), replace = FALSE)
    gen <- gen[s,]
    gsd_df <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_lfmm_full(gen, gsd_df, loci_df)
    result <- data.frame(sampstrat = "full", nsamp = nrow(gsd_df), full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("LFMM_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_lfmm(subgen, subgsd_df, loci_df, K = full_result$K)
        
        #save and format new result
        sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, sub_result)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  #end pdf()
  dev.off()
  
  return(result)
  
  gc()
}
)

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)
write.csv(stats_out, "outputs/LFMM_results.csv")

