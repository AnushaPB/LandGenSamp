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


run_lfmm_trait <- function(gen, gsd_df, loci_df, K = NULL){
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #PCA to determine number of latent factors
  #if K is not specified it is calculated based on PCs
  if(is.null(K)){
    pc <- prcomp(gen)
    par(pty="s",mfrow=c(1,1))
    eig <- pc$sdev[1:100]^2
    #estimate number of latent factors using quick.elbow (see general functions for description of how this function works)
    #this is a crude way to determine the number of latent factors that is based on an arbitrary "low" value 
    #(low defaults to 0.08, but this was too high imo so I changed it t0 0.05)
    K <- quick.elbow(eig, low = 0.05, max.pc = 0.9)
    plot(eig, xlab = 'PC', ylab = "Variance explained")
    abline(v = K, col= "red", lty="dashed")
  }
  
  
  #gen matrix
  genmat = as.matrix(gen)
  #env matrix
  env1mat = as.matrix(gsd_df$z1)
  env2mat = as.matrix(gsd_df$z2)
  envmat = cbind(env1mat, env2mat)
  
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
  #Identify LFMM cand loci (P)
  lfmm_loci1 <- which(pvalues[,1] < 0.05) 
  #Identify negatives
  lfmm_neg1 <- which(!(pvalues[,1] < 0.05))
  #get confusion matrix values
  #True Positives
  TP1 <- sum(lfmm_loci1 %in% loci_trait1)
  #False Positives
  FP1 <- sum(lfmm_loci1 %notin% loci_trait1)
  #True Negatives
  TN1 <- sum(lfmm_neg1 %notin% loci_trait1)
  #False Negatives
  FN1 <- sum(lfmm_neg1 %in% loci_trait1)
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,2] < 0.05) 
  #Identify negatives
  lfmm_neg2 <- which(!(pvalues[,2] < 0.05))
  #True Positives
  TP2 <- sum(lfmm_loci2 %in% loci_trait2)
  #False Positives
  FP2 <- sum(lfmm_loci2 %notin% loci_trait2)
  #True Negatives
  TN2 <- sum(lfmm_neg2 %notin% loci_trait2)
  #False Negatives
  FN2 <- sum(lfmm_neg2 %in% loci_trait2)
  
  #stats for all loci 
  lfmm_loci <- c(lfmm_loci1, lfmm_loci2)
  #calc confusion matrix
  TP <- TP1 + TP2
  FP <- FP1 + FP2
  TN <- TN1 + TN2
  FN <- FN1 + FN2
  
  #calc True Positive Rate
  TPRCOMBO <- TP/(TP + FN)
  #calc True Negative Rate
  TNRCOMBO <- TN/(TN + FP)
  #calc False Discovery Rate 
  FDRCOMBO <- FP/(FP + TP)
  #calc False Positive Rate
  FPRCOMBO <- FP/(FP + TN)
  
  return(data.frame(K = K,
                    TPRCOMBO = TPRCOMBO, 
                    TNRCOMBO = TNRCOMBO,
                    FDRCOMBO = FDRCOMBO, 
                    FPRCOMBO = FPRCOMBO,
                    TOTALN = length(lfmm_loci), 
                    TOTALTP = TP, 
                    TOTALFP = FP, 
                    TOTALTN = TN,
                    TOTALFN = FN))
}

#register cores
cores <- 10
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

system.time(
res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  library("stringr")
  
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
    full_result <- run_lfmm_trait(gen_2k, gsd_df_2k, loci_df)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = 2000, full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LFMM_trait/LFMM_trait_sitesampling_results_",paramset,".csv")
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
        #calculate env/z values by site
        sitegsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN=mean)[,-1]) 
        
        #run analysis using subsample
        #sub_result <- run_lfmm(subgen, subgsd_df, loci_df, K = full_result$K)
        sub_result <- run_lfmm_trait(sitegen, sitegsd_df, loci_df, K = full_result$K)
        
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

write.csv(res_lfmm, "outputs/LFMM_trait/lfmm_trait_sitesampling_results.csv", row.names = FALSE)

