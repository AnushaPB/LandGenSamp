
run_lfmm <- function(gen, gsd_df, loci_df, K = NULL){
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  
  #if K is not specified it is calculated based on a tracy widom test
  if(is.null(K)){
    K <- get_K_tw(gen)
  }
  
  
  #gen matrix
  genmat = as.matrix(gen)
  #env matrix
  env1mat = as.matrix(gsd_df$env1)
  env2mat = as.matrix(gsd_df$env2)
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
  lfmm_loci1 <- which(pvalues[,"env1"] < 0.05) 
  #Identify negatives
  lfmm_neg1 <- which(pvalues[,"env1"] >= 0.05 | is.na(pvalues[,"env1"]))
  #check length makes sense
  stopifnot(length(lfmm_loci1) + length(lfmm_neg1) == ncol(gen))
  
  #get confusion matrix values
  #True Positives
  TP1 <- sum(lfmm_loci1 %in% loci_trait1)
  #False Positives
  FP1 <- sum(lfmm_loci1 %notin% loci_trait1)
  #True Negatives
  TN1 <- sum(lfmm_neg1 %notin% loci_trait1)
  #False Negatives
  FN1 <- sum(lfmm_neg1 %in% loci_trait1)
  #check sum makes sense
  stopifnot(sum(TP1, FP1, TN1, FN1) == ncol(gen))
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,"env2"] < 0.05) 
  #Identify negatives
  lfmm_neg2 <- which(pvalues[,"env2"] >= 0.05 | is.na(pvalues[,"env2"]))
  #check length makes sense
  stopifnot(length(lfmm_loci2) + length(lfmm_neg2) == ncol(gen))
  
  #True Positives
  TP2 <- sum(lfmm_loci2 %in% loci_trait2)
  #False Positives
  FP2 <- sum(lfmm_loci2 %notin% loci_trait2)
  #True Negatives
  TN2 <- sum(lfmm_neg2 %notin% loci_trait2)
  #False Negatives
  FN2 <- sum(lfmm_neg2 %in% loci_trait2)
  #check length makes sense
  stopifnot(sum(TP2, FP2, TN2, FN2) == ncol(gen))
  
  #stats for all loci 
  lfmm_loci <- c(lfmm_loci1, lfmm_loci2)
  #calc confusion matrix
  TP <- TP1 + TP2
  FP <- FP1 + FP2
  TN <- TN1 + TN2
  FN <- FN1 + FN2
  #check sum makes sense
  stopifnot(sum(TP, FP, TN, FN) == 2*ncol(gen))
  
  #calc True Positive Rate (i.e. Sensitivity)
  TPRCOMBO <- TP/(TP + FN)
  #calc True Negative Rate (i.e. Specificity)
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

get_K_tw <- function(gen){
  # run pca
  pc <- prcomp(gen)
  
  # get eig
  eig <- pc$sdev^2
  
  # run tracy widom test
  # NOTE: 	
  # the critical point is a numeric value corresponding to the significance level. 
  # If the significance level is 0.05, 0.01, 0.005, or 0.001, 
  # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. 
  # The default is 2.0234.
  tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = 0.9793)
  
  # get K based on number of significant eigenvalues
  K <- tw_result$SigntEigenL
  
  plot(eig)
  
  return(K)
}

run_sub <- function(sampcombo, i, params, gen, gsd_df, mode = "ind"){
  
  # get nsamp and sampstrat
  nsamp <- sampcombos[1]
  sampstrat <- sampcombos[2]
  
  #subsample from data based on sampling strategy and number of samples
  subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
  subgen <- gen[subIDs,]
  subgsd_df <- gsd_df[subIDs,]
  
  if(mode == "site"){
    #get sites
    siteIDs <- get_sites(params[i,], params, sampstrat, nsite)
    #confirm that number of sites matches number of sample IDs
    stopifnot(length(subIDs) == length(siteIDs))
    #calculate allele frequency by site (average)
    subgen <- data.frame(aggregate(subgen, list(siteIDs), FUN=mean)[,-1])
    #calculate env values by site
    subgsd_df <- data.frame(aggregate(subgsd_df, list(siteIDs), FUN=mean)[,-1]) 
  }
  
  #run analysis using subsample
  sub_result <- run_lfmm(subgen, subgsd_df, loci_df, K = NULL)
  
  #save and format new result
  sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsamp, sub_result)
  
  #bind results
  return(sub_result)
}
