
run_lfmm_helper <- function(x, gen, gsd_df, loci_df){
  return(run_lfmm(gen, gsd_df, loci_df, K = NULL, K_selection = x["K_selection"], method = x["method"]))
}

run_lfmm <- function(gen, gsd_df, loci_df, K = NULL, K_selection = "tracy.widom", method = "ridge"){
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  
  #if K is not specified, it is automatically calculated
  if(is.null(K)) K <- get_K(gen, K_selection = K_selection) 
  
  #gen matrix
  genmat = as.matrix(gen)
  #env matrix
  env1mat = as.matrix(gsd_df$env1)
  env2mat = as.matrix(gsd_df$env2)
  envmat = cbind(env1mat, env2mat)
  
  #BOTH ENV
  #run model
  if (method == "ridge") lfmm_mod <- lfmm_ridge(genmat, envmat, K = K)
  if (method == "lasso") lfmm_mod <- lfmm_lasso(genmat, envmat, K = K)
  
  #performs association testing using the fitted model:
  pv <- lfmm_test(Y = genmat, 
                  X = envmat, 
                  lfmm = lfmm_mod, 
                  calibrate = "gif")
  
  # correct pvals and get confusion matrix stats
  p05 <- purrr::map_dfr(c("none", "fdr", "holm", "bonferroni"), calc_confusion, pv, loci_trait1, loci_trait2, alpha = 0.05)
  p10 <- purrr::map_dfr(c("none", "fdr", "holm", "bonferroni"), calc_confusion, pv, loci_trait1, loci_trait2, alpha = 0.10)
  pdf <- rbind.data.frame(p05, p10)
  df <- data.frame(K = K, K_method = K_selection, lfmm_method = method, pdf)
  
  return(df)
}

calc_confusion <- function(padj, pv, loci_trait1, loci_trait2, alpha = 0.05){
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  # adjust pvalues (or passs through if padj = "none")
  pvalues <-  data.frame(env1 = p.adjust(pv$calibrated.pvalue[,1], method = padj),
                         env2 = p.adjust(pv$calibrated.pvalue[,2], method = padj))
  
  #env1 candidate loci
  #Identify LFMM cand loci (P)
  lfmm_loci1 <- which(pvalues[,"env1"] < alpha) 
  #Identify negatives
  lfmm_neg1 <- which(pvalues[,"env1"] >= alpha | is.na(pvalues[,"env1"]))
  #check length makes sense
  stopifnot(length(lfmm_loci1) + length(lfmm_neg1) == nrow(pvalues))
  
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
  stopifnot(sum(TP1, FP1, TN1, FN1) == nrow(pvalues))
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci2 <- which(pvalues[,"env2"] < alpha) 
  #Identify negatives
  lfmm_neg2 <- which(pvalues[,"env2"] >= alpha | is.na(pvalues[,"env2"]))
  #check length makes sense
  stopifnot(length(lfmm_loci2) + length(lfmm_neg2) == nrow(pvalues))
  
  #True Positives
  TP2 <- sum(lfmm_loci2 %in% loci_trait2)
  #False Positives
  FP2 <- sum(lfmm_loci2 %notin% loci_trait2)
  #True Negatives
  TN2 <- sum(lfmm_neg2 %notin% loci_trait2)
  #False Negatives
  FN2 <- sum(lfmm_neg2 %in% loci_trait2)
  #check length makes sense
  stopifnot(sum(TP2, FP2, TN2, FN2) == nrow(pvalues))
  
  #stats for all loci 
  lfmm_loci <- c(lfmm_loci1, lfmm_loci2)
  #calc confusion matrix
  TP <- TP1 + TP2
  FP <- FP1 + FP2
  TN <- TN1 + TN2
  FN <- FN1 + FN2
  #check sum makes sense
  stopifnot(sum(TP, FP, TN, FN) == 2*nrow(pvalues))
  
  #calc True Positive Rate (i.e. Sensitivity)
  TPRCOMBO <- TP/(TP + FN)
  #calc True Negative Rate (i.e. Specificity)
  TNRCOMBO <- TN/(TN + FP)
  #calc False Discovery Rate 
  FDRCOMBO <- FP/(FP + TP)
  #calc False Positive Rate 
  FPRCOMBO <- FP/(FP + TN)
  
  # Calculate empirical pvalues (I THINK - CHECK THIS)
  null1 <- pvalues$env1[-loci_trait1]
  emp1 <- sapply(pvalues$env1[loci_trait1], function(x){mean(x > null1, na.rm = TRUE)})
  emp1_TPR <- sum(emp1 < alpha, na.rm  = TRUE)
  null2 <- pvalues$env2[-loci_trait2]
  emp2 <- sapply(pvalues$env2[loci_trait2], function(x){mean(x > null2, na.rm = TRUE)})
  emp2_TPR <- sum(emp2 < alpha, na.rm  = TRUE)
  
  return(data.frame(padj = padj,
                    alpha = alpha,
                    TPRCOMBO = TPRCOMBO, 
                    TNRCOMBO = TNRCOMBO,
                    FDRCOMBO = FDRCOMBO, 
                    FPRCOMBO = FPRCOMBO,
                    TOTALN = length(lfmm_loci), 
                    TOTALTP = TP, 
                    TOTALFP = FP, 
                    TOTALTN = TN,
                    TOTALFN = FN,
                    emp1_TPR = emp1_TPR,
                    emp2_TPR = emp2_TPR,
                    emp1_mean = mean(emp1, na.rm = TRUE),
                    emp2_mean = mean(emp2, na.rm = TRUE)))
}


# function to determine K
get_K <- function(gen, coords = NULL, k_selection = "find.clusters", ...){
  
  if(k_selection == "tracy.widom") K <- get_K_tw(gen)
  
  if(k_selection == "quick.elbow") K <- get_K_elbow(gen)
  
  if(k_selection == "find.clusters") K <- get_K_fc(gen)
  
  return(K)
}

# Determine best K using find.clusters
get_K_fc <- function(gen, max.n.clust = NULL, perc.pca = 90){
  if(is.null(max.n.clust)) max.n.clust <- nrow(gen) - 1
  fc <- adegenet::find.clusters(gen,  pca.select = "percVar", perc.pca = perc.pca, choose.n.clust = FALSE, criterion = "diffNgroup", max.n.clust = max.n.clust)
  K <- max(as.numeric(fc$grp))
  return(K)
}

# Function to determine best K based on elbow
get_K_elbow <- function(gen){
  # run pca
  pc <- prcomp(gen)
  
  # get eig
  eig <- pc$sdev^2
  # estimate number of latent factors using quick.elbow (see general functions for description of how this function works)
  # this is a crude way to determine the number of latent factors that is based on an arbitrary "low" value 
  K <- quick.elbow(eig, low = 0.08, max.pc = 0.7)
  
  par(pty = "s",mfrow = c(1,1))
  plot(eig, xlab = 'PC', ylab = "Variance explained")
  abline(v = K, col = "red", lty = "dashed")
  
  return(K)
}

# Function to determine best K using tracy widom test
get_K_tw <- function(gen, maxK = NULL){
  # run pca
  pc <- prcomp(gen)
  
  # get eig
  eig <- pc$sdev^2
  
  # reduce K values if max is provided
  if(!is.null(maxK)){eig <- eig[1:maxK]}
  
  # run tracy widom test
  # NOTE: 	
  # the critical point is a numeric value corresponding to the significance level. 
  # If the significance level is 0.05, 0.01, 0.005, or 0.001, 
  # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly. 
  # The default is 2.0234.
  
  tryCatch(tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = 0.9793), 
           error = function(e) {
             err <<- conditionMessage(e)
             write.table(err, "error_msg.txt")
             write.csv(eig, "eig_error.csv")
             
             print(err)
             print(eig)
             print(dim(gen))
             
             stop(err)})
  
  
  # get K based on number of significant eigenvalues
  K <- tw_result$SigntEigenL
  
  # if none are significant, return 1
  if(K == 0) K <- 1
  
  plot(eig)
  abline(v = K)
  
  return(K)
}
