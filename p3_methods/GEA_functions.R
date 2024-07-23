############
#   LFMM   #
############

# run LFMM2
run_lfmm <- function(gen, gsd_df, K_selection = c("tess"), lfmm_method = c("ridge", "lasso"), K = NULL, loci_df = NULL, Kvals = 1:9, maf = c(0, 0.05)){
  
  # get adaptive loci
  if (is.null(loci_df)) loci_df <- get_loci()
  
  # For full K: if K is provided, run with that K
  if (!is.null(K)) return(run_lfmm_helper(gen = gen, gsd_df = gsd_df, loci_df = loci_df, K = K, maf = 0.05))
  
  # if no K is provided, run tests to get_K
  result <-
    expand.grid(K_selection = K_selection, lfmm_method = lfmm_method, maf = maf) %>%
    pmap(\(K_selection, lfmm_method, maf) run_lfmm_helper(K_selection = K_selection, lfmm_method = lfmm_method,
                                                gen = gen, gsd_df = gsd_df, loci_df = loci_df, Kvals = Kvals, maf = maf)) %>%
    bind_rows()
  
  return(result)
}

# helper for run_lfmm
run_lfmm_helper <- function(gen, gsd_df, loci_df, K = NULL, K_selection = "tess", lfmm_method = "ridge", Kvals = 1:9, maf = 0){
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 
  loci_trait2 <- loci_df$trait2 
  
  # maf filter
  if (maf > 0) {
    # get colnames of loci
    loci_trait1 <- colnames(gen)[loci_trait1]
    loci_trait2 <- colnames(gen)[loci_trait2]
    
    # perform maf filtering
    gen_maf <- map_dbl(1:ncol(gen), ~mean(gen[,.x], na.rm = TRUE))/2
    gen <- gen[,gen_maf > maf & gen_maf < (1 - maf)]
    
    # reassign loci names with col numbers
    loci_trait1 <- which(colnames(gen) %in% loci_trait1)
    loci_trait2 <- which(colnames(gen) %in% loci_trait2)
  }
  
  #if K is not specified, it is automatically calculated
  if (is.null(K)) K <- get_K(gen, coords = gsd_df[,c("x", "y")], K_selection = K_selection, Kvals = Kvals) 
  
  #gen matrix
  genmat = as.matrix(gen)
  
  #env matrix
  env1mat = as.matrix(gsd_df$env1)
  env2mat = as.matrix(gsd_df$env2)
  envmat = cbind(env1mat, env2mat)
  
  #BOTH ENV
  #run model
  if (lfmm_method == "ridge") lfmm_mod <- tryCatch(lfmm::lfmm_ridge(genmat, envmat, K = K), error = function(x) NULL)
  if (lfmm_method == "lasso") lfmm_mod <- tryCatch(lfmm::lfmm_lasso(genmat, envmat, K = K), error = function(x) NULL)
  if (is.null(lfmm_mod)) return(data.frame(K_method = K_selection, lfmm_method = lfmm_method, NULL_mod = TRUE))
  
  # correct pvals and get confusion matrix stats
  # change padjustment, alpha level, minimum minor allele frequency, and whether to test for all adaptive loci (i.e. combined trait1/2)
  combos <-
    expand.grid(
      padj = c("none", "fdr", "holm", "bonferroni"),
      sig = c(0.05, 0.10),
      all = c(TRUE, FALSE)
    )
  
  result <-
    pmap(
      combos,
      \(padj, sig, maf, all) return(lfmm_calc_confusion(
        padj = as.character(padj),
        genmat = genmat,
        envmat = envmat,
        lfmm_mod = lfmm_mod,
        loci_trait1 = loci_trait1,
        loci_trait2 = loci_trait2,
        sig = sig,
        all = all
      ))
    )
  
  df <- data.frame(K_factor = K, K_method = K_selection, lfmm_method = lfmm_method, NULL_mod = FALSE, maf = maf, bind_rows(result))
  
  return(df)
}



# determine K
get_K <- function(gen, coords = NULL, K_selection = "find.clusters", Kvals = 1:9, ...){
  
  if(K_selection == "tracy.widom") K <- get_K_tw(gen)
  
  if(K_selection == "quick.elbow") K <- get_K_elbow(gen)
  
  if(K_selection == "find.clusters") K <- get_K_fc(gen)
  
  if(K_selection == "tess") K <- get_K_tess(gen, coords, Kvals = Kvals)
  
  return(K)
}

# Determine best K using find.clusters
get_K_fc <- function(gen, max.n.clust = 9, perc.pca = 70){
  if(is.null(max.n.clust)) max.n.clust <- nrow(gen) - 1
  if((nrow(gen) - 1) < max.n.clust) max.n.clust <- nrow(gen) - 1
  fc <- adegenet::find.clusters(gen,  pca.select = "percVar", perc.pca = perc.pca, choose.n.clust = FALSE, criterion = "diffNgroup", max.n.clust = max.n.clust)
  K <- max(as.numeric(fc$grp))
  return(K)
}

# determine best K based on elbow
get_K_elbow <- function(gen){
  # run pca
  pc <- prcomp(gen)
  
  # get eig
  eig <- pc$sdev^2
  # estimate number of latent factors using quick.elbow (see general functions for description of how this function works)
  # this is a crude way to determine the number of latent factors that is based on an arbitrary "low" value 
  K <- quick.elbow(eig, low = 0.08, max.pc = 0.7)
  
  #par(pty = "s",mfrow = c(1,1))
  #plot(eig, xlab = 'PC', ylab = "Variance explained")
  #abline(v = K, col = "red", lty = "dashed")
  
  return(K)
}

# determine best K using tracy widom test
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
  if (K == 0) K <- 1
  
  #plot(eig)
  #abline(v = K)
  
  return(K)
}

# determine best K based on TESS
get_K_tess <- function(gen, coords, Kvals = 1:9, tess_method = "projected.ls", ploidy = 2){
  # coordinates must be a matrix
  coords <- as.matrix(coords)
  
  # Run tess for all K values
  tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy)
  
  # Plot CV results and mark the K-value that is automatically selected
  #plot(tess3_obj, pch = 19, col = "blue",
       #xlab = "Number of ancestral populations",
       #ylab = "Cross-validation score")
  
  # Get best K value
  K <-  bestK(tess3_obj, Kvals)
  
  return(K)
}

#source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
bestK <- function(tess3_obj, Kvals){
  ce <- list()
  for(k in Kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for(k in Kvals) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-max(Kvals)]
  slope <- exp(-diff) - 1
  # K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  return(K)
}

# calculate confusion matrix stats for LFMM
lfmm_calc_confusion <- function(padj, genmat, envmat, lfmm_mod, loci_trait1, loci_trait2, sig = 0.05, all = FALSE){
  
  #performs association testing using the fitted model:
  pv <- lfmm_test(Y = genmat, 
                  X = envmat, 
                  lfmm = lfmm_mod, 
                  calibrate = "gif")
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  # adjust pvalues (or passs through if padj = "none")
  pvalues <-  data.frame(env1 = p.adjust(pv$calibrated.pvalue[,1], method = padj),
                         env2 = p.adjust(pv$calibrated.pvalue[,2], method = padj))
  rownames(pvalues) <- colnames(genmat)
  
  if (all) {
    loci_trait <- c(loci_trait1, loci_trait2)
    
    #env1 candidate loci
    #Identify LFMM cand loci (P)
    lfmm_loci <- unique(c(which(pvalues[,"env1"] < sig), which(pvalues[,"env2"] < sig)))
    #Identify negatives
    lfmm_neg <- (1:nrow(pvalues))[!(1:nrow(pvalues) %in% lfmm_loci)]
    #check length makes sense
    stopifnot(length(lfmm_loci) + length(lfmm_neg) == nrow(pvalues))
    
    #get confusion matrix values
    #True Positives
    TP <- sum(lfmm_loci %in% loci_trait)
    #False Positives
    FP <- sum(lfmm_loci %notin% loci_trait)
    #True Negatives
    TN <- sum(lfmm_neg %notin% loci_trait)
    #False Negatives
    FN <- sum(lfmm_neg %in% loci_trait)
    #check sum makes sense
    stopifnot(sum(TP, FP, TN, FN) == nrow(pvalues))
    
    #calc True Positive Rate (i.e. Sensitivity)
    TPRCOMBO <- TP/(TP + FN)
    #calc True Negative Rate (i.e. Specificity)
    TNRCOMBO <- TN/(TN + FP)
    #calc False Discovery Rate 
    FDRCOMBO <- FP/(FP + TP)
    # correct FDR if 0 is in the denom
    if ((FP + TP == 0) & FP == 0) FDRCOMBO <- 0
    #calc False Positive Rate 
    FPRCOMBO <- FP/(FP + TN)
    
    df <- 
      data.frame(padj = padj,
                 sig = sig,
                 TPRCOMBO = TPRCOMBO, 
                 TNRCOMBO = TNRCOMBO,
                 FDRCOMBO = FDRCOMBO, 
                 FPRCOMBO = FPRCOMBO,
                 TOTALN = length(lfmm_loci), 
                 TOTALTP = TP, 
                 TOTALFP = FP, 
                 TOTALTN = TN,
                 TOTALFN = FN, 
                 all = TRUE)
    
  } else {
    # calculate roc and pr 
    pr1 <- PRROC::pr.curve(scores.class0 = na.omit(pvalues$env1[loci_trait1]), scores.class1 = na.omit(pvalues$env1[-loci_trait1]))$auc.integral
    pr2 <- PRROC::pr.curve(scores.class0 = na.omit(pvalues$env2[loci_trait2]), scores.class1 = na.omit(pvalues$env2[-loci_trait2]))$auc.integral
    roc1 <- PRROC::roc.curve(scores.class0 = na.omit(pvalues$env1[loci_trait1]), scores.class1 = na.omit(pvalues$env1[-loci_trait1]))$auc
    roc2 <- PRROC::roc.curve(scores.class0 = na.omit(pvalues$env2[loci_trait2]), scores.class1 = na.omit(pvalues$env2[-loci_trait2]))$auc
    
    #env1 candidate loci
    #Identify LFMM cand loci (P)
    lfmm_loci1 <- which(pvalues[,"env1"] < sig) 
    #Identify negatives
    lfmm_neg1 <- which(pvalues[,"env1"] >= sig | is.na(pvalues[,"env1"]))
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
    lfmm_loci2 <- which(pvalues[,"env2"] < sig) 
    #Identify negatives
    lfmm_neg2 <- which(pvalues[,"env2"] >= sig | is.na(pvalues[,"env2"]))
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
    lfmm_loci <- unique(c(lfmm_loci1, lfmm_loci2))
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
    # correct FDR if 0 is in the denom
    if ((FP + TP == 0) & FP == 0) FDRCOMBO <- 0
    #calc False Positive Rate 
    FPRCOMBO <- FP/(FP + TN)
    
    # Calculate empirical pvalues (I THINK - CHECK THIS)
    # Get B values (fixed effect)
    Bvalues <-  data.frame(env1 = abs(lfmm_mod$B[,1]), env2 = abs(lfmm_mod$B[,2]))
    null1 <- Bvalues$env1[-loci_trait1]
    emp1 <- sapply(Bvalues$env1[loci_trait1], function(x){mean(x > null1, na.rm = TRUE)})
    emp1_mean <- mean(emp1, na.rm = TRUE)
    null2 <- Bvalues$env2[-loci_trait2]
    emp2 <- sapply(Bvalues$env2[loci_trait2], function(x){mean(x > null2, na.rm = TRUE)})
    emp2_mean <- mean(emp2, na.rm = TRUE)
    EMPCOMBO <- mean(c(emp1_mean, emp2_mean), na.rm = TRUE)
    
    df <- 
      data.frame(padj = padj,
                 sig = sig,
                 TPRCOMBO = TPRCOMBO, 
                 TNRCOMBO = TNRCOMBO,
                 FDRCOMBO = FDRCOMBO, 
                 FPRCOMBO = FPRCOMBO,
                 TOTALN = length(lfmm_loci), 
                 TOTALTP = TP, 
                 TOTALFP = FP, 
                 TOTALTN = TN,
                 TOTALFN = FN,
                 emp1_mean = emp1_mean,
                 emp2_mean = emp2_mean,
                 EMPCOMBO = EMPCOMBO,
                 pr1 = pr1,
                 pr2 = pr2,
                 roc1 = roc1,
                 roc2 = roc2,
                 all = FALSE)
  }
  
  return(df)
}


############
#   RDA    #
############

# run RDA
run_rda <- function(gen, gsd_df, loci_df = NULL, correctPC = c(TRUE, FALSE), maf = c(0, 0.05)){
  # get loci
  if (is.null(loci_df)) loci_df <- get_loci()
  
  # run RDA and pRDA
  result <-
    expand.grid(correctPC = correctPC, maf = maf) %>%
    pmap(\(correctPC, maf) run_rda_helper(gen = gen, gsd_df = gsd_df, loci_df = loci_df, correctPC = correctPC, maf = maf)) %>%
    bind_rows()
  
  return(result)
}

# helper function for run_rda
run_rda_helper <- function(gen, gsd_df, loci_df, sig, correctPC, maf){
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 
  loci_trait2 <- loci_df$trait2 
  
  # maf filter
  if (maf > 0){
    # get colnames of loci
    loci_trait1 <- colnames(gen)[loci_trait1]
    loci_trait2 <- colnames(gen)[loci_trait2]
    
    # perform maf filtering
    gen_maf <- map_dbl(1:ncol(gen), ~mean(gen[,.x], na.rm = TRUE))/2
    gen <- gen[,gen_maf > maf & gen_maf < (1 - maf)]
    
    # reassign loci names with col numbers
    loci_trait1 <- which(colnames(gen) %in% loci_trait1)
    loci_trait2 <- which(colnames(gen) %in% loci_trait2)
  }
  
  # Run RDA
  mod <- rda(gen ~ gsd_df$env1 + gsd_df$env2, scale = F)
  
  # correct PC
  if(correctPC){
    pcres <- prcomp(gen)
    gsd_df$PC1 <-  pcres$x[,1]
    gsd_df$PC2 <-  pcres$x[,2]
    mod <- rda(gen ~ env1 + env2 +  Condition(PC1 + PC2), data = gsd_df, scale = F)
  }
  
  # load scores and get pvalues
  naxes <- ncol(mod$CCA$v)
  rdadapt_env <- rdadapt(mod, naxes)
  
  # P-values
  pv <- rdadapt_env$p.values
  
  # r values
  rv <- 
    rda_cor(gen, gsd_df[,c("env1", "env2")])  %>% 
    dplyr::select(r, snp, var) %>% mutate(r = abs(r)) %>% 
    pivot_wider(names_from = "var", values_from = "r") %>%
    mutate(var = apply(.[,c("env1", "env2")], 1, which.max))
  
  ## make into assignments for each var based on which has the highest correlation
  
  # NOTE: https://github.com/Capblancq/RDA-landscape-genomics used a bonferonni correction
  # correct pvals and get confusion matrix stats
  # change padjustment, alpha level, minimum minor allele frequency, and whether to test for all adaptive loci (i.e. combined trait1/2)
  combos <-
    expand.grid(
      padj = c("none", "fdr", "holm", "bonferroni"),
      sig = c(0.05, 0.10),
      all = c(TRUE, FALSE)
    )
  
  result <-
    pmap(
      combos,
      \(padj, sig, all) return(rda_calc_confusion(
        padj = as.character(padj),
        pv = pv,
        rv = rv,
        loci_trait1 = loci_trait1,
        loci_trait2 = loci_trait2,
        sig = sig,
        all = all
      ))
    )
  
  
  df <- data.frame(correctPC = correctPC, maf = maf, bind_rows(result))
  
  return(df)
}


# Genotype-environment correlation test
rda_cor <- function(gen, var){
  cor_df <- purrr::map_dfr(colnames(gen), rda_cor_env_helper, gen, var)
  rownames(cor_df) <- NULL
  colnames(cor_df) <- c("r", "p", "snp", "var")
  return(cor_df)
}

# Helper function for rda_cor_test
rda_cor_env_helper <- function(snp_name, snp_df, env){
  cor_df <- data.frame(t(apply(env, 2, rda_cor_helper, snp_df[,snp_name])))
  cor_df$snp <- snp_name
  cor_df$env <- colnames(env)
  return(cor_df)
}

# Helper function for rda_cor_test
rda_cor_helper <- function(envvar, snp){
  if(sum(!is.na(envvar)) < 3 | sum(!is.na(snp)) < 3) return(c(r = NA, p = NA))
  # kendall is used instead of pearson because it is non-parameteric and doesn't require vars to be continuous
  mod <- stats::cor.test(envvar, snp, alternative = "two.sided", method = "kendall", na.action = "na.omit")
  pvalue <- mod$p.value
  r <- mod$estimate
  results <- c(r, pvalue)
  names(results) <- c("r", "p")
  return(results)
}

# calculate confusion matrix stats for RDA
rda_calc_confusion <- function(padj = "fdr", sig = 0.05, all = FALSE, pv, rv, loci_trait1, loci_trait2){
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  # adjust pvalues (or pass through if padj = "none")
  pvalues <- data.frame(p = p.adjust(pv, method = padj))
  
  # calculate TPR/FDR based only on RDA p-values (used for comparison with TPR/FDR that includes r)
  rda_loci <- which(pvalues$p < sig)
  loci_trait <- c(loci_trait1, loci_trait2)
  TPRCOMBO_rdap <- sum(rda_loci %in% loci_trait)/length(loci_trait)
  FDRCOMBO_rdap <- sum(rda_loci %notin% loci_trait)/length(rda_loci)

  if (all) {
    #Identify rda cand loci (P)
    rda_loci1 <- which(pvalues$p < sig & rv$var == 1)
    rda_loci2 <- which(pvalues$p < sig & rv$var == 2) 
    
    rda_loci <- unique(c(rda_loci1, rda_loci2))
    
    #Identify negatives
    rda_neg <- which(!(1:nrow(pvalues) %in% rda_loci))
    #check length makes sense
    stopifnot(length(rda_loci) + length(rda_neg) == nrow(pvalues))
    
    #get confusion matrix values
    #True Positives
    TP <- sum(rda_loci %in% loci_trait)
    #False Positives
    FP <- sum(rda_loci %notin% loci_trait)
    #True Negatives
    TN <- sum(rda_neg %notin% loci_trait)
    #False Negatives
    FN <- sum(rda_neg %in% loci_trait)
    #check sum makes sense
    stopifnot(sum(TP, FP, TN, FN) == nrow(pvalues))
    
    #calc True Positive Rate (i.e. Sensitivity)
    TPRCOMBO <- TP/(TP + FN)
    #calc True Negative Rate (i.e. Specificity)
    TNRCOMBO <- TN/(TN + FP)
    #calc False Discovery Rate 
    FDRCOMBO <- FP/(FP + TP)
    # correct FDR if 0 is in the denom
    if ((FP + TP == 0) & FP == 0) FDRCOMBO <- 0
    #calc False Positive Rate 
    FPRCOMBO <- FP/(FP + TN)
    
    df <- 
      data.frame(padj = padj,
                 sig = sig,
                 TPRCOMBO_rdap = TPRCOMBO_rdap,
                 FDRCOMBO_rdap = FDRCOMBO_rdap,
                 TPRCOMBO = TPRCOMBO, 
                 TNRCOMBO = TNRCOMBO,
                 FDRCOMBO = FDRCOMBO, 
                 FPRCOMBO = FPRCOMBO,
                 TOTALN = length(rda_loci), 
                 TOTALTP = TP, 
                 TOTALFP = FP, 
                 TOTALTN = TN,
                 TOTALFN = FN, 
                 all = TRUE)
    
  } else {
    
    #env1 candidate loci
    #Identify rda cand loci (P)
    rda_loci1 <- which(pvalues$p < sig & rv$var == 1) 
    #Identify negatives
    rda_neg1 <- which(!(pvalues$p < sig & rv$var == 1))
    #check length makes sense
    stopifnot(length(rda_loci1) + length(rda_neg1) == nrow(pvalues))
    
    #get confusion matrix values
    #True Positives
    TP1 <- sum(rda_loci1 %in% loci_trait1)
    #False Positives
    FP1 <- sum(rda_loci1 %notin% loci_trait1)
    #True Negatives
    TN1 <- sum(rda_neg1 %notin% loci_trait1)
    #False Negatives
    FN1 <- sum(rda_neg1 %in% loci_trait1)
    #check sum makes sense
    stopifnot(sum(TP1, FP1, TN1, FN1) == nrow(pvalues))
    
    #env2 candidate loci
    #Identify rda cand loci
    rda_loci2 <- which(pvalues$p < sig & rv$var == 2) 
    #Identify negatives
    rda_neg2 <- which(!(pvalues$p < sig & rv$var == 2))
    #check length makes sense
    stopifnot(length(rda_loci2) + length(rda_neg2) == nrow(pvalues))
    
    #True Positives
    TP2 <- sum(rda_loci2 %in% loci_trait2)
    #False Positives
    FP2 <- sum(rda_loci2 %notin% loci_trait2)
    #True Negatives
    TN2 <- sum(rda_neg2 %notin% loci_trait2)
    #False Negatives
    FN2 <- sum(rda_neg2 %in% loci_trait2)
    #check length makes sense
    stopifnot(sum(TP2, FP2, TN2, FN2) == nrow(pvalues))
    
    #stats for all loci 
    rda_loci <- unique(c(rda_loci1, rda_loci2))
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
    if ((FP + TP == 0) & FP == 0) FDRCOMBO <- 0
    #calc False Positive Rate 
    FPRCOMBO <- FP/(FP + TN)
    
    # calculate roc and pr 
    pr <- PRROC::pr.curve(scores.class0 = na.omit(pvalues$p[c(loci_trait1, loci_trait2)]), scores.class1 = na.omit(pvalues$p[-c(loci_trait1, loci_trait2)]))$auc.integral
    roc <- PRROC::roc.curve(scores.class0 = na.omit(pvalues$p[c(loci_trait1, loci_trait2)]), scores.class1 = na.omit(pvalues$p[-c(loci_trait1, loci_trait2)]))$auc
    
    df <- 
      data.frame(padj = padj,
                 sig = sig,
                 TPRCOMBO_rdap = TPRCOMBO_rdap,
                 FDRCOMBO_rdap = FDRCOMBO_rdap,
                 TPRCOMBO = TPRCOMBO, 
                 TNRCOMBO = TNRCOMBO,
                 FDRCOMBO = FDRCOMBO, 
                 FPRCOMBO = FPRCOMBO,
                 TOTALN = length(rda_loci), 
                 TOTALTP = TP, 
                 TOTALFP = FP, 
                 TOTALTN = TN,
                 TOTALFN = FN,
                 pr = pr,
                 roc = roc, 
                 all = FALSE)
  
  }
  
  return(df)
}

# Function to conduct a RDA based genome scan from Capblancq & Forester 2021
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd
rdadapt <- function(rda,K)
{
  zscores <- rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
