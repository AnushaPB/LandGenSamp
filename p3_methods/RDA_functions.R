
############
#   RDA    #
############
run_rda <- function(gen, gsd_df, loci_df, nloci = 10000, sig = 0.05){
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #Run RDA
  mod <- rda(gen[,1:nloci] ~ gsd_df$env1 + gsd_df$env2, scale=T)
  
  #plot(mod, type="n", scaling=3)
  #points(mod, display="species", pch=20, cex=2, col="gray32", scaling=3) 
  #text(mod, scaling=3, display="bp", col="#0868ac", cex=1)     
  
  #Get RSQ
  RsquareAdj(mod)
  
  #Plot screeplot
  screeplot(mod)
  
  #load scores and get pvalues
  naxes <- ncol(mod$CCA$v)
  rdadapt_env <- rdadapt(mod, naxes)
  
  # P-values
  pv <- rdadapt_env$p.values
  
  # correct pvals and get confusion matrix stats
  p05 <- purrr::map_dfr(c("none", "fdr", "holm"), calc_confusion, pv, adaptive_loci, neutral_loci, alpha = 0.05)
  p10 <- purrr::map_dfr(c("none", "fdr", "holm"), calc_confusion, pv, adaptive_loci, neutral_loci, alpha = 0.10)
  df <- rbind.data.frame(p05, p10)
  
  return(df)
}


calc_confusion <- function(padj, pv, adaptive_loci, neutral_loci, alpha = 0.05){
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  # adjust pvalues (or passs through if padj = "none")
  pvalues <- data.frame(env = p.adjust(pv, method = padj))
  
  #env candidate loci
  rda_loci <- which(pvalues[, "env"] < alpha)
  
  #get confusion matrix values
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  #True Positives
  TP <- sum(adaptive_loci %in% rda_loci)
  #False Positives
  FP <- sum(neutral_loci %in% rda_loci)
  #True Negatives
  TN <- sum(neutral_loci %notin% rda_loci)
  #False Negatives
  FN <- sum(adaptive_loci %notin% rda_loci)
  
  #calc True Positive Rate (i.e. Sensitivity)
  TPRCOMBO <- TP/(TP + FN)
  #calc True Negative Rate (i.e. Specificity)
  TNRCOMBO <- TN/(TN + FP)
  #calc False Discovery Rate 
  FDRCOMBO <- FP/(FP + TP)
  #calc False Positive Rate 
  FPRCOMBO <- FP/(FP + TN)
  
  # Calculate relative pvalues
  null <- pvalues$env[-adaptive_loci]
  emp <- sapply(pvalues$env[adaptive_loci], function(x){mean(x > null, na.rm = TRUE)})
  emp_TPR <- sum(emp < alpha, na.rm  = TRUE)
  
  return(data.frame(padj = padj,
                    alpha = alpha,
                    TPRCOMBO = TPRCOMBO, 
                    TNRCOMBO = TNRCOMBO,
                    FDRCOMBO = FDRCOMBO, 
                    FPRCOMBO = FPRCOMBO,
                    TOTALN = length(rda_loci), 
                    TOTALTP = TP, 
                    TOTALFP = FP, 
                    TOTALTN = TN,
                    TOTALFN = FN,
                    emp_TPR = emp_TPR,
                    emp_mean = mean(emp, na.rm = TRUE)))
}

# Function to conduct a RDA based genome scan from Capblancq & Forester 2021
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd
# NOTE: GO THROUGH THIS CODE
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
