
############
#   RDA    #
############
run_rda <- function(gen, gsd_df, loci_df, nloci = 10000, sig = 0.05, correctPC = FALSE){
  
  # get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  # Run RDA
  mod <- rda(gen[,1:nloci] ~ gsd_df$env1 + gsd_df$env2, scale=F)
  
  # correct PC
  if(correctPC){
    pcres <- prcomp(gen)
    gsd_df$PC1 <-  pcres$x[,1]
    gsd_df$PC2 <-  pcres$x[,2]
    mod <- rda(gen[,1:nloci] ~ env1 + env2 +  Condition(PC1 + PC2), data = gsd_df, scale = F)
  }
  
  #plot(mod, type="n", scaling=3)
  #points(mod, display="species", pch=20, cex=2, col="gray32", scaling=3) 
  #text(mod, scaling=3, display="bp", col="#0868ac", cex=1)     
  
  # Get RSQ
  RsquareAdj(mod)
  
  # Plot screeplot
  screeplot(mod)
  
  # load scores and get pvalues
  naxes <- ncol(mod$CCA$v)
  rdadapt_env <- rdadapt(mod, naxes)
  
  # P-values
  pv <- rdadapt_env$p.values
  
  # r values
  rv <- rda_cor(gen, gsd_df[,c("env1", "env2")])
  
  # correct pvals and get confusion matrix stats
  # NOTE: https://github.com/Capblancq/RDA-landscape-genomics used a bonferonni correction
  p05 <- purrr::map_dfr(c("none", "fdr", "holm", "bonferroni"), calc_confusion, pv, rv, loci_trait1, loci_trait2, sig = 0.05)
  p10 <- purrr::map_dfr(c("none", "fdr", "holm", "bonferroni"), calc_confusion, pv, rv, loci_trait1, loci_trait2, sig = 0.10)
  df <- rbind.data.frame(p05, p10)
  
  return(df)
}


#' Genotype-environment correlation test
#'
#' @param gen dosage matrix
#' @param var dataframe with predictor variables
#'
#' @return dataframe with r and p-values from correlation test
#' @export
rda_cor <- function(gen, var){
  cor_df <- purrr::map_dfr(colnames(gen), rda_cor_env_helper, gen, var)
  rownames(cor_df) <- NULL
  colnames(cor_df) <- c("r", "p", "snp", "var")
  return(cor_df)
}

#' Helper function for rda_cor_test
#'
#' @export
#' @noRd
rda_cor_env_helper <- function(snp_name, snp_df, env){
  cor_df <- data.frame(t(apply(env, 2, rda_cor_helper, snp_df[,snp_name])))
  cor_df$snp <- snp_name
  cor_df$env <- colnames(env)
  return(cor_df)
}

#' Helper function for rda_cor_test
#'
#' @export
#' @noRd
rda_cor_helper <- function(envvar, snp){
  if(sum(!is.na(envvar)) < 3 | sum(!is.na(snp)) < 3) return(c(r = NA, p = NA))
  # set to kendall because it is non-parameteric 
  mod <- stats::cor.test(envvar, snp, alternative = "two.sided", method = "kendall", na.action = "na.omit")
  pvalue <- mod$p.value
  r <- mod$estimate
  results <- c(r, pvalue)
  names(results) <- c("r", "p")
  return(results)
}

calc_confusion <- function(padj, pv, rv, loci_trait1, loci_trait2, sig = 0.05){
  
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  
  # adjust pvalues and bind with rv (or pass through if padj = "none")
  pvalues <-  data.frame(p = p.adjust(pv, method = padj))

  #env1 candidate loci
  #Identify rda cand loci (P)
  rv1 <- rv[rv$var == "env1", ]
  stopifnot(nrow(rv1) == nrow(pvalues))
  rda_loci1 <- which(pvalues$p < sig & rv1$p < sig) 
  #Identify negatives
  rda_neg1 <- which(!(pvalues$p < sig & rv1$p < sig))
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
  rv2 <- rv[rv$var == "env2", ]
  stopifnot(nrow(rv2) == nrow(pvalues))
  rda_loci2 <- which(pvalues$p < sig & rv2$p < sig) 
  #Identify negatives
  rda_neg2 <- which(!(pvalues$p < sig & rv2$p < sig))
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
  rda_loci <- c(rda_loci1, rda_loci2)
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

  # calculate roc and pr 
  pr <- PRROC::pr.curve(scores.class0 = na.omit(pvalues$p[c(loci_trait1, loci_trait2)]), scores.class1 = na.omit(pvalues$p[-c(loci_trait1, loci_trait2)]))$auc.integral
  roc <- PRROC::roc.curve(scores.class0 = na.omit(pvalues$p[c(loci_trait1, loci_trait2)]), scores.class1 = na.omit(pvalues$p[-c(loci_trait1, loci_trait2)]))$auc
  
  return(data.frame(padj = padj,
                    sig = sig,
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
                    roc = roc))
}

# Function to conduct a RDA based genome scan from Capblancq & Forester 2021
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd
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
