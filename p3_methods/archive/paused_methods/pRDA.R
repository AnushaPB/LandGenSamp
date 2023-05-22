set.seed(42)

library(here) #paths
library(vegan) #RDA
library(vcfR)  #read VCF files
#parallel
library(foreach)
library(doParallel)
#RDA
#devtools::install_github("jdstorey/qvalue")
library(qvalue)
library(robust)

#read in general functions and objects
source("general_functions.R")


############
#   RDA    #
############


run_rda <- function(gen, gsd_df, loci_df, K=NULL){

  #define nloci
  nloci = ncol(gen)
  
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #create formula such that it can change
  f <- as.formula("gen[, 1:nloci] ~ gsd_df$env1 + gsd_df$env2")
  mod <- rda(f)
  
  #### Function to conduct a RDA based genome scan from Capblancq & Forester 2021
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
  
  ## Running the function with two axes (two env variables)
  rdadapt_env <- rdadapt(mod, 2)
   
  ## P-values threshold after FDR correction (different from Capblancq & Forester)
  pvalues <- p.adjust(rdadapt_env$p.values, method="fdr")
  #Capblancq include a step where they only take pvalues with highest loading for each contig to deal with LD (not applied here)
  
  ## Identifying the loci that are below the p-value threshold
  ##ADD THIS: . Selected loci can then be tested for their association with proposed environmental drivers using simple correlations  with  allele  frequencies,  or  permutations  (e.g.,  Pavlova  et al., 2013).
  #Identify rda cand loci (P)
  rda_loci <- which(pvalues < 0.05) 
  #Identify negatives
  rda_neg <- which(pvalues >= 0.05 | is.na(pvalues))
  #for readibility, just negates the in function
  `%notin%` <- Negate(`%in%`)
  #get confusion matrix values
  #True Positives
  TP <- sum(rda_loci %in% adaptive_loci)
  #False Positives
  FP <- sum(rda_loci %notin% adaptive_loci)
  #True Negatives
  TN <- sum(rda_neg %notin% adaptive_loci)
  #False Negatives
  FN <- sum(rda_neg %in% adaptive_loci)
  
  #calc True Positive Rate
  TPRCOMBO <- TP/(TP + FN)
  #calc True Negative Rate
  TNRCOMBO <- TN/(TN + FP)
  #calc False Discovery Rate 
  FDRCOMBO <- FP/(FP + TP)
  #calc False Positive Rate
  FPRCOMBO <- FP/(FP + TN)
  
  #VARIANCE EXPLAINED CODE (placeholder - might add in later)
  ## Full model
  #full <- as.formula(paste("gen[, 1:nloci] ~ gsd_df$env1 + gsd_df$env2 + gsd_df$x + gsd_df$y +", pcform))
  #pRDAfull <- rda(full)
  #RsquareAdj(pRDAfull)
  #anova(pRDAfull)
  ## Pure climate model
  #env <- as.formula(paste("gen[, 1:nloci] ~ gsd_df$env1 + gsd_df$env2 + Condition(gsd_df$x + gsd_df$y +", pcform, ")"))
  #pRDAenv<- rda(env)
  #RsquareAdj(pRDAenv)
  #anova(pRDAenv)
  ## Pure neutral population structure model  
  #struct <- as.formula(paste("gen[, 1:nloci] ~ ", pcform, "+ Condition(gsd_df$env1 + gsd_df$env2 + gsd_df$x + gsd_df$y)"))
  #pRDAstruct <- rda(struct),  Variables)
  #RsquareAdj(pRDAstruct)
  #anova(pRDAstruct)
  ##Pure geography model 
  #struct <- as.formula(paste("gen[, 1:nloci] ~ gsd_df$x + gsd_df$y + Condition(gsd_df$env1 + gsd_df$env2 +", pcform, ")"))
  #pRDAgeog <- rda(AllFreq ~ Longitude + Latitude + Condition(MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + PC1 + PC2 + PC3),  Variables)
  #RsquareAdj(pRDAgeog)
  #anova(pRDAgeog)
  
  return(data.frame(
                    TPRCOMBO = TPRCOMBO, 
                    TNRCOMBO = TNRCOMBO,
                    FDRCOMBO = FDRCOMBO, 
                    FPRCOMBO = FPRCOMBO,
                    TOTALN = length(rda_loci), 
                    TOTALTP = TP, 
                    TOTALFP = FP, 
                    TOTALTN = TN,
                    TOTALFN = FN))
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-5) #not to overload your computer
registerDoParallel(cl)

res_rda <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {  
# res_rda <- foreach(i=1:2, .combine=rbind) %dopar% {
  library("vcfR")
  library("vegan")
  library("qvalue")
  library("robust")
  library("here")
  library("stringr")

  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  loci_filepath <- create_filepath(i, params = params, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run RDA
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_rda(gen_2k, gsd_df_2k, loci_df)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = nrow(gsd_df), full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/RDA/RDA_results_",paramset,".csv")
    write.csv(result, csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_rda(subgen, subgsd_df, loci_df)
        
        #save and format new result
        sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsamp, sub_result)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, sub_result)
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  return(result)
  
  gc()
  
}

#stop cluster
stopCluster(cl)

write.csv(res_rda, "outputs/RDA/rda_results.csv", row.names = FALSE)
