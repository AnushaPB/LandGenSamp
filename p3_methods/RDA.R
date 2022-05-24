set.seed(42)

library(here) #paths
library(vegan) #RDA
library(vcfR)  #read VCF files
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")


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
  mod <- rda(gen[, 1:nloci] ~ gsd_df$env1 + gsd_df$env2, scale=T)
  
  #Get RSQ
  RsquareAdj(mod)
  
  #Plot screeplot
  screeplot(mod)
  
  #load scores and get pvalues
  naxes <- ncol(mod$CCA$v)
  rdadapt_env <- rdadapt(final_mod, naxes)
  
  # P-values threshold after FDR correction (different from Capblancq & Forester 2021)
  pvalues <- p.adjust(rdadapt_env$p.values, method = padj_method)
 
  ## Identifying the loci that are below the p-value threshold
  # NEED TO ADD STEP TO FIGURE OUT SIGNIFICANCE OF ENV VARS
  #Identify rda cand loci (P)
  rda_loci <- which(pvalues < sig) 
  
  #Calc True Positive Rate
  TP <- sum(rda_loci %in% adaptive_loci)
  TPR <- TP/length(adaptive_loci)
  
  #Calc False Discovery Rate
  FD <- sum(neutral_loci %in% rda_loci)
  FDR <- FD/length(rda_loci)
  
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
  
  
  return(data.frame(TPR = TPR, 
                    FDR = FDR,
                    TP = TP,
                    FP = FP,
                    TN = TN,
                    FN = FN))
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

#register cores
cores <- 10
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

res_rda <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("vcfR")
  library("vegan")
  library("here")
  
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

write.csv(res_rda, "outputs/rda_results.csv", row.names = FALSE)
