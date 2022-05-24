set.seed(42)
#paths
library("here") 
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")
library("AssocTests")

#read in general functions and objects
source("general_functions.R")

##########
#  LFMM  #
##########


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
    K <- get_K(gen, k_selection = "quick.elbow")
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


# function to determine K
get_K <- function(gen, coords = NULL, k_selection = "quick.elbow", Kvals = Kvals, ...){
  
  if(k_selection == "tracy.widom"){K <- get_K_tw(gen)}
  
  if(k_selection == "quick.elbow"){K <- get_K_elbow(gen)}
  
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
  K <- quick.elbow(eig, low = 0.08, max.pc = 0.9)
  
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
  tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = 3.2724)
  
  # get K based on number of significant eigenvalues
  K <- tw_result$SigntEigenL
  
  plot(eig)
  abline(v = K)
  
  return(K)
}


#register cores
cores <- 10
cl <- makeCluster(cores)
#not to overload your computer
registerDoParallel(cl)

system.time(
res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  library("stringr")
  library("AssocTests")

  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #create pdf to store plots
  # pdf(paste0("outputs/LFMM/plots/lfmm_plots_",paramset,".pdf"))
  
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
    full_result <- run_lfmm(gen_2k, gsd_df_2k, loci_df, K = NULL)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = 2000, full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LFMM/LFMM_results_",paramset,".csv")
    write.csv(result, csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_lfmm(subgen, subgsd_df, loci_df, K = full_result$K)
        
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
  
  #end pdf()
  # dev.off()
  
  return(result)
  
  gc()
}
)


#stop cluster
stopCluster(cl)

write.csv(res_lfmm, "outputs/lfmm_results.csv", row.names = FALSE)

