
set.seed(42)

library("here") #paths
library("vcfR")
#to install LEA:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("LEA")
library("LEA") #LEA

#for plotting (from LEA)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#kriging
library("automap")
library("raster")
library("rgdal")

#parallel
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")

#########
#  LEA  #
#########
run_lea_full <- function(gen, gsd_df, loci_df, paramset){
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #remove adaptive loci
  gen <- gen[,-adaptive_loci]
  
  #create gen matrix
  gen <- as.matrix(gen)
  

  #create temporary file with genotypes
  #!names must be unique for running in parallel!
  write.geno(gen, here("data", paste0(paramset,"_temp_genotypes.geno")))
  
  #Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
  
  #Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
  #Code for testing multiple K values:
  maxK <- 10
  obj.snmf <- snmf(here("data", paste0(paramset,"_temp_genotypes.geno")), K = 1:maxK, ploidy = 2, entropy = T, alpha = 100, project = "new")
  
  #determining best K and picking best replicate for best K (source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
  ce <- list()
  for(k in 1:maxK) ce[[k]] <- cross.entropy(obj.snmf, K=k)
  ce.K <- c()
  for(k in 1:maxK) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-maxK]
  slope <- exp(-diff) - 1
  #K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  
  
  plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
  abline(v = K, col = "red", lty = "dashed")

  #Get Qmatrix
  pred_admix <- Q(obj.snmf, K = K) 
  
  #export full krig admix and Qmatrix
  full_admix <- pred_admix
  rownames(full_admix) <- rownames(gen)
  write.csv(full_admix, paste0("outputs/LEA/",paramset,"_qmat.csv"))
  
  #remove project
  remove.snmfProject(here("data", paste0(paramset,"_temp_genotypes.geno")))
  
  results <- data.frame(popK = K, cor = NA, rmse = NA)
  
  return(results)
  
}

run_lea <- function(gen, gsd_df, loci_df, K, full_admix){
  #get adaptive loci
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #remove adaptive loci
  gen <- gen[,-adaptive_loci]
  
  #create gen matrix
  gen <- as.matrix(gen)
  
  #create temporary file with genotypes
  write.geno(gen, here("data", paste0(paramset,"_temp_genotypes.geno")))
  
  #Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
  #Code for running one K value
  obj.snmf = snmf(here("data", paste0(paramset,"_temp_genotypes.geno")), K = K, ploidy = 2, entropy = T, alpha = 100, project = "new")
  
  #Get Qmatrix
  pred_admix <- Q(obj.snmf, K = K) 
  
  #remove project
  remove.snmfProject(here("data", paste0(paramset,"_temp_genotypes.geno")))
  
  ############################
  #  Admix Coeff Comparison  #
  ############################
  
  rownames(pred_admix) <- rownames(full_admix)
  colnames(pred_admix) <- paste0("pred",1:K)
  colnames(full_admix) <- paste0("full",1:K)
  
  #TEMP TO DEAL WITH NAS
  full_admix <- full_admix[complete.cases(full_admix),]
  pred_admix <- pred_admix[row.names(full_admix),]
  
  rmse <- c()
  for(k in 1:K){
    rmse[k] <- sqrt(mean((pred_admix[,k] - full_admix[,k])^2)) #NEED TO FIGURE OUT HOW TO MAKE SURE K=1 in true is same as K=1 in pred
  }
  
  #correlation between pred and true admix 
  sub_cor <- cor(pred_admix, full_admix)
  #use the absolute value of the correlation (in case the layers are flipped) (e.g. so a correlation of -1 would be treated like a correlation of 1) (THINK THROUGH THIS)
  sub_cor <- abs(sub_cor)
  
  results <- data.frame(popK = K, cor = mean(sub_cor), rmse = mean(rmse))
  
  return(results)
  

}


#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

res_lea <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  library("LEA")
  library("automap")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #save plots
  pdf(paste0("outputs/LEA/LEA_plots_",paramset))
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  loci_filepath <- create_filepath(i, params = params, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LEA
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    #TEMP FIX: SAMPLE INCLUDES IDS USED IN SUB SAMPLING
    allsubIDs <- c()
    #get list of subIDs
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        allsubIDs <- c(allsubIDs, subIDs)
      }
    }
    allsubIDs <- unique(allsubIDs)
    
    leftoverIDs <- setdiff(row.names(gsd_df), allsubIDs)
    
    #subsample full data randomly from subIDS (plus extra to hit 2000)
    if(length(allsubIDs) < 2000){
      s <- sample(leftoverIDs, 2000 - length(allsubIDs), replace = FALSE)
      s <- c(allsubIDs, s)
      gen_2k <- gen[s,]
      gsd_df_2k <- gsd_df[s,]
    } else {
      gen_2k <- gen[allsubIDs,]
      gsd_df_2k <- gsd_df[allsubIDs,]
    }
    
    
    #run model on full data set
    full_result <- run_lea_full(gen_2k, gsd_df_2k, loci_df, paramset)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = nrow(gsd_df_2k), full_result)
    
    #get admix and krig_admix from full data
    full_admix <- read.csv(paste0("outputs/LEA/",paramset,"_qmat.csv"), row.names = 1)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LEA/LEA_results_",paramset,".csv")
    write.csv(result, csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        subfull_admix <- full_admix[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_lea(subgen, subgsd_df, loci_df, K = full_result$popK, full_admix = subfull_admix)
        
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
  
  dev.off()
  
  return(result)
  
  gc()
  
}

#stop cluster
stopCluster(cl)

write.csv(res_lea, "outputs/LEA/lea_results.csv", row.names = FALSE)
