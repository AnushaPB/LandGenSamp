library(here) #paths
library(vwgan) #RDA
library(vcfR)  #read VCF files
#parallel
library(foreach)
library(doParallel)

#read in general functions and objects
source("general_functions.R")


############
#   RDA    #
############
run_rda <- function(gen, gsd_df, loci_df){

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
  
  #SIGNIFICANCE CALCS (decide later if these are worth calculating)
  #Determine significance of full model
  #signif.full <- anova.cca(mod, parallel = getOption("mc.cores")) # default is permutation=999
  #print(signif.full)
  
  #Determine significance of axes (variables)
  #signif.axis <- anova.cca(mod, by="axis", parallel = getOption("mc.cores"))
  #print(signif.axis)
  
  #Look at VIF
  #vif.cca(mod)
  
  #load scores
  load.rda <- scores(mod, choices=c(1:2), display="species")  #Choices are RDAs (2 vars = 2 RDAs max)
  
  #OUTLIER FUNCTION
  outliers <- function(x,z){
    lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
    x[x < lims[1] | x > lims[2]]               # locus names in these tails
  }
  
  #Define z (default to 2 or 3)
  z = 2
  cand1 <- outliers(load.rda[,1], z = z) 
  cand2 <- outliers(load.rda[,2], z = z)
  
  #Determine number of candidate loci
  ncand <- length(cand1) + length(cand2)
  ncand
  
  #Create dataframes for each 
  cand1 <- cbind.data.frame(rep(1, times=length(cand1)), names(cand1), unname(cand1))
  cand2 <- cbind.data.frame(rep(2, times=length(cand2)), names(cand2), unname(cand2))
  
  colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")
  
  #combine into one DF
  cand <- rbind(cand1, cand2)
  cand$snp <- as.character(cand$snp)
  
  #Remove X and change to numeric for comparison
  rda_loci <- as.numeric(gsub("0_", "", cand$snp))
  
  #Calc True Positive Rate
  TP <- sum(adaptive_loci %in% rda_loci)
  TPR <- TP/length(adaptive_loci)
  
  #Calc False Discovery Rate
  FD <- sum(neutral_loci %in% rda_loci)
  FDR <- FD/length(rda_loci)
  
  return(data.frame(TPR = TPR, FDR = FDR))
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_rda <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  loci_filepath <- create_filepath(i, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run RDA
  if(skip_to_next == FALSE){
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    loci_df <- get_data(i, "loci")
    
    #run model on full data set
    full_result <- run_rda(gen, gsd_df, loci_df)
    result <- data.frame(sampstrat = "full", nsamp = nrow(gsd_df), full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("RDA_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_rda(subgen, subgsd_df, loci_df)
        
        #save and format new result
        sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, sub_result)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
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

stats_out <- cbind.data.frame(params, res_rda)
write.csv(stats_out, "rda_results.csv")
