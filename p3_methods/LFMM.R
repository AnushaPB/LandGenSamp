library("here") #paths
library("vegan") #RDA
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")

#################
#   Test Data   #
#################
#define nloci 
nloci = 10000

#read in geospatial data
file_path = here("data","mod-10k_K1_phi10_m1_seed1_H50_r60_it--1_t-500_spp-spp_0.csv")
gsd_df <- read.csv(file_path)
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
head(gsd_df)

#read in genetic data
file_path = here("data","mod-10k_K1_phi10_m1_seed1_H50_r60_it--1_t-500_spp-spp_0.vcf")
vcf <- read.vcfR(file_path)
x <- vcfR2genlight(vcf) #CHECK THIS
gen <- as.matrix(x)

#create gea_df
gea_df <- data.frame(gen,
                     x = gsd_df$x,
                     y = gsd_df$y, 
                     env1 = gsd_df$env1,
                     env2 = gsd_df$env2)



loci_trait1 <- c(1731,4684,4742,6252) + 1 #add one to convert from python to R indexing
loci_trait2 <- c(141,1512,8481,9511) + 1 #add one to convert from python to R indexing
adaptive_loci <- c(loci_trait1, loci_trait2)
neutral_loci <- c(1:nloci)[-adaptive_loci]

#FOR TESTING USE SAMPLE OF 100 (!!!COMMENT OUT!!!)
set.seed(42)
s <- sample(1:nrow(gea_df),1000)
gea_df <- gea_df[s,]

palz <- magma(100)
par(pty="s",mfrow=c(1,2))
tmpcol<- palz[as.numeric(cut(gea_df$env1,breaks = 100))]
plot(gea_df$x, gea_df$y, col=tmpcol, pch = 19, cex=1.5, main = "env1", xlab="", ylab="", box=TRUE)
tmpcol<- palz[as.numeric(cut(gea_df$env2,breaks = 100))]
plot(gea_df$x, gea_df$y, col=tmpcol, pch = 19, cex=1.5, main = "env2", xlab="", ylab="", box=TRUE)





##########
#  Data  #
##########
#Get gen data
get_gen <- function(filepath){
  vcf <- read.vcfR(filepath)
  x <- vcfR2genlight(vcf) #CHECK THIS
  gen <- as.matrix(x)
  return(gen)
}

#Get geospatial data
get_gsd <- function(filepath){
  gsd_df <- read.csv(filepath)
  gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
  gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
  return(gsd_df)
}

#Create dataframe with all variable combos
params <- expand.grid(K = c(2,5),
            phi = c(0.1,0.5),
            m = c(0.25,0.5,1.0),
            seed = c(1,2,3),
            H = c(0.05,0.5),
            r = c(0.3, 0.6))

#define nloci 
nloci = 10000

##########
#  LFMM  #
##########

run_lfmm <- function(gen_filepath, gsd_filepath, loci_filepath){
  
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #get adaptive loci
  loci_df <- read.csv(loci_filepath)
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  #PCA to determine number of latent factors
  pc <- prcomp(gen)
  par(pty="s",mfrow=c(1,1))
  plot(pc$sdev[1:100]^2, xlab = 'PC', ylab = "Variance explained")
  K <- 10 #NUMBER OF LATENT FACTORS (NEED TO MODIFY TO MAKE AUTO)
  
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
  #Identify LFMM cand loci
  lfmm_loci <- which(pvalues[,1] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci %in% loci_trait1)
  TPR1 <- TP/length(loci_trait1)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci %in% neutral_loci) + sum(lfmm_loci %in% loci_trait2)
  FDR1 <- FD/length(lfmm_loci)
  
  #env2 candidate loci
  #Identify LFMM cand loci
  lfmm_loci <- which(pvalues[,2] < 0.05) 
  #calc True Positive Rate
  TP <- sum(lfmm_loci %in% loci_trait2)
  TPR2 <- TP/length(loci_trait2)
  #calc False Discovery Rate 
  FD <- sum(lfmm_loci %in% neutral_loci) + sum(lfmm_loci %in% loci_trait1)
  FDR2 <- FD/length(lfmm_loci)

  #PLOT TO CHECK RESULTS
  
  par(mfrow=c(1,2))
  plot(-log10(pvalues[,1]), 
       pch = 19, 
       cex = .2, 
       xlab = "SNP", ylab = "-Log P",
       col = "grey",
       main = "env1")
  points(loci_trait1, 
         -log10(pvalues[,1])[loci_trait1], 
         col = "red", 
         cex = 1.5)
  abline(h = -log10(0.05), col="red", lty=2)
  
  plot(-log10(pvalues[,2]), 
       pch = 19, 
       cex = .2, 
       xlab = "SNP", ylab = "-Log P",
       col = "grey",
       main = "env2")
  points(loci_trait2, 
         -log10(pvalues[,2])[loci_trait2], 
         col = "red", 
         cex = 1.5)
  abline(h = -log10(0.05), col="red", lty=2)
  
  return(data.frame(TPR1 = TPR1, FDR1 = FDR1, TPR2 = TPR2, FDR2 = FDR2))
}



#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library(here)
  
  gen_filepath <- here("data",paste0("mod-K",params[i,"K"],
                                     "_phi",params[i,"phi"]*100,
                                     "_m",params[i,"m"]*100,
                                     "_seed",params[i,"seed"],
                                     "_H",params[i,"H"]*100,
                                     "_r",params[i,"r"]*100,
                                     "it--1_t-500_spp-spp_0.vcf"))
  
  gsd_filepath <- here("data",paste0("mod-K",params[i,"K"],
                                     "_phi",params[i,"phi"]*100,
                                     "_m",params[i,"m"]*100,
                                     "_seed",params[i,"seed"],
                                     "_H",params[i,"H"]*100,
                                     "_r",params[i,"r"]*100,
                                     "it--1_t-500_spp-spp_0.csv"))
  
  loci_filepath <- here("data",paste0("nnloci_",params[i,"K"],
                                      "_phi",params[i,"phi"]*100,
                                      "_m",params[i,"m"]*100,
                                      "_seed",params[i,"seed"],
                                      "_H",params[i,"H"]*100,
                                      "_r",params[i,"r"]*100,
                                      ".csv"))
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(exists(loci_filepath) == FALSE | gen_filepath == FALSE | gsd_filepath == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
                      print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_lfmm(gen_filepath = gen_filepath,
                       gsd_filepath = gsd_filepath,
                       loci_filepath = loci_filepath)
  }
  
  return(result)
  
}
  
#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)
write.csv(stats_out, "LFMM_results.csv")

