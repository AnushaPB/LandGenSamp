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

#parallel
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")


############
#   TEST   #
############

#DEFINE NLOCI
nloci = 1000 #CHANGE LATER TO 10k
#NEEDS TO BE MODIFIED FOR FUTURE REAL DATA
gea_df <- read.csv(here("data","gea_m0.5_phi0.5_H0.5_k10_t100_df.csv"))
gea_df <- gea_df[,-1] #PROB CAN REMOVE ONCE GNX SCRIPTS ARE CORRECTED
colnames(gea_df) <- c(paste0("X",1:nloci), colnames(gea_df)[(nloci+1):ncol(gea_df)]) # CHANGE FROM BASE 0 TO BASE 1

loci_df <- read.csv(here("data","loci_m0.5_phi0.5_H0.5_k10_t100_df.csv"))
loci_df <- data.frame(trait0 = loci_df$trait0)
adaptive_loci <- which(loci_df$trait0 == 1) #CURRENTLY DEFINED FOR ONE TRAIT
neutral_loci <- which(loci_df$trait0 == 0) #CURRENTLY DEFINED FOR ONE TRAIT

#FOR TESTING USE SAMPLE OF 100 (!!!COMMENT OUT!!!)
set.seed(42)
s <- sample(1:nrow(gea_df),100)
loci_df <- loci_df[s,neutral_loci]
gea_df <- gea_df[s,]


#########
#  LEA  #
#########
run_lea <- function(gen_filepath, gsd_filepath, loci_filepath){
  
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #get adaptive loci
  loci_df <- read.csv(loci_filepath)
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
    
  
  #remove adaptive loci
  gen <- gen[,-adaptive_loci]
  
  #create gen matrix
  gen <- as.matrix(gen)
  
  #create temporary file with genotypes
  write.geno(gen, here("data","temp_genotypes.geno"))
  
  #Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
  
  #define K based on "true" value
  #MODIFY LATER TO READ FROM FILE (incorporate into main function?)
  K = 4
  
  #Code for testing one K value
  obj.snmf = snmf(here("data","temp_genotypes.geno"), K = K, ploidy = 2, entropy = T, alpha = 100, project = "new")
  
  #Get Qmatrix
  pred_admix <- Q(obj.snmf, K = K) 
  
  pred_krig_admix <- stack()
  for(k in 1:K){
    #krig admix proportions
    krig_df = data.frame(x = gsd_df$x,
                         y = gsd_df$y, 
                         prop = pred_admix[,k])
    
    coordinates(krig_df)=~x+y
    
    x.range <- c(0,40)
    y.range <- c(0,40)
    
    krig_df.grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], len=40),  
                               y=seq(from=y.range[1], to=y.range[2], len=40))
    coordinates(krig_df.grd) <- ~x+y
    gridded(krig_df.grd) <- TRUE
    #extent(krig_df.grid) #FIGURE OUT WHY EXTENT ISNT (0,40,0,40)
    
    krig_res <- autoKrige(prop ~ 1, krig_df, krig_df.grd)
    pred_krig_admix <- stack(pred_krig_admix, raster(krig_res$krige_output))
  }
  
  plot(pred_krig_admix)
  
  #remove project
  remove.snmfProject("data/temp_genotypes.snmfProject")
  
  
  ###########################
  #  Admix Krig Comparison  #
  ###########################
  
  true_krig_admix <- stack("data/test_pred_krig_admix.tif") #MODIFY LATER DEPENDING ON FILE FORMAT
  
  #correlation between all layers in stacks
  full_cor <- layerStats(stack(pred_krig_admix, true_krig_admix), "pearson")$'pearson correlation coefficient'
  #use the absolute value of the correlation (in case the layers are flipped) (e.g. so a correlation of -1 would be treated like a correlation of 1) (THINK THROUGH THIS)
  full_cor <- abs(full_cor)
  
  #subset of correlation matrix to just compare pred admix to true admix (pred admix is rows, true admix is columns)
  sub_cor <- full_cor[1:K,1:K+K]
  
  #this loop works by first identifying the max correlation between the layers, which presumably should be the K layers that correspond to each other
  #then it removes those layers from the matrix and calculates the next maximum correlation (e.g. the next corresponding K layers)
  #it repeats this process until the final iteration
  #this should get the correlations between the sampe K layers (hopefully)
  krigcor <- c()
  for(k in 1:K){
    #pull out max correlation
    krigcor <- c(krigcor, max(sub_cor))
    #identify index of max correlation
    loc <- which(sub_cor == max(sub_cor), arr.ind=TRUE)
    #remove layers corresponding to max correlation until the last iteration of the loop (when there will be only one number remaining)
    #pred admix is rows, true admix is columns
    if(k != K){sub_cor <- sub_cor[-loc[,"row"], -loc[,"col"]]} 
  }
  print(krigcor)
  
  
  ############################
  #  Admix Coeff Comparison  #
  ############################
  
  true_admix <- read.csv("data/test_pred_admix.csv") #MODIFY LATER DEPENDING ON FILE FORMAT
  #NEED TO WRITE CODE SUCH THAT THE ADMIX COEFFS FROM THE "TRUE" MATRIX ARE FROM THE SAME INDIVIDUALS OF THE SUBSAMPLE (either modify the input file, or use the individual IDS?)
  
  colnames(pred_admix) <- paste0("pred",1:K)
  colnames(true_admix) <- paste0("true",1:K)
  
  rmse <- c()
  for(k in 1:K){
    rmse[k] <- sqrt(mean((pred_admix[,k] - true_admix[,k])^2)) #NEED TO FIGURE OUT HOW TO MAKE SURE K=1 in true is same as K=1 in pred
  }
  
  #correlation between pred and true admix 
  sub_cor <- cor(pred_admix, true_admix)
  #use the absolute value of the correlation (in case the layers are flipped) (e.g. so a correlation of -1 would be treated like a correlation of 1) (THINK THROUGH THIS)
  sub_cor <- abs(sub_cor)
  #this loop works by first identifying the max correlation between the layers, which presumably should be the K layers that correspond to each other
  #then it removes those layers from the matrix and calculates the next maximum correlation (e.g. the next corresponding K layers)
  #it repeats this process until the final iteration
  #this should get the correlations between the sampe K layers (hopefully)
  admixcor <- c()
  for(k in 1:K){
    #pull out max correlation
    admixcor <- c(admixcor, max(sub_cor))
    #identify index of max correlation
    loc <- which(sub_cor == max(sub_cor), arr.ind=TRUE)
    #remove layers corresponding to max correlation until the last iteration of the loop (when there will be only one number remaining)
    #pred admix is rows, true admix is columns
    if(k != K){sub_cor <- sub_cor[-loc[,"row"], -loc[,"col"]]} 
  }
  print(admixcor)
  
  results <- data.frame(krigcor = mean(krigcor), admixcor = mean(admixcor))
  
  return(results)
  

}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lea <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library(here)
  
  gen_filepath <- create_filepath(i, "gen")
  
  gsd_filepath <- create_filepath(i, "gsd")
  
  loci_filepath <- create_filepath(i, "loci")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(exists(loci_filepath) == FALSE | exists(gen_filepath) == FALSE | exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_lea(gen_filepath = gen_filepath,
                      gsd_filepath = gsd_filepath,
                      loci_filepath = loci_filepath)
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)
write.csv(lea_out, "LEA_results.csv")
