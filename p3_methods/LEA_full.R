library("here") #paths
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

gen <- as.matrix(gen)
#create temporary file with genotypes
write.geno(gen, here("data","temp_genotypes.geno"))

#Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
#Code for testing mulitple K values:
obj.snmf = snmf(here("data","temp_genotypes.geno"), K = 1:10, ploidy = 2, entropy = T, alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
#Code for testing one K value
obj.snmf = snmf(here("data","temp_genotypes.geno"), K = K, ploidy = 2, entropy = T, alpha = 100, project = "new")

#define K based on "true" value
#MODIFY LATER
K = 3


#Get Qmatrix
pred_admix <- Q(obj.snmf, K = K) 


pred_krig_admix <- stack()
for(i in 1:K){
  #krig admix proportions
  krig_df = data.frame(x = gea_df$x,
                       y = gea_df$y, 
                       prop = pred_admix[,i])
  
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

############################
#  Admix Coeff Comparison  #
############################

true_admix <- read.csv() #MODIFY LATER DEPENDING ON FILE FORMAT

rmse <- c()
for(i in 1:K){
  rmse[i] <- sqrt(mean((pred_admix[,i] - true_admix[,i])^2)) #NEED TO FIGURE OUT HOW TO MAKE SURE K=1 in true is same as K=1 in pred
}

#ADD OUTPUT

###########################
#  Admix Krig Comparison  #
###########################

true_krig_admix <- stack() #MODIFY LATER DEPENDING ON FILE FORMAT

krigcor <- c()
for(i in 1:K){
  krigcor[i] <- layerStats(stack(pred_krig_admix[[i]], true_krig_admix[[i]]), "pearson")$'pearson correlation coefficient'[1,2]
}

remove.snmfProject("data/temp_genotypes.snmfProject")
#ADD OUTPUT\
