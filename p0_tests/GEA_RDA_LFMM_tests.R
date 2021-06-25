library("here") #paths
library("vegan") #RDA
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")


############
#   Data   #
############
#define nloci 
nloci = 10000

#read in geospatial data
#file_path = here("data","mod-10k_K1_phi10_m1_seed1_H50_r60_it--1_t-1000_spp-spp_0.csv")
#file_path = here("data","mod-10k_K1_phi10_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.csv")
#file_path = here("data","mod-10k_K1_phi10_m1_seed1_H50_r60_it--1_t-500_spp-spp_0.csv")
#file_path = here("data","mod-10k_K5_phi50_m1_seed1_H50_r60_it--1_t-300_spp-spp_0.csv")
file_path = here("data","mod-10k_K5_phi10_m1_seed1_H50_r60_it--1_t-300_spp-spp_0.csv")
#file_path = here("data","mod-10k_K5_phi10_m0.25_seed1_H50_r60_it--1_t-300_spp-spp_0.csv")
gsd_df <- read.csv(file_path)
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
head(gsd_df)

#read in genetic data
#file_path = here("data","mod-10k_K1_phi10_m1_seed1_H50_r60_it--1_t-1000_spp-spp_0.vcf")
#file_path = here("data","mod-10k_K1_phi10_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.vcf")
#file_path = here("data","mod-10k_K1_phi10_m1_seed1_H50_r60_it--1_t-500_spp-spp_0.vcf")
#file_path = here("data","mod-10k_K5_phi50_m1_seed1_H50_r60_it--1_t-300_spp-spp_0.vcf")
file_path = here("data","mod-10k_K5_phi10_m1_seed1_H50_r60_it--1_t-300_spp-spp_0.vcf")
#file_path = here("data","mod-10k_K5_phi10_m0.25_seed1_H50_r60_it--1_t-300_spp-spp_0.vcf")
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
gen <- gen[s,]


palz <- magma(100)
par(pty="s",mfrow=c(1,2))
tmpcol<- palz[as.numeric(cut(gea_df$env1,breaks = 100))]
plot(gea_df$x, gea_df$y, col=tmpcol, pch = 19, cex=1.5, main = "env1", xlab="", ylab="", box=TRUE)
tmpcol<- palz[as.numeric(cut(gea_df$env2,breaks = 100))]
plot(gea_df$x, gea_df$y, col=tmpcol, pch = 19, cex=1.5, main = "env2", xlab="", ylab="", box=TRUE)



########
# LFMM #
########

#PCA to determine number of latent factors
#pc <- prcomp(gea_df[,1:nloci])
#par(pty="s",mfrow=c(1,1))
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
K <- 10 #NUMBER OF LATENT FACTORS (NEED TO MODIFY TO MAKE AUTO)

#gen matrix
genmat = as.matrix(gen)
#env matrix
env1mat = as.matrix(gea_df[,"env1"])
env2mat = as.matrix(gea_df[,"env2"])
envmat = cbind(env1mat, env2mat)


##############
#RIDGE METHOD#
##############

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

#OUTPUT RESULTS
#TBD - need to decide on file structure

#PLOT TO CHECK RESULTS

par(mfrow=c(1,2))
plot(-log10(pvalues[,1]), 
     pch = 19, 
     cex = .5, 
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
     cex = .5, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey",
     main = "env2")
points(loci_trait2, 
       -log10(pvalues[,2])[loci_trait2], 
       col = "red", 
       cex = 1.5)
abline(h = -log10(0.05), col="red", lty=2)



#ENV1
#run model
lfmm_mod <- lfmm_ridge(genmat, env1mat, K = K)


#performs association testing using the fitted model:
pv <- lfmm_test(Y = genmat, 
                X = env1mat, 
                lfmm = lfmm_mod, 
                calibrate = "gif")

#adjust pvalues
pvalues <- data.frame(env1=p.adjust(pv$calibrated.pvalue[,1], method="fdr"))


#PLOT TO CHECK RESULTS

plot(-log10(pvalues[,1]), 
     pch = 19, 
     cex = .5, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey",
     main = "env1")
points(loci_trait1, 
       -log10(pvalues[,1])[loci_trait1], 
       col = "red", 
       cex = 1.5)
abline(h = -log10(0.05), col="red", lty=2)




#ENV2
#run model
lfmm_mod <- lfmm_ridge(genmat, env2mat, K = K)


#performs association testing using the fitted model:
pv <- lfmm_test(Y = genmat, 
                X = env2mat, 
                lfmm = lfmm_mod, 
                calibrate = "gif")

#adjust pvalues
pvalues <- data.frame(env1=p.adjust(pv$calibrated.pvalue[,1], method="fdr"))


#PLOT TO CHECK RESULTS

plot(-log10(pvalues[,1]), 
     pch = 19, 
     cex = .5, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey",
     main = "env2")
points(loci_trait2, 
       -log10(pvalues[,1])[loci_trait2], 
       col = "red", 
       cex = 1.5)
abline(h = -log10(0.05), col="red", lty=2)












