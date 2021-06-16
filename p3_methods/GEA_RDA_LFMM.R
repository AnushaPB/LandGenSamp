library("here") #paths
library("vegan") #RDA
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM



############
#   Data   #
############
#define nloci 
nloci = 10000

#read in geospatial data
file_path = ""
gsd_df <- read.csv(file_path)
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=\\[)[^,]+(?=,)')) #CHECK THIS
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) #CHECK THIS
head(gsd_df)

#read in genetic data
file_path = ""
vcf <- read.vcfR(file_path)
x <- vcfR2genlight(vcf) #CHECK THIS
gen <- as.matrix(x)

#create gea_df
gea_df <- data.frame(gen,
                     x = gsd_df$x,
                     y = gsd_df$y, 
                     env1 = gsd_df$env1,
                     env2 = gsd_df$env2)

#read in adaptive loci
file_path = ""
loci_df <- read.csv(file_path)
loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
adaptive_loci <- c(loci_trait1, loci_trait2)
neutral_loci <- c(1:nloci)[-adaptive_loci]

#FOR TESTING USE SAMPLE OF 100 (!!!COMMENT OUT!!!)
set.seed(42)
s <- sample(1:nrow(gea_df),100)
gea_df <- gea_df[s,]


############
#   RDA    #
############

#Run RDA
mod <- rda(gea_df[, 1:nloci] ~ env1 + env2, data=gea_df, scale=T)

#Get RSQ
RsquareAdj(mod)

#Plot screeplot
screeplot(mod)

#Determine significance of full model
signif.full <- anova.cca(mod, parallel = getOption("mc.cores")) # default is permutation=999
signif.full

#Determine significance of axes (variables)
signif.axis <- anova.cca(mod, by="axis", parallel = getOption("mc.cores"))
signif.axis

#Look at VIF
#vif.cca(mod)

#load scores
load.rda <- scores(mod, choices=c(1:2), display="species")  #Choices are RDAs (2 vars = 2 RDAs max)

#OUTLIER FUNCTION
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#Define z (default to 3)
z = 3
cand1 <- outliers(load.rda[,1], z = z) 
cand2 <- outliers(load.rda[,2], z = z)

#Determine number of candidate loci
ncand <- length(cand1) + length(cand2)
ncand

#Create dataframes for each 
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")

#combine into one DF
cand <- rbind(cand1, cand2)
cand$snp <- as.character(cand$snp)

#Remove X and change to numeric for comparison
rda_loci <- as.numeric(gsub("X", "", cand$snp))

#Calc True Positive Rate
TP <- sum(adaptive_loci %in% rda_loci)
TPR <- TP/length(adaptive_loci)

#Calc False Discovery Rate
FD <- sum(neutral_loci %in% rda_loci)
FDR <- FD/length(rda_loci)

#OUTPUT RESULTS
#TBD - need to decide on file structure

########
# LFMM #
########

#PCA to determine number of latent factors
pc <- prcomp(gea_df[,1:nloci])
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
K <- 6 #NUMBER OF LATENT FACTORS (NEED TO MODIFY TO MAKE AUTO)

#gen matrix
genmat = as.matrix(gea_df[,1:nloci])
#env matrix
env1mat = as.matrix(gea_df[,"env1"])
env2mat = as.matrix(gea_df[,"env2"])
envmat = cbind(env1mat, env2mat)

#run model
lfmm_mod <- lfmm_ridge(genmat, envmat, K = K)

#performs association testing using the fitted model:
pv <- lfmm_test(Y = genmat, 
                X = envmat, 
                lfmm = lfmm_mod, 
                calibrate = "gif")

#define pvalue (simple correction - MODIFY? very strict)
padj = 0.05/nloci

#Identify LFMM cand loci
lfmm_loci <- which(pv$calibrated.pvalue < padj) 

#calc True Postive Rate
TP <- sum(lfmm_loci %in% adaptive_loci)
TPR <- TP/length(adaptive_loci)

#calc False Discovery Rate 
FD <- sum(lfmm_loci %in% neutral_loci)
FDR <- FD/length(lfmm_loci)

#OUTPUT RESULTS
#TBD - need to decide on file structure

#PLOT TO CHECK RESULTS
pvalues <- pv$calibrated.pvalue[,1]
plot(-log10(pvalues), 
     pch = 19, 
     cex = .2, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey")
points(adaptive_loci, 
       -log10(pvalues)[adaptive_loci], 
       col = "red", 
       cex = 1.5)
