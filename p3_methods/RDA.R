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

#read in adaptive loci
file_path = here("data","nnloci_10k_K1_phi10_m1_seed1_H50_r60.csv")
loci_df <- read.csv(file_path)
loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
adaptive_loci <- c(loci_trait1, loci_trait2)
neutral_loci <- c(1:nloci)[-adaptive_loci]

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
rda_loci <- as.numeric(gsub("X0_", "", cand$snp))

#Calc True Positive Rate
TP <- sum(adaptive_loci %in% rda_loci)
TPR <- TP/length(adaptive_loci)

#Calc False Discovery Rate
FD <- sum(neutral_loci %in% rda_loci)
FDR <- FD/length(rda_loci)

#OUTPUT RESULTS
#TBD - need to decide on file structure

