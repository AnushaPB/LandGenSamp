## GLMM ##
library(spaMM)
library(statgenGWAS)
#BiocManager::install("LEA")
library("LEA") #LEA

#for plotting (from LEA)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")


#define nloci 
nloci = 10000

#read in geospatial data
file_path = here("data","mod-10k_K1_phi10_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.csv")
gsd_df <- read.csv(file_path)
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
head(gsd_df)

#read in genetic data
file_path = here("data","mod-10k_K1_phi10_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.vcf")
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
#create temporary file with genotypes
write.geno(gen,here("data","temp_genotypes.geno"))

#define K based on "true" value
#MODIFY LATER
K = 20

#Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
#Code for testing mulitple K values:
#obj.snmf = snmf(here("data","temp_genotypes.geno"), K = 1:20, ploidy = 2, entropy = T, alpha = 100, project = "new")
#plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
#Code for testing one K value
obj.snmf = snmf(here("data","temp_genotypes.geno"), K = K, ploidy = 2, entropy = T, alpha = 100, project = "new")

#Get Qmatrix
Qmat <- Q(obj.snmf, K = K) 
#get pops
pops <- colnames(Qmat)[apply(Qmat,1,which.max)]

#make kinship matrix (method = c("astle", "IBS", "vanRaden"))
Kmatrix <- kinship(gen, method = "astle")

#genotypes
AF.pca <- gen[,adaptive_loci]

pVal <- c()
coeffs <- list()
logLikes <- list()
AICs <- list()
for(i in 1:ncol(AF.pca)){
  # Data frame for GLMM analysis of each SMV
  df <- data.frame(pop=1:nrow(gea_df), env1 = gea_df$env1, env2 = gea_df$env2, SNP=AF.pca[,i])
  
  # Fit GLMM using spaMM
  models <- list()
  models[[1]] <- corrHLfit(SNP ~ env1 + corrMatrix(1|pop), data = df,
                           corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  models[[2]] <- corrHLfit(SNP ~ env2 + corrMatrix(1|pop), data = df,
                           corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  models[[3]] <- corrHLfit(SNP ~ env1 + env2 + corrMatrix(1|pop), data = df,
                           corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")

  
  null.fit <- corrHLfit(SNP ~ 1 + corrMatrix(1|pop), data = df,
                        corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  LLs <- sapply(models, logLik)
  crit.scores <- sapply(models, AIC)
  cond.AICs <- crit.scores[2,]
  best <- which(cond.AICs == min(cond.AICs))
  pVal[i] <- 1-pchisq(2*(LLs[best]-logLik(null.fit)), df=1)
  coeffs[[i]] <- models[[best]]$fixef
  logLikes[[i]] <- LLs
  AICs[[i]] <- cond.AICs
}

results.tab <- matrix(nrow=ncol(AF.pca), ncol=3)
for(z in 1:ncol(AF.pca)){
  results.tab[z, 1] <- pVal[z]
  results.tab[z, 2] <- coeffs[[z]]["env1"]
  results.tab[z, 3] <- coeffs[[z]]["env2"]
}
results.tab <- cbind(results.tab, p.adjust(pVal, method="fdr"))
colnames(results.tab) = c("p", "env1", "env2", "p.adj")
results.tab
