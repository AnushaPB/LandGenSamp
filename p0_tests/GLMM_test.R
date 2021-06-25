## GLMM ##
library(spaMM)
library(statgenGWAS)
#BiocManager::install("LEA")
library("LEA") #LEA

#for plotting (from LEA)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#parallel
library(foreach)
library(doParallel)


#define nloci 
nloci = 10000

#read in geospatial data
file_path = here("data","mod-10k_K5_phi10_m0.25_seed1_H50_r60_it--1_t-300_spp-spp_0.csv")
gsd_df <- read.csv(file_path)
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
head(gsd_df)

#read in genetic data
file_path = here("data","mod-10k_K5_phi10_m0.25_seed1_H50_r60_it--1_t-300_spp-spp_0.vcf")
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


#make kinship matrix (method = c("astle", "IBS", "vanRaden"))
Kmatrix <- kinship(gen, method = "astle")
s <- sample(1:nrow(gea_df),50)
ks <- Kmatrix[s,s]
diag(ks) <- NA
heatmap(ks, keep.dendro=FALSE)

#genotypes
AF.pca <- gen[,c(1,2,3,4,5,adaptive_loci)]


# create class which holds multiple results for each loop iteration.
# Each loop iteration populates three properties: $result1 and $result2 and $result3
# For a great tutorial on S3 classes, see: 
# http://www.cyclismo.org/tutorial/R/s3Classes.html#creating-an-s3-class
multiResultClass <- function(pVal=NULL,coeffs=NULL)
{
  me <- list(
    pVal = pVal,
    coeffs = coeffs
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}


#register cores
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


modres <- foreach(i=1:ncol(AF.pca)) %dopar% {
  library(spaMM)
  result <- multiResultClass()
  
  # Data frame for GLMM analysis of each SMV
  df <- data.frame(pop=1:nrow(AF.pca), env1 = gea_df$env1, env2 = gea_df$env2, SNP=AF.pca[,i])
  
  # Fit GLMM using spaMM
  models <- list()
  models[[1]] <- corrHLfit(SNP ~ env1 + corrMatrix(1|pop), data = df, corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  models[[2]] <- corrHLfit(SNP ~ env2 + corrMatrix(1|pop), data = df, corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  #models[[3]] <- corrHLfit(SNP ~ env1 + env2 + corrMatrix(1|pop), data = df, corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  null.fit <- corrHLfit(SNP ~ 1 + corrMatrix(1|pop), data = df,
                        corrMatrix=Kmatrix, HLmethod="ML", family="gaussian")
  
  LLs <- sapply(models, logLik)
  crit.scores <- sapply(models, AIC)
  cond.AICs <- crit.scores[2,]
  best <- which(cond.AICs == min(cond.AICs))
  result$pVal <- 1-pchisq(2*(LLs[best]-logLik(null.fit)), df=1)
  result$coeffs <- models[[best]]$fixef
  return(result)
}


#stop cluster
stopCluster(cl)

results.tab <- matrix(nrow=ncol(AF.pca), ncol=3)
for(z in 1:ncol(AF.pca)){
  results.tab[z, 1] <- modres[[z]]$pVal
  results.tab[z, 2] <- modres[[z]]$coeffs["env1"]
  results.tab[z, 3] <- modres[[z]]$coeffs["env2"]
}

pVal <- c()
for(i in 1:length(modres)){pVal[i] <- modres[[i]]$pVal}

results.tab <- cbind(results.tab, p.adjust(pVal, method="fdr"))
colnames(results.tab) = c("p", "env1", "env2", "p.adj")
results.tab

par(mfrow=c(1,2))
plot(1:nrow(results.tab), -log10(results.tab[,"p.adj"]))
abline(h=-log10(0.05), col="red", lty=2)

plot(1:nrow(results.tab), results.tab[,"env1"])
points(1:nrow(results.tab), results.tab[,"env2"], col = "black", pch=19)
abline(h=0, col="black", lty=1)

