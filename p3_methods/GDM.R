library("here") #paths
library("gdm") #GDM

############
#   Data   #
############
#DEFINE NLOCI
nloci <- 1000

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
s <- sample(1:nrow(gea_df),1000)
loci_df <- loci_df[s,]
gea_df <- gea_df[s,]


###########
#   GDM   #
###########

#CALCULATING PC BASED GENETIC DISTANCE
#convert gen data to matrix
gen <- as.matrix(gea_df[,1:nloci])
#perform PCA
pc <- prcomp(gen)
#Calculate PC distance based on first three PCs (?MODIFY?)
pc_dist <- as.matrix(dist(pc$x[,1:3], diag = TRUE, upper = TRUE))

#Format gdm dataframe
site <- 1:nrow(pc_dist) #vector of sites
gdmGen <- cbind(site, pc_dist) #bind vector of sites with gen distances
gdmPred <- data.frame(site = site, Longitude = gea_df$x, Latitude = gea_df$y, env = gea_df$env)
gdmData <- formatsitepair(gdmGen, bioFormat = 3, predData = gdmPred, XColumn = "Longitude", YColumn = "Latitude", siteCol = "site")
#SCALE DISTANCE FROM 0 to 1 if max(distance) >1 (gdm only works for 0<vals<1) (?MODIFY?)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
gdmData$distance <- range01(gdmData$distance) 

#run GDM
gdm.model <- gdm(gdmData, geo = TRUE)

# Sum coefficients for each predictor (each has 3 splines)
coeffs <- function(gdm.model){
  coefSums <- c()
  for (i in 1:length(gdm.model$predictors)){
    j <- (i * 3) - 2
    coefSums[i] <- sum(gdm.model$coefficients[j:(j+2)])
  }
  
  # Add those values to a simple data frame
  coeffs <- data.frame(predictor = gdm.model$predictors, coefficient = coefSums)
  return(coeffs)
}

predictors <- coeffs(gdm.model)
predictors

#OUTPUT RESULTS
#TBD - need to decide on file structure

