library(gdm)
library(tidyverse)

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

# Scale genetic distances from 0 to 1 for GDM
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Data
gsd_df <- read.csv("gdm_problem_example_gsd.csv")
gen <- read.csv("gdm_problem_example_gen.csv")

# Format data for GDM  
gendist <- as.matrix(dist(gen))
site <- 1:nrow(gendist) #vector of sites
gdmGen <- cbind(site, gendist) #bind vector of sites with gen distances

# Test 1: geo = TRUE, no geodist matrix provided
gdmPred <-
  data.frame(
    site = site,
    Longitude = gsd_df$x,
    Latitude = gsd_df$y,
    env1 = gsd_df$env1,
    env2 = gsd_df$env2
  )

gdmData <-
  formatsitepair(
      gdmGen,
      bioFormat = 3,
      predData = gdmPred,
      XColumn = "Longitude",
      YColumn = "Latitude",
      siteColumn = "site"
    )
  
#scale distance from 01
gdmData$distance <- range01(gdmData$distance) 
  
#run GDM
gdm.model1 <- gdm(gdmData, geo = TRUE)
coeffs1 <- coeffs(gdm.model1)

# Test 2: geo = FALSE, geodist matrix provided separately
gdmPred <-
  data.frame(
    site = site,
    Longitude = gsd_df$x,
    Latitude = gsd_df$y,
    env1 = gsd_df$env1,
    env2 = gsd_df$env2
  )


geoDist <- cbind(site, as.matrix(dist(gsd_df[,c("x", "y")], method = "euclidean", diag = TRUE, upper = TRUE)))
gdmData <-
  formatsitepair(
    gdmGen,
    bioFormat = 3,
    predData = gdmPred,
    XColumn = "Longitude",
    YColumn = "Latitude",
    siteColumn = "site", 
    distPreds = list(geo = geoDist)
  )

#scale distance from 01
gdmData$distance <- range01(gdmData$distance) 

gdm.model2 <- gdm(gdmData, geo = FALSE)
coeffs2 <- coeffs(gdm.model2)
coeffs2$predictor <- c("env1", "env2", "Geographic")

# Test 3: geo = FALSE, geodist matrix provided separately and env matrices provided seperately
# model combo
env1Dist <- cbind(site, as.matrix(dist(gsd_df[,"env1"])))
env2Dist <- cbind(site, as.matrix(dist(gsd_df[,"env2"])))
geoDist <- cbind(site, as.matrix(dist(gsd_df[,c("x", "y")])))
gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, REMOVE = rep(1, nrow(gsd_df)))

gdmData <-
  formatsitepair(
    gdmGen,
    bioFormat = 3,
    predData = gdmPred,
    XColumn = "Longitude",
    YColumn = "Latitude",
    siteColumn = "site",
    distPreds = list(env1 = env1Dist, env2 = env2Dist, geo = geoDist)
  )

#remove placeholder column
gdmData <- gdmData[,!grepl("*REMOVE*", colnames(gdmData))]

#scale distance from 01
gdmData$distance <- range01(gdmData$distance) 

#run GDM
gdm.model3 <- gdm(gdmData, geo = FALSE)
coeffs3 <- coeffs(gdm.model3)
coeffs3$predictor <- c("env1", "env2", "Geographic")


# Test 4: geo = FALSE, geodist matrix provided separately and env matrices provided seperately and order flipped
# model combo

env1Dist <- cbind(site, as.matrix(dist(gsd_df[,"env1"], method = "euclidean", diag = TRUE, upper = TRUE)))
env2Dist <- cbind(site, as.matrix(dist(gsd_df[,"env2"], method = "euclidean", diag = TRUE, upper = TRUE)))
geoDist <- cbind(site, as.matrix(dist(gsd_df[,c("x", "y")], method = "euclidean", diag = TRUE, upper = TRUE)))
gdmPred <- data.frame(site = site, Longitude = gsd_df$x, Latitude = gsd_df$y, REMOVE = rep(1, nrow(gsd_df)))

gdmData <-
  formatsitepair(
    gdmGen,
    bioFormat = 3,
    predData = gdmPred,
    XColumn = "Longitude",
    YColumn = "Latitude",
    siteColumn = "site",
    distPreds = list(env2 = env2Dist, env1 = env1Dist, geo = geoDist)
  )

#remove placeholder column
gdmData <- gdmData[,!grepl("*REMOVE*", colnames(gdmData))]

#scale distance from 01
gdmData$distance <- range01(gdmData$distance) 

#run GDM
gdm.model4 <- gdm(gdmData, geo = FALSE)
coeffs4 <- coeffs(gdm.model4)
coeffs4$predictor <- c("env2", "env1", "Geographic")

# Comparison of results
coeffs1$test <- "geo = TRUE, no distance matrices"
coeffs2$test <- "geo as a distance matrix,\nenv as values"
coeffs3$test <- "geo and env as distance matrices"
coeffs4$test <- "geo and env as distance matrices\norder switched"

df <- bind_rows(coeffs1, coeffs2, coeffs3, coeffs4)
ggplot(data = df, aes(y = test, x = coefficient)) +
  geom_col( position = position_dodge(width = 0.9)) +
  facet_wrap(~predictor) +
  ylab("")

vegan::mantel(dist(gsd_df[,"env1"]), dist(gsd_df[,"env2"]))
vegan::mantel(dist(gsd_df[,"env1"]), dist(gsd_df[,c("x", "y")]))
vegan::mantel(dist(gsd_df[,"env2"]), dist(gsd_df[,c("x", "y")]))
cor(gsd_df$env1, gsd_df$env2)
