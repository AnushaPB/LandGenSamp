library("here") #paths
library("conStruct") #conStruct
library("automap") #kriging

#parallel
library("parallel")
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")



############
#   Data   #
############

#define nloci 
nloci = 10000

#read in geospatial data
file_path = here("data","mod-10k_K2_phi50_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.csv")
gsd_df <- read.csv(file_path)
#Extract env values from lists
gsd_df$env1 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=,)')) 
gsd_df$env2 <- as.numeric(stringr::str_extract(gsd_df$e, '(?<=, )[^,]+(?=\\])')) 
#Extract phenotype values from lists
gsd_df$z1 <- as.numeric(stringr::str_extract(gsd_df$z, '(?<=\\[)[^,]+(?=,)'))
gsd_df$z2 <- as.numeric(stringr::str_extract(gsd_df$z, '(?<=, )[^,]+(?=\\])')) 
head(gsd_df)

#read in genetic data
file_path = here("data","mod-10k_K2_phi50_m25_seed1_H50_r60_it--1_t-500_spp-spp_0.vcf")
vcf <- read.vcfR(file_path)
#convert to matrix
x <- vcfR2genlight(vcf) 
gen <- as.matrix(x)

#FOR TESTING USE SUBSAMPLE(!!!COMMENT OUT!!!)
set.seed(42)
s <- sample(1:nrow(gea_df),1000)
gsd_df <- gsd_df[s,]
gen <- gen[s,]


###############
#  conStruct  #
###############

#prep data
gen <- as.matrix(gen)
if(max(gen) != 1){gen <- gen*0.5} #multiply by 0.5 to turn 0/1/2 coding to 0/0.5/1 coding
coords <- as.matrix(gsd_df[,c("x","y")])
geoDist <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))

#cross validation to determine K
cl <- makeCluster(4)
registerDoParallel(cl)

my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 8,
                         K = 1:4,
                         freqs = gen,
                         data.partitions = NULL,
                         geoDist = geoDist,
                         coords = coords,
                         n.iter = 1e3,
                         make.figs = TRUE,
                         save.files = FALSE,
                         parallel = FALSE,
                         n.nodes = 4,
                         prefix="temp")

stopCluster(cl)

#define K based on "true" value
#MODIFY LATER
K = 4

#run construct
results <- conStruct(spatial = TRUE,
                    K = K,
                    freqs = gen,
                    geoDist = geoDist,
                    coords = coords,
                    prefix = "temp")

#get admix proportions
pred_admix <- results$chain_1$MAP$admix.proportions
make.structure.plot(admix.proportions = pred_admix)


pred_krig_admix <- stack()
for(i in 1:K){
  #krig admix proportions
  krig_df = data.frame(x = gsd_df$x,
                       y = gsd_df$y, 
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
