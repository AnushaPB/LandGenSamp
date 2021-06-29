library("here") #paths
library("conStruct") #conStruct
library("automap") #kriging

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
library(parallel)
library(foreach)
library(doParallel)



#prep data
gen <- as.matrix((gea_df[,1:nloci]))
if(max(gen) != 1){gen <- gen*0.5} #multiply by 0.5 to turn 0/1/2 coding to 0/0.5/1 coding
coords <- as.matrix(gea_df[,c("x","y")])
geoDist <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))

#cross validation to determine K
cl <- makeCluster(4, type="FORK")
registerDoParallel(cl)

my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 8,
                         K = 1:3,
                         freqs = gen,
                         data.partitions = NULL,
                         geoDist = geoDist,
                         coords = coords,
                         n.iter = 1e3,
                         make.figs = TRUE,
                         save.files = FALSE,
                         parallel = TRUE,
                         n.nodes = 4)

stopCluster(cl)

#define K based on "true" value
#MODIFY LATER
K = 3


#run construct
results <- conStruct(spatial = TRUE,
                    K = K,
                    freqs = gen,
                    geoDist = geoDist,
                    coords = coords,
                    prefix = "spK2")

#get admix proportions
pred_admix <- results$chain_1$MAP$admix.proportions
make.structure.plot(admix.proportions = pred_admix)


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


#ADD OUTPUT
