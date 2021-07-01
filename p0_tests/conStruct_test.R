library("here") #paths
library("conStruct") #conStruct
library("automap") #kriging
library("parallel")
library("foreach")
library("doParallel")
library("vcfR")

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
s <- sample(1:nrow(gsd_df),100)
gsd_df <- gsd_df[s,]
gen <- gen[s,]

#adaptive loci
loci_trait1 <- c(1731,4684,4742,6252) + 1 #add one to convert from python to R indexing
loci_trait2 <- c(141,1512,8481,9511) + 1 #add one to convert from python to R indexing
adaptive_loci <- c(loci_trait1, loci_trait2)
neutral_loci <- c(1:nloci)[-adaptive_loci]

###############
#  conStruct  #
###############

#prep data
#prep data
gen <- as.matrix(gen[,])
if(max(gen) != 1){gen <- gen*0.5} #multiply by 0.5 to turn 0/1/2 coding to 0/0.5/1 coding
coords <- as.matrix(gsd_df[,c("x","y")])
geoDist <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))

#cross validation to determine K
#register cores
cores=detectCores()
cl <- makeCluster(cores[1]-2,type="PSOCK") #not to overload your computer
registerDoParallel(cl)

my.xvals <- x.validation(train.prop = 0.9,
                         n.reps = 8,
                         K = 1:3,
                         freqs = gen[,sample(neutral_loci,1000)],
                         data.partitions = NULL,
                         geoDist = geoDist,
                         coords = coords,
                         prefix = "example",
                         n.iter = 100, #change back to 1e3
                         make.figs = TRUE,
                         save.files = FALSE,
                         parallel = TRUE,
                         n.nodes = cores[1]-2)

stopCluster(cl)

#define K based on "true" value
#MODIFY LATER
K = 2


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


pred_krig_admix <- raster::stack()
for(i in 1:K){
  #krig admix proportions
  krig_df = data.frame(x = gsd_df$x,
                       y = gsd_df$y, 
                       prop = pred_admix[,i])
  
  coordinates(krig_df)=~x+y
  
  x.range <- c(0,40)
  y.range <- c(0,40)
  
  krig_df.grd <- expand.grid(x=seq(from=0, to=40, len=40),  
                             y=seq(from=0, to=40, len=40))
  coordinates(krig_df.grd) <- ~x+y
  gridded(krig_df.grd) <- TRUE
  #extent(krig_df.grid) #FIGURE OUT WHY EXTENT ISNT (0,40,0,40)
  
  krig_res <- autoKrige(prop ~ 1, krig_df, krig_df.grd)
  pred_krig_admix <- stack(pred_krig_admix, raster(krig_res$krige_output))
}

plot(pred_krig_admix)
plot(env1)
plot(env2)


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
