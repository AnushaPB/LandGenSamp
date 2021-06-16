library("here") #paths
library("conStruct") #conStruct
library("automap") #kriging

############
#   Data   #
############

#DEFINE NLOCI
nloci = 1000
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
s <- sample(1:nrow(gea_df),100)
loci_df <- loci_df[s, neutral_loci]
gea_df <- gea_df[s,]


###############
#  conStruct  #
###############

#define K based on "true" value
#MODIFY LATER
K = 3

#prep data
gen <- as.matrix((gea_df[,1:nloci]))
if(max(gen) != 1){gen <- gen*0.5} #multiply by 0.5 to turn 0/1/2 coding to 0/0.5/1 coding
coords <- as.matrix(gea_df[,c("x","y")])
geoDist <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))

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
