library("here") #paths
#to install LEA uncomment the lines below:
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("LEA")
library("LEA") #LEA
library("vcfR")

#for plotting (from LEA)
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

#kriging
library("automap")
library("raster")

#parallel
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")


#########
#  LEA  #
#########
run_lea_full <- function(gen_filepath, gsd_filepath, loci_filepath){
  
  #Read in data
  gen <- get_gen(gen_filepath)
  gsd_df <- get_gsd(gsd_filepath)
  
  #get adaptive loci
  loci_df <- read.csv(loci_filepath)
  loci_trait1 <- loci_df$trait1 + 1 #add one to convert from python to R indexing
  loci_trait2 <- loci_df$trait2 + 1 #add one to convert from python to R indexing
  adaptive_loci <- c(loci_trait1, loci_trait2)
  neutral_loci <- c(1:nloci)[-adaptive_loci]
  
  
  #remove adaptive loci
  gen <- gen[,-adaptive_loci]
  
  #create gen matrix
  gen <- as.matrix(gen)
  
  #create temporary file with genotypes
  write.geno(gen, here("data","temp_genotypes.geno"))
  
  #Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
  
  #Estimate admixture coefficients using sparse Non-Negative Matrix Factorization algorithms,
  #Code for testing multiple K values:
  obj.snmf <- snmf(here("data","temp_genotypes.geno"), K = 1:20, ploidy = 2, entropy = T, alpha = 100, project = "new")
  plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19, xlim = c(0,30), ylim = c(0.89,1))
  #define K based on "true" value
  #MODIFY LATER TO READ FROM FILE (incorporate into main function?)
  K = 5
  
  #Code for testing one K value
  obj.snmf = snmf(here("data","temp_genotypes.geno"), K = K, ploidy = 2, entropy = T, alpha = 100, project = "new")
  
  #Get Qmatrix
  pred_admix <- Q(obj.snmf, K = K) 
  
  #make grid for kriging
  x.range <- c(0,40)
  y.range <- c(0,40)
  krig_df.grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], len=40),  
                             y=seq(from=y.range[1], to=y.range[2], len=40))
  coordinates(krig_df.grd) <- ~x+y
  gridded(krig_df.grd) <- TRUE
  #print(extent(krig_df.grd)) #FIGURE OUT WHY EXTENT ISNT (0,41,0,41)
  
  pred_krig_admix <- stack()
  for(k in 1:K){
    #krig admix proportions
    krig_df = data.frame(x = gsd_df$x,
                         y = gsd_df$y, 
                         prop = pred_admix[,k])
    
    coordinates(krig_df)=~x+y
    
    krig_res <- autoKrige(prop ~ 1, krig_df, krig_df.grd)
    krig_raster <- raster(krig_res$krige_output)
    print(extent(krig_raster))
    extent(krig_raster) <- c(0,40,0,40)
    pred_krig_admix <- stack(pred_krig_admix, krig_raster)
  }
  
  #convert all values greater than 1 to 1 and all values less than 0 to 0
  pred_krig_admix[pred_krig_admix > 1] <- 1
  pred_krig_admix[pred_krig_admix < 0] <- 0
  
  rbw <- turbo(K)
  par(pty="s", mar=rep(0,4), oma=rep(0,4))
  plot(1, type="n", xlab="", ylab="", axes = NULL, xlim=c(0,40), ylim=c(0,40))
  for(l in 1:K){
    cols <- c(rgb(1,1,1,0), rbw[l])
    kpal <- colorRampPalette(cols, interpolate="linear")
    plot(pred_krig_admix[[l]], zlim=c(0.5,1), col=kpal(100), add=TRUE, legend=FALSE, xlim=c(0,41), ylim =c(0,41))
  }
    
  #remove project
  remove.snmfProject("data/temp_genotypes.snmfProject")

  results <- data.frame(krigcor = mean(krigcor), admixcor = mean(admixcor))
  
  return(results)
  
  
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lea <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library(here)
  
  gen_filepath <- create_filepath(i, "gen")
  
  gsd_filepath <- create_filepath(i, "gsd")
  
  loci_filepath <- create_filepath(i, "loci")
  
  #skip iteration if file does not exist
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    result <- run_lea_full(gen_filepath = gen_filepath,
                      gsd_filepath = gsd_filepath,
                      loci_filepath = loci_filepath)
  }
  
  return(result)
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)
write.csv(lea_out, "LEA_results.csv")