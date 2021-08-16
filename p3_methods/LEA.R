
set.seed(42)

library("here") #paths
library("vcfR")
#to install LEA:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("LEA")
library("LEA") #LEA

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
run_lea_full <- function(gen, gsd_df, loci_df, paramset){
  #get adaptive loci
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
  maxK <- 20
  obj.snmf <- snmf(here("data","temp_genotypes.geno"), K = 1:maxK, ploidy = 2, entropy = T, alpha = 100, project = "new")
  
  #determining best K and picking best replicate for best K (source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
  ce <- list()
  for(k in 1:maxK) ce[[k]] <- cross.entropy(obj.snmf, K=k)
  ce.K <- c()
  for(k in 1:maxK) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-maxK]
  slope <- exp(-diff) - 1
  #K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  
  
  plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)
  abline(v = K, col = "red", lty = "dashed")

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
  
  points(gsd_df$x, gsd_df$y)
  
  #export full krig admix and Qmatrix
  full_krig_admix <- pred_krig_admix
  full_admix <- pred_admix
  writeRaster(full_krig_admix, paste0(paramset,"_krig_admix.tif"), format="GTiff", overwrite = TRUE)
  write.csv(full_admix, paste0(paramset,"_qmat.csv"), row.names = FALSE)
  
  #remove project
  remove.snmfProject("data/temp_genotypes.snmfProject")
  
  results <- data.frame(K = K)
  
  return(results)
  
}

run_lea <- function(gen, gsd_df, loci_df, K, full_krig_admix, full_admix){
  #get adaptive loci
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
  #Code for running one K value
  obj.snmf = snmf(here("data","temp_genotypes.geno"), K = K, ploidy = 2, entropy = T, alpha = 100, project = "new")
  
  #Get Qmatrix
  pred_admix <- Q(obj.snmf, K = K) 
  
  pred_krig_admix <- stack()
  for(k in 1:K){
    #krig admix proportions
    krig_df = data.frame(x = gsd_df$x,
                         y = gsd_df$y, 
                         prop = pred_admix[,k])
    
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
  
  #remove project
  remove.snmfProject("data/temp_genotypes.snmfProject")
  
  
  ###########################
  #  Admix Krig Comparison  #
  ###########################
  #correlation between all layers in stacks
  full_cor <- layerStats(stack(pred_krig_admix, full_krig_admix), "pearson")$'pearson correlation coefficient'
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
  colnames(pred_admix) <- paste0("pred",1:K)
  colnames(full_admix) <- paste0("full",1:K)
  
  rmse <- c()
  for(k in 1:K){
    rmse[k] <- sqrt(mean((pred_admix[,k] - full_admix[,k])^2)) #NEED TO FIGURE OUT HOW TO MAKE SURE K=1 in true is same as K=1 in pred
  }
  
  #correlation between pred and true admix 
  sub_cor <- cor(pred_admix, full_admix)
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
  
  return(results)
  

}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

res_lea <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100)
  
  #save plots
  pdf(paste0("outputs/LEA/LEA_plots_",paramset))
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, "gen")
  gsd_filepath <- create_filepath(i, "gsd")
  loci_filepath <- create_filepath(i, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
    print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LEA
  if(skip_to_next == FALSE){
    gen <- get_data(i, "gen")
    gsd_df <- get_data(i, "gsd")
    loci_df <- get_data(i, "loci")
    
    #subsample full data randomly
    s <- sample(2000, nrow(gsd_df), replace = FALSE)
    gen <- gen[s,]
    gsd_df <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_lea_full(gen, gsd_df, loci_df, paramset)
    result <- data.frame(sampstrat = "full", nsamp = nrow(gsd_df), full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LEA/LEA_results_",paramset,".csv")
    write.csv(data.frame(params[i,], result), csv_file, row.names = FALSE)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        #subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        #run analysis using subsample
        sub_result <- run_lea(subgen, subgsd_df, loci_df, K = full_result$K, full_krig_admix = full_krig_admix, full_admix = full_admix)
        
        #save and format new result
        sub_result <- data.frame(sampstrat = sampstrat, nsamp = nsamp, sub_result)
        
        #export data to csv (temp)
        csv_df <- read.csv(csv_file)
        csv_df <- rbind(csv_df, data.frame(params[i,], sub_result))
        write.csv(csv_df, csv_file, row.names = FALSE)
        
        #bind results
        result <- rbind.data.frame(result, sub_result)
      }
    }
  }
  
  dev.off()
  
  return(result)
  
  gc()
  
}

#stop cluster
stopCluster(cl)

stats_out <- cbind.data.frame(params, res_lfmm)
write.csv(lea_out, "outputs/LEA/LEA_results.csv")
