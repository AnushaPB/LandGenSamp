set.seed(42)
#paths
library("here") 
#to install LFMM:
#devtools::install_github("bcm-uga/lfmm")
library("lfmm") #LFMM
library("vcfR")
library("foreach")
library("doParallel")
library("AssocTests")
library("purrr")

#read in general functions and objects
source("general_functions.R")

##########
#  LFMM  #
##########


#register cores
cores <- 25
cl <- makeCluster(cores)
#not to overload your computer
registerDoParallel(cl)

system.time(
res_lfmm <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("vcfR")
  library("lfmm")
  library("stringr")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  #create pdf to store plots
  # pdf(paste0("outputs/LFMM/plots/lfmm_plots_",paramset,".pdf"))
  
  #skip iteration if files do not exist
  gen_filepath <- create_filepath(i, params = params, "gen")
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  loci_filepath <- create_filepath(i, params = params, "loci")
  skip_to_next <- FALSE
  if(file.exists(loci_filepath) == FALSE | file.exists(gen_filepath) == FALSE | file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:"); print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #run LFMM
  if(skip_to_next == FALSE){
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    #subsample full data randomly
    s <- sample(nrow(gsd_df), 2000, replace = FALSE)
    gen_2k <- gen[s,]
    gsd_df_2k <- gsd_df[s,]
    
    #run model on full data set
    full_result <- run_lfmm(gen_2k, gsd_df_2k, loci_df, K = NULL)
    result <- data.frame(params[i,], sampstrat = "full", nsamp = 2000, full_result)
    
    #write full datafile (temp)
    csv_file <- paste0("outputs/LFMM/LFMM_results_",paramset,".csv")
    write.csv(result, csv_file, row.names = FALSE)
    
    #run subset data
    sampcombos <- expand.grid(sampstrats, npts)
    sub_result <- map_dfr(sampcombos, run_sub, i, params, gen, gsd_df)
    
    #bind results
    result <- rbind.data.frame(result, sub_result)
    
    #export data to csv (temp)
    csv_df <- read.csv(csv_file)
    csv_df <- rbind(csv_df, sub_result)
    write.csv(csv_df, csv_file, row.names = FALSE)
  }
  
  #end pdf()
  # dev.off()
  
  return(result)
  
  gc()
}
)


#stop cluster
stopCluster(cl)

write.csv(res_lfmm, "outputs/LFMM/lfmm_results.csv", row.names = FALSE)

