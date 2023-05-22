
library("here")
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")

#register cores
cores <- 20
cl <- makeCluster(cores)
registerDoParallel(cl)

system.time(
  foreach(i = 1:nrow(params), .packages = c("here", "vcfR", "lfmm", "stringr", "AssocTests", "adegenet", "purrr", "dplyr")) %dopar% {
    
    datadir = here(dirname(getwd()), "p1_gnxsims", "parallel", "LGS_data")
    
    #set of parameter names in filepath form
    paramset <- paste0("K",params[i,"K"],
                       "_phi",params[i,"phi"]*100,
                       "_m",params[i,"m"]*100,
                       "_seed",params[i,"seed"],
                       "_H",params[i,"H"]*100,
                       "_r",params[i,"r"]*100)
    
    gen <- get_data(i, params = params, "gen")
    filepath <- paste0(datadir, "/mod-", paramset,"_it-", params[i,"it"], "_t-1000_spp-spp_0.rda")
    save(gen, file = dilepath)
    
    return(NULL)
    
    gc()
  }
)

#stop cluster
stopCluster(cl)

