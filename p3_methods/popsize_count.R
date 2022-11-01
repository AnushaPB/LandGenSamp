set.seed(42)

library("here") 
library("foreach")
library("doParallel")

#read in general functions and objects
source("general_functions.R")

#register cores
cores <- 1
cl <- makeCluster(cores)
#not to overload your computer
registerDoParallel(cl)

system.time(
res_popsize <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("stringr")
  
  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])

  #skip iteration if files do not exist
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  skip_to_next <- FALSE
  if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
                      print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  #get pop size
  if(skip_to_next == FALSE){
    gsd_df <- get_data(i, params = params, "gsd")
    result <- data.frame(params[i,], popsize = nrow(gsd_df))
  }
  
  #end pdf()
  # dev.off()
  
  return(result)
  
  gc()
}
)

#stop cluster
stopCluster(cl)

write.csv(res_popsize, "outputs/popsize_results.csv", row.names = FALSE)

# code for analysis
df <- read.csv("outputs/popsize_results.csv")

#remove empty rows (there shouldn't be any, but just in case I do a partial run)
df <- df[complete.cases(df),]

#create result df
tdf <- matrix(ncol = 8)

#ttests
for(i in c("m", "phi", "H", "r", "K")){
  tt <- t.test(df$popsize ~ df[,i], var.equal=TRUE)
  mean1 <- tt$estimate[[1]]
  mean2 <- tt$estimate[[2]]
  t.value  <- tt[1]
  dfs <- tt[2]
  conf.int1 <- tt$conf.int[1]
  conf.int2 <- tt$conf.int[2]
  p.value <- tt[3]
  t.vect <- cbind(param=i, mean1, mean2, t.value, dfs, conf.int1, conf.int2, p.value)
  tdf <- rbind(tdf, t.vect)
}

results_popsize <- tdf[-1,]

write.csv(results_popsize, "popsize_ttest_results.csv")

