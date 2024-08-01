set.seed(42)

library("here") #paths
library("foreach")
library("doParallel")
library("dplyr")
library("tidyr")
#read in general functions and objects
source(here("general_functions.R"))

#######################################
#  Genotype-Environment Correlations  #
#######################################

#register cores
cl <- makeCluster(20)
registerDoParallel(cl)

system.time(
results <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
  library("here")
  library("dplyr")
  library("purrr")

  #set of parameter names in filepath form (for creating temp files)
  paramset <- paste0("K",params[i,"K"],
                     "_phi",params[i,"phi"]*100,
                     "_m",params[i,"m"]*100,
                     "_seed",params[i,"seed"],
                     "_H",params[i,"H"]*100,
                     "_r",params[i,"r"]*100,
                     "_it",params[i,"it"])
  
  
  # Skip iteration if files do not exist
  gsd_filepath <- create_filepath(i, params = params, "gsd")
  skip_to_next <- FALSE
  if(file.exists(gsd_filepath) == FALSE){skip_to_next <- TRUE}
  if(skip_to_next) { print("File does not exist:")
                      print(params[i,]) } 
  if(skip_to_next) { result <- NA } 
  
  if(skip_to_next == FALSE){
    gsd_df <- get_data(i, params = params, "gsd")
    gen <- get_data(i, params = params, "dos")
    
    # Calculate genotype-env correlation
    cor1 <- map(gen, ~cor.test(.x, gsd_df$env1))
    cor2 <- map(gen, ~cor.test(.x, gsd_df$env2))
    p1 <- map_dbl(cor1, "p.value")
    p2 <- map_dbl(cor2, "p.value")
    r1 <- map_dbl(cor1, "estimate")
    r2 <- map_dbl(cor2, "estimate")

    # Calculate correlation between environmental and neutral SNPs
    r <- cor(gen[,1:8], gen[,-c(1:8)])
    # There are NA values but I am dropping them because they represent fixed alleles. I could also count them as 1s theoretically since they are identical.
    rgen <- mean(r, na.rm = TRUE)
    r[is.na(r)] <- 1
    rgen1 <- mean(r)

    result <- 
      data.frame(params[i,], 
      sampstrat = "full", 
      cor1_p = p1,
      cor2_p = p2,
      cor1_r = r1,
      cor2_r = r2,
      rgen = rgen,
      rgen1 = rgen1,
      adaptive = c(rep(TRUE, 8), rep(FALSE, 10000)))

    result <-
      result %>%
      mutate(cor1_sig = cor1_p < 0.05, cor2_sig = cor2_p < 0.05, cor_sig = (cor1_sig + cor2_sig)/2)
  }
  
  return(result)
  
  gc()
}
)

# Stop cluster
stopCluster(cl)

write.csv(results, here("p3_methods", "outputs", "genotype_environment_results.csv"), row.names = FALSE)

