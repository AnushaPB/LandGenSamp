


set.seed(42)

library("here") #paths
library("foreach")
library("doParallel")
library("dplyr")

#read in general functions and objects
source("general_functions.R")

################
#  Cline Test  #
################

prop_cline <- function(gen, loci_df, gsd_df, sig = 0.05){
  res1 <- data.frame(gen[,(loci_df$trait0 + 1)]) %>% purrr::map_dbl(~ cor.test(., gsd_df$env1, method = "kendall")$p.value)
  res2 <- data.frame(gen[,(loci_df$trait1 + 1)]) %>% purrr::map_dbl(~ cor.test(., gsd_df$env1, method = "kendall")$p.value)
  res <- c(res1, res2)
  # note: use sum/length instead of mean because you want NAs to count as the cline not being detected
  return(sum(res < sig, na.rm = TRUE)/length(res))
}

#register cores
cores <- detectCores()
cl <- makeCluster(cores[1]-3) 
registerDoParallel(cl)

system.time(
  res_cline <- foreach(i=1:nrow(params), .combine=rbind) %dopar% {
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    #calculate prop of clines
    pc <- prop_cline(gen, loci_df, gsd_df)
    
    result <- data.frame(params[i,], prop_cline = pc)
    
    return(result)
    
    gc()
  }
)

#stop cluster
stopCluster(cl)

write.csv(res_cline, "outputs/cline_results.csv", row.names = FALSE)
