set.seed(42)

library("here")
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
  res_cline <- foreach(i=1:nrow(params), .combine=rbind, .packages = c("here", "vcfR", "dplyr", "purrr")) %dopar% {
    gen <- get_data(i, params = params, "gen")
    gsd_df <- get_data(i, params = params, "gsd")
    loci_df <- get_data(i, params = params, "loci")
    
    # calculate prop of clines for full dataset
    pc <- prop_cline(gen, loci_df, gsd_df)
    
    result <- data.frame(params[i,], sampstrat = "full", nsamp = nrow(gsd_df), prop_cline = pc)
    
    for(nsamp in npts){
      for(sampstrat in sampstrats){
        # subsample from data based on sampling strategy and number of samples
        subIDs <- get_samples(params[i,], params, sampstrat, nsamp)
        subgen <- gen[subIDs,]
        subgsd_df <- gsd_df[subIDs,]
        
        # calculate prop of clines for sub dataset
        pc <- prop_cline(gen = subgen, loci_df = loci_df, gsd_df = subgsd_df)
        
        # save and format new result
        sub_result <- data.frame(params[i,], sampstrat = sampstrat, nsamp = nsamp, prop_cline = pc)
        
        #bbind results
        result <- bind_rows(result, sub_result)
        
      }
    }
    
    return(result)
    
    gc()
  }
)

#stop cluster
stopCluster(cl)

write.csv(res_cline, "outputs/cline_results.csv", row.names = FALSE)
