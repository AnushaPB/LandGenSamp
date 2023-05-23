set.seed(42)

#parallel
library(furrr)

#read in general functions and objects
source("general_functions.R")
source("gdm_functions.R")

#register cores
future::plan(future::multisession, workers = 2)

res_gdm <- furrr::future_map(1:nrow(params), \(i) {
  #skip iteration if files do not exist
  skip_to_next <- skip_check(i, params)
  if(skip_to_next) { result <- NA } 
  
  #run GDM
  tryCatch({
    if(skip_to_next == FALSE){
      gen <- get_data(i, params = params, "gen")
      gsd_df <- get_data(i, params = params, "gsd")
      
      #subsample full data randomly
      s <- sample(nrow(gsd_df), 1000, replace = FALSE)
      gen_2k <- gen[s,]
      gsd_df_2k <- gsd_df[s,]
      
      #run model on full data set
      full_result <- run_gdm(gen_2k, gsd_df_2k, distmeasure = "euc")
      fullratio <- (abs(full_result$env1_coeff) + abs(full_result$env2_coeff))/abs(full_result$geo_coeff)
      result <- data.frame(params[i,], 
                           sampstrat = "full", 
                           nsamp = 2000, 
                           full_result, 
                           ratio = fullratio,
                           env1_err = NA, 
                           env2_err = NA, 
                           geo_err = NA,
                           ratio_err = NA)
      result <- sapply(result, as.character)
      
      for(nsamp in npts){
        for(sampstrat in sampstrats){
          #subsample from data based on sampling strategy and number of samples
          subIDs <- get_samples(params[i,], sampstrat, nsamp)
          subgen <- gen[subIDs,]
          subgsd_df <- gsd_df[subIDs,]
          
          #run analysis using subsample
          sub_result <- run_gdm(subgen, subgsd_df, distmeasure = "euc")
          
          #calculate err if not null
          if(sub_result$env1_coeff == "NULL" | sub_result$env2_coeff == "NULL" | sub_result$geo_coeff == "NULL"){
            subratio <- "NULL"
            env1_err <- "NULL"
            env2_err <- "NULL"
            geo_err <- "NULL"
            ratio_err <- "NULL"
          } else {
            subratio <- (abs(sub_result$env1_coeff) + abs(sub_result$env2_coeff))/abs(sub_result$geo_coeff)
            env1_err <- err_coeff(full_result$env1_coeff, sub_result$env1_coeff)
            env2_err <- err_coeff(full_result$env2_coeff, sub_result$env2_coeff)
            geo_err <- err_coeff(full_result$geo_coeff, sub_result$geo_coeff)
            ratio_err <- err_coeff(fullratio, subratio)
          }
          
          
          #save and format new result
          sub_result <- data.frame(params[i,], 
                                   sampstrat = sampstrat, 
                                   nsamp = nsamp, 
                                   sub_result, 
                                   ratio = subratio,
                                   env1_err = env1_err, 
                                   env2_err = env2_err, 
                                   geo_err = geo_err, 
                                   ratio_err = ratio_err)
          
          #bind results 
          sub_result <- sapply(sub_result, as.character)
          result <- bind_rows(result, sub_result)
        }
      }
    }
  },
  error = function(e) {
    err <<- conditionMessage(e)
    write.table(err, "error_msg.txt")
    write.csv(gen_2k, "error_gen_2k.csv", row.names = FALSE)
    write.csv(gsd_2k_df, "error_gsd_df_2k.csv", row.names = FALSE)
    write.csv(subgen, "error_subgen.csv", row.names = FALSE)
    write.csv(subgsd_df, "error_subgsd_df.csv", row.names = FALSE)
    write.csv(gen, "error_gen.csv", row.names = FALSE)
    write.csv(gsd_df, "error_gsd_df.csv", row.names = FALSE)
    
    message(err)
    
    stop(err)
  }
  )
  
  return(result)
}, 
.options = furrr::furrr_options(seed = TRUE, packages = c("vcfR", "gdm", "adegenet", "stringr", "dplyr", "here", "purrr")),
.progress = TRUE
)

res_gdm <- res_gdm %>% bind_rows()

write.csv(res_gdm, "outputs/gdm_indsampling_results.csv", row.names = FALSE)

