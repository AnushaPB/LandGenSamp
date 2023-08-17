library(tidyverse)
library(future)
library(future.apply)

# Helper functions

#' Calls out to `free` to get total system memory used
sys_used <- function() {
  .f <- system2("free", "-b", stdout = TRUE)
  as.numeric(unlist(strsplit(.f[2], " +"))[3])
}

#' Write time, and memory usage to log file in CSV format
#' @param .f the file to write to 
#' @param .id identifier for the row to be written
mem_string <- function(.f, .id) {
  .s <- paste(.id, Sys.time(), sys_used(), Sys.getpid(), sep = ",")
  write_lines(.s, .f, append = TRUE)
}

# Inputs
fake_inputs <- 1:16
nsim <- 100
nrows <- 1e6

log_file <- "future_mem_leak_log.csv"
if (fs::file_exists(log_file)) fs::file_delete(log_file)

test_cases <- list(
  list(
    name = "multisession",
    plan = multisession
  ),
  list(
    name = "multicore",
    plan = multicore
  )
)

# Test code

for (.t in test_cases) {
  plan(.t$plan, workers = 4)
  
  # loop over subsets of the data 
  final_out <- future_lapply(fake_inputs, function(.i) {
    # loop over simulations
    out <- future_lapply(1:nsim, function(.j) {
      # in real life this would be doing simulations, 
      # but here we just create "results" using rnorm()
      res <- data.frame(
        id = rep(.j, nrows),
        col1 = rnorm(nrows) * .i,
        col2 = rnorm(nrows) * .i,
        col3 = rnorm(nrows) * .i,
        col4 = rnorm(nrows) * .i,
        col5 = rnorm(nrows) * .i,
        col6 = rnorm(nrows) * .i
      )
      
      # write memory usage to file
      mem_string(log_file, .t$name)
      
      # in real life I would write res to file to read in later, but here we
      # only return head of df so we know the returned value isn't filling up memory
      res %>% slice_head(n = 10) 
    })
  })
  
  # clean up any leftover objects before testing the next plan
  try(rm(final_out))
  try(rm(out))
  try(rm(res))
}