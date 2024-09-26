library(here)
library(tidyverse)
library(furrr)
source(here("general_functions.R"))

future::plan(future::multisession, workers = 20)
furrr::future_walk(
  1:nrow(params),
  \(i) {

    # Create the file path and check if it exists
    file_path <- create_filepath(i, params, type = "gen")
    if (!file.exists(file_path)) warning(paste("does not exist:", file_path))

    # Get the VCF
    gen <- get_data(i, params = params, "gen")

    # Create a new file path
    file_name <- basename(file_path)
    new_file_path <- gsub("mod-(.*?)_", "dos-\\1_", file_name)
    new_file_path <- gsub(".vcf", ".csv", new_file_path)
    new_file_path <- here::here("p1_gnxsims", "gnx", "dosage", new_file_path)

    write.csv(gen, new_file_path, row.names = TRUE)
  }, 
  .options = furrr_options(seed = TRUE, packages = get_packages()), 
  .progress = TRUE
)
future::plan("sequential")

# check files exist
result <- purrr::map(
  1:nrow(params),
  \(i) {
    file_path <- create_filepath(i, params, type = "dos")
    return(data.frame(file_path = file_path, exists = file.exists(file_path)))
  }
)
df <- result %>% bind_rows() %>% filter(!exists)
print(df)
