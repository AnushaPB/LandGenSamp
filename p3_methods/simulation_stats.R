# Load required libraries
library(here)
library(purrr)
library(tidyr)
library(dplyr)
library(readr)

# Source a custom function file
source(here("general_functions.R"))

# Define a function to generate file paths based on parameters
make_path <- function(K, phi, m, seed, H, r, it){
  mod_name <- paste0("mod-K", K, "_phi", phi*100, "_m", m*100, "_seed", seed, "_H", H*100, "_r", r*100)
  folder_name <- here("p1_gnxsims", "gnx", "LGS_data", paste0("GNX_", mod_name), paste0("it-", it), "spp-spp_0")
  path <- here(folder_name, paste0(mod_name, "_it-", it, "_spp-spp_0_OTHER_STATS.csv"))
  return(path)
}

# Generate file paths using the 'params' object
files <- pmap_chr(params, make_path)

# Check if the files exist
files_exist <- map_lgl(files, file.exists)

# Create a data frame with file paths, parameters, and existence status
df <- data.frame(file = files, params, files_exist)

# Check if all files exist
stopifnot(all(df$files_exist))

# Define a function to read a CSV file and add the file path as a column
read <- function(x){
  df <- read.csv(x)
  df$file = x
  return(df)
}

# Read the existing files into a list of data frames
pdf <- 
 df %>%
 dplyr::filter(files_exist == TRUE) %>%
 dplyr::pull(file) %>%
 map(read)

# Combine the list of data frames into a single data frame
pdf2 <- 
  bind_rows(pdf) %>%
  left_join(df, by = "file")

# Use map_dfr to iterate over vars and bind the results together row-wise to get low/high summaries
stats <- 
  purrr::map(c("H", "m", "phi", "K", "r"), ~{
    ggdf <- 
      pdf2 %>%
      group_by_at(c(.x, "t")) %>%
      summarize(
        sd_fit = sd(mean_fit, na.rm =TRUE), 
        mean_fit = mean(mean_fit, na.rm =TRUE), 
        sd_Nt = sd(Nt, na.rm = TRUE),
        mean_Nt = mean(Nt, na.rm = TRUE), 
        min_Nt = min(Nt, na.rm = TRUE),
        max_Nt = max(Nt, na.rm = TRUE),
        .groups = 'drop') %>%
      mutate(
        ymin_fit = mean_fit - sd_fit, ymax_fit = mean_fit + sd_fit,
        ymin_Nt = mean_Nt - sd_Nt, ymax_Nt = mean_Nt + sd_Nt
        ) %>%
      mutate(param = .x) 

    ggdf$level <- factor(case_when(ggdf[[.x]] == min(ggdf[[.x]]) ~ "Low", ggdf[[.x]] == max(ggdf[[.x]]) ~ "High"), levels = c("Low", "High"))
    ggdf <- ggdf %>% select(-all_of(.x))
  }) %>%
  bind_rows() %>%
  mutate(old_param = param, param = case_when(
    param == "H" ~ "Spatial autocorrelation",
    param == "m" ~ "Migration",
    param == "phi" ~ "Selection strength",
    param == "K" ~ "Population size",
    param == "r" ~ "Environmental correlation"
  ))

write_csv(stats, here("p3_methods", "outputs", "simulation_stats.csv"))
