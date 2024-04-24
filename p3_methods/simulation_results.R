# Load required libraries
library(here)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Source a custom function file
source(here("general_functions.R"))

# Print the 'params' object
params

# List files in a specific directory matching a pattern
list.files(here("p1_gnxsims/gnx/LGS_data/stats"), pattern = "OTHER_STATS")

# Define a function to generate file paths based on parameters
make_path <- function(K, phi, m, seed, H, r, it){
  path <- paste0("mod-K", K, "_phi", phi*100, "_m", m*100, "_seed", seed, "_H", H*100, "_r", r*100, "_it-", it,"_spp-spp_0_OTHER_STATS.csv")
  path <- here("p1_gnxsims/gnx/LGS_data/stats", path)
  return(path)
}

# Generate file paths using the 'params' object
files <- pmap_chr(params, make_path)

# Check if the files exist
files_exist <- map_lgl(files, file.exists)

# Create a data frame with file paths, parameters, and existence status
df <- data.frame(file = files, params, files_exist)

# Check if all files exist
all(!df$files_exist)

# Find the unique combinations of parameters for missing files
missing <- 
  df %>% 
  dplyr::filter(files_exist == FALSE) %>%
  distinct(K, phi, m, seed, H, r) %>% 
  head()

# Print summary of missing files
summary(missing)

# Print unique values of parameters for missing files
unique(missing$K)
unique(missing$phi)
unique(missing$m)
unique(missing$seed)
unique(missing$H)
unique(missing$r)

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
  left_join(df, by = "file") %>%
  group_by(t, H, phi) %>%
  summarize(mean_fit = mean(mean_fit, na.rm =TRUE))

# Create a scatter plot using the summarized data
ggplot(pdf2, aes(x = t, y = mean_fit, col = factor(H), pch = factor(phi))) +
  geom_point()