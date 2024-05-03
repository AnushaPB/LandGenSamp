# Load required libraries
library(here)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

# Source a custom function file
source(here("general_functions.R"))

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
  distinct(K, phi, m, seed, H, r, it) 

# There are 439 missing files
nrow(missing)

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

# Combine the list of data frames into a single data frame
# List of variables to group by
vars <- c("H", "m", "phi", "K", "r")

# Use map_dfr to iterate over vars and bind the results together row-wise
results <- 
  purrr::map(vars, ~{
    ggdf <- 
      pdf2 %>%
      group_by_at(c(.x, "t")) %>%
      summarize(
        mean_fit = mean(mean_fit, na.rm =TRUE), 
        sd_fit = mean(mean_fit, na.rm =TRUE), 
        mean_Nt = mean(Nt, na.rm = TRUE), 
        sd_Nt = sd(Nt, na.rm = TRUE),
        .groups = 'drop') %>%
      mutate(ymin = mean_fit - sd_fit, ymax = mean_fit + sd_fit) %>%
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


ggplot(drop_na(results, mean_fit), aes(x = t, y = mean_fit, col = level)) +
  geom_line() +
  labs(x = "Timepoint", y = "Mean fitness", col = "Parameter level") +
  #geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = level, col = NA), alpha = 0.2) +
  facet_wrap(~param, ncol = 1) +
  theme_classic() +
  theme(strip.background = element_blank())

pdf2 %>% filter(t == 0) %>% count()

ggplot(drop_na(results, mean_Nt), aes(x = t, y = mean_Nt, col = level)) +
  geom_line() +
  labs(x = "Timepoint", y = "Mean popiulation size", col = "Parameter level") +
  facet_wrap(~param, ncol = 1) +
  theme_classic() +
  theme(strip.background = element_blank())

results_total <-
  pdf2 %>%
  group_by(t) %>%
  summarize(mean_fit = mean(mean_fit, na.rm =TRUE), sd_fit = mean(mean_fit, na.rm =TRUE), .groups = 'drop') %>%
  drop_na(mean_fit) %>%
  mutate(ymin = mean_fit - sd_fit, ymax = mean_fit + sd_fit)

ggplot(results_total, aes(x = t, y = mean_fit)) +
  geom_line() +
  labs(x = "Timepoint", y = "Mean fitness") +
  #geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  theme_classic() +
  theme(strip.background = element_blank())

results_total <-
  pdf2 %>%
  group_by(t, K) %>%
  summarize(mean_Nt = mean(Nt, na.rm =TRUE), sd_Nt = mean(Nt, na.rm =TRUE), .groups = 'drop') %>%
  mutate(ymin = mean_Nt - sd_Nt, ymax = mean_Nt + sd_Nt) %>%
  mutate(`Population size level` = case_when(K == 1 ~ "Low", K == 2 ~ "High"))

ggplot(results_total, aes(x = t, y = mean_Nt, col = `Population size level`)) +
  geom_line() +
  labs(x = "Timepoint", y = "Mean population size") +
  #geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2) +
  theme_classic() +
  theme(strip.background = element_blank())

