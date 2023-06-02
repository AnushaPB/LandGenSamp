set.seed(42)

purrr::map(get_packages(), library, character.only = TRUE)

#read in general functions and objects
source(here("general_functions.R"))
source(here("p4_analysis", "general_run_functions.R"))

fst_calc <- function(i, params){
  
  filepath <- create_filepath(i, params = params, "dos")
  if (!file.exists(filepath)){ print("File does not exist:"); print(params[i,]); return(NA)}
  
  # Get data
  gen <- get_data(i, params = params, "dos")
  gsd_df <- get_data(i, params = params, "gsd")
  
  # Create a blank raster with 100 rows and 100 columns
  r <- raster(nrow = 10, ncol = 10)
  y <- init(r, fun=1:ncell(r))
  extent(y) <- extent(0, 100, -100, 0)

  set.seed(245)
  s <- sample(nrow(gsd_df), 2000)
  pop <- extract(y, gsd_df[s, c("x", "y")])
  fst <- fst.dosage(gen[s,], pop = pop)
  return(fst["All"])
}

future::plan(future::multisession, workers = 1)
fst <- future_map_dbl(1:nrow(params), ~fst_calc(.x, params),
                      .progress = TRUE, 
                      .options = furrr_options(seed = TRUE, packages = get_packages()))

future::plan("sequential")
fst_df <- data.frame(params, fst = fst)

write.csv(fst_df, here("p3_methods/outputs/fst_results.csv"), row.names = FALSE)


###
library(purrr)
library(dplyr)

# Specify the column names for which you want to run linear models
columns <- grep("_TPR", names(mmrr_df), value = TRUE)

# Function to run linear models and extract coefficients and p-values
run_lm <- function(column_name, df) {
  model <- lm(as.formula(paste(column_name, "~ fst")), data = df)
  coef <- coef(summary(model))[, "Estimate"]["fst"]
  pval <- coef(summary(model))[, "Pr(>|t|)"]["fst"]
  tibble(column_name = column_name, coefficient = coef, pvalue = pval) %>%  
    mutate(across(c(coefficient, pvalue), round, 2))
}

# Apply the function to the specified columns using purrr's map_df function
result_df <- map_df(columns, ~run_lm(.x, df = mmrr_df))

# Print the resulting data frame
print(result_df)

# Specify the column names for which you want to run linear models
columns <- c("TPRCOMBO", "FDRCOMBO")

# Apply the function to the specified columns using purrr's map_df function
result_df <- map_df(columns, ~run_lm(.x, df = lfmm_df))

# Print the resulting data frame
print(result_df)

