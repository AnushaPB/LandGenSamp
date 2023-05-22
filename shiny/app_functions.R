get_app_data <- function(method, sampling) {
  df <- read_csv(here("shiny", "app_data", paste0(method, "_", sampling, ".csv")))
  # convert parameters to factors
  df <- var_to_fact(df) 
  return(df)
}