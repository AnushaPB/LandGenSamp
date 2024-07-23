library(here)

source(here("p3_methods", "general_run_functions.R"))
combos <- expand.grid(method = c("lfmm", "rda", "mmrr", "gdm"), sampling = c("individual", "site"))

pmap(combos, ~run_method(method = .x, sampling = .y))

run_method("mmrr", "site", ncores = 20)
