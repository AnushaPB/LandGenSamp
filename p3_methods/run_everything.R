
combos <- expand.grid(method = c("lfmm", "rda", "mmrr", "gdm"), sampling = c("individual", "site"))

pmap(combos, ~run_method(method = .x, sampling = .y))
