library(here)
library(tidyverse)
source(here("general_functions.R"))

# Test of many simulations for one iteration
get_stats <- function(K, phi, m, H, r){
  wdir <- here("p1_gnxsims", "gnx")
  folder_name <- paste0("ttestall_may2/GNX_mod-ttestall_K", K, "_phi", phi*100, "_m", m*100, "_seed1_H", H*100, "_r", r*100)
  file_name <- paste0("mod-ttestall_K", K, "_phi", phi*100, "_m", m*100, "_seed1_H", H*100, "_r", r*100, "_it-0_spp-spp_0_OTHER_STATS.csv")
  path <- here(wdir,  paste0(folder_name, "/it-0/spp-spp_0/", file_name))
  if (!file.exists(path)) {warning(paste0("File does not exist: ", path)); return(NULL)}
  df <- 
    read_csv(path) %>% 
    mutate(K = K, phi = phi, m = m, H = H, r = r)
  return(df)
}

combos <- expand.grid(K = c(1, 2), phi = c(0.5, 1), m = c(0.25, 1), H = c(0.05, 0.5), r = c(0.3, 0.6))

ls <- 
  pmap(combos, get_stats) %>% 
  compact()
  
df <- 
  ls %>%
  bind_rows() %>%
  mutate(group = paste0("K", K, "_phi", phi, "_m", m, "_H", H, "_r", r))

length(unique(df$group))/32*100

mean_df <-
  df %>%
  group_by(t) %>%
  drop_na(mean_fit) %>%
  summarise(mean_fit = mean(mean_fit, na.rm = TRUE), .groups = "drop")

# test for stationarity
adf_test <- function(subdf, nsteps = NULL, var = "mean_fit", timepoints = NULL){
  # Creating a time series object
  ts <- 
    subdf %>% 
    pull(.data[[var]]) %>%
    ts()

  # Testing for stationarity of different time points
  # Note: used na.omit() instead of drop_na(mean_fit) because then the time series would be shorter than the original data
  # and the time points would not line up correctly
  if (is.null(timepoints)) timepoints <- map(seq(0, 3500, nsteps), function(x) x:(x+nsteps))
  p_value <- map_dbl(timepoints, ~tseries::adf.test(na.omit(ts[.x]), alternative = "stationary")$p.value)
  p_df <- data.frame(t_start = map_dbl(timepoints, 1), t_end = map_dbl(timepoints, ~.x[length(.x)]), p = p_value)
  
  ts_df <- 
    p_df %>%
    bind_cols(unique(select(subdf, K, phi, m, H, r))) %>%
    mutate(nsteps = as.character(nsteps))

  return(ts_df)
}

# Test for stationarity in mean fitness
st_df <-
  map(ls, adf_test, nsteps = 1000) %>% 
  bind_rows() %>%
  mutate(nsteps = as.character(nsteps)) %>%
  mutate(stationarity = p < 0.05) 

count_st <-
  st_df %>%
  group_by(t_start, t_end, nsteps) %>%
  summarize(stationarity = sum(stationarity)/n() * 100) %>%
  drop_na()

ggplot() +
  geom_rect(data = count_st, aes(xmin = t_start, xmax = t_end, ymin = -Inf, ymax = Inf, fill = stationarity), alpha = 0.9) +
  geom_vline(data = count_st, aes(xintercept = t_end), lty = "dashed", alpha = 0.5, col = "white") +
  geom_line(data = drop_na(df, mean_fit), aes(x = t, y = mean_fit, group = group), alpha = 0.3) +
  geom_line(data = mean_df, aes(x = t, y = mean_fit), lwd = 1) +
  labs(x = "Timepoint", y = "Mean fitness", fill = "% at\nstationarity") +
  scale_fill_viridis_c(option = "mako", end = 1, begin = 0, direction = -1, limits = c(15, 100)) +
  theme_classic() +
  theme(strip.background = element_blank())

# Test with population size
st_df <-
  bind_rows(map(ls, adf_test, nsteps = 1000, var = "Nt")) %>% 
  mutate(nsteps = as.character(nsteps)) %>%
  mutate(stationarity = p < 0.05) 

count_st <-
  st_df %>%
  group_by(t_start, t_end, nsteps) %>%
  summarize(stationarity = sum(stationarity)/n() * 100) %>%
  drop_na()

ggplot() +
  geom_rect(data = count_st, aes(xmin = t_start, xmax = t_end, ymin = -Inf, ymax = Inf, fill = stationarity), alpha = 0.9) +
  geom_vline(data = count_st, aes(xintercept = t_end), lty = "dashed", alpha = 0.5, col = "white") +
  geom_line(data = drop_na(df, Nt), aes(x = t, y = Nt, group = group), alpha = 0.3) +
  labs(x = "Timepoint", y = "Nt", fill = "% at\nstationarity") +
  facet_wrap(~nsteps, ncol = 1) +
  scale_fill_viridis_c(option = "mako", end = 1, begin = 0.1, direction = -1, limits = c(15, 100)) +
  theme_classic() +
  theme(strip.background = element_blank())

# JP test

library(nlraa)
library(minpack.lm)

# test for https://gradcylinder.org/post/linear-plateau/
jp_test <- function(subdf, nsteps = NULL, var = "mean_fit", timepoints = NULL){

  subdf <- subdf %>% drop_na(mean_fit)
  fit <- nlsLM(formula = mean_fit ~ SSlinp(t, a, b, jp), data = subdf)
  jp <- 
    summary(fit)$coefficients %>% 
    as.data.frame() %>%
    mutate(name = rownames(.)) %>%
    filter(name == "jp") %>%
    pull(Estimate) 
  
  jp_df <- 
    data.frame(jp = jp, pred = predict(fit), obs = subdf$mean_fit, resids = residuals(fit), t = subdf$t) %>%
    bind_cols(unique(select(subdf, K, phi, m, H, r))) 

  return(jp_df)
}

jp_df <-
  map(ls, jp_test) %>%
  bind_rows() %>%
  mutate(group = paste0("K", K, "_phi", phi, "_m", m, "_H", H, "_r", r))

mean_jp <-
  jp_df %>%
  summarize(mean = mean(jp, na.rm = TRUE))

ggplot() +
  geom_vline(data = jp_df, aes(xintercept = jp), lty = "dashed", alpha = 0.5, col = "red") +
  geom_line(data =jp_df, aes(x = t, y = obs, group = group), alpha = 0.3) +
  #geom_line(data =jp_df, aes(x = t, y = pred, group = group), alpha = 0.3, col = "red") +
  geom_vline(data = mean_jp, aes(xintercept = mean), lwd = 1.5, col = "red") +
  labs(x = "Timepoint", y = "Mean fitness", fill = "% at\nstationarity") +
  theme_classic() +
  theme(strip.background = element_blank())

ggplot() +
  geom_vline(data = jp_df, aes(xintercept = jp), lty = "dashed", alpha = 0.5, col = "red") +
  geom_point(data =jp_df, aes(x = t, y = resids, group = group), alpha = 0.1, cex = 0.3) +
  geom_vline(data = mean_jp, aes(xintercept = mean), lwd = 1.5, col = "red") +
  labs(x = "Timepoint", y = "Model residuals", fill = "% at\nstationarity") +
  theme_classic() +
  theme(strip.background = element_blank())


# Pi test
get_pi <- function(){
  pi_df <- 
    map(
      list.files(here("p1_gnxsims", "gnx", "pi")), 
      ~read_csv(here("p1_gnxsims", "gnx", "pi", .x)) %>% 
      mutate(file = .x, t = as.numeric(t), pi = as.numeric(pi))
      ) %>%
    bind_rows() %>%
    separate(
      file, 
      into = c("ttestall", "K", "phi", "m", "seed", "H", "r", "it", "csv"), 
      sep = "_", 
      extra = "merge", 
      fill = "right"
    ) %>%
    mutate(
      K = parse_number(K),
      phi = parse_number(phi) / 100,
      m = parse_number(m) / 100,
      seed = parse_number(seed),
      H = parse_number(H) / 100,
      r = parse_number(r) / 100,
      it = parse_number(it)
    ) %>%
    select(-ttestall, -csv) %>%
    mutate(group = paste0("K", K, "_phi", phi, "_m", m, "_H", H, "_r", r, "_seed", seed, "_it", it)) 

  return(pi_df)
}

pi_df <- get_pi()

ggplot(pi_df) + 
  geom_line(aes(x = t, y = pi, group = group)) +
  theme_classic() +
  labs(x = "Timepoint", y =  expression(pi))
