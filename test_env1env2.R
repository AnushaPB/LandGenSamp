a <- read.csv("p3_methods/outputs/lfmm_indsampling_results.csv")
library(tidyverse)
b <- 
  a %>% 
  filter(
    K == 2, 
    phi == 1,
    m == 1,
    H == 0.5,
    r == 0.3,
    sampstrat == "rand", 
    nsamp == 225, 
    H == 0.5, 
    K_method == "tess", 
    lfmm_method == "ridge",
    sig == 0.05,
    padj == "fdr",
    maf == 0.05,
    all == FALSE)

nrow(b)/3/10

names(b)

d <- b %>% pivot_longer(c(emp1_mean, emp2_mean))
ggplot(d) +
  geom_point(aes(x = name, y = value))

                    
a <- read.csv("p3_methods/outputs/mmrr_indsampling_results.csv")
library(tidyverse)
b <- 
  a %>% 
  filter(
    K == 2, 
    phi == 1,
    m == 1,
    H == 0.5,
    r == 0.3,
    sampstrat == "rand", 
    nsamp == 225, 
    H == 0.5)


d <- b %>% pivot_longer(c(env1_coeff, env2_coeff))
ggplot(d) +
  geom_boxplot(aes(x = name, y = value))



ge_df <- read.csv(here("p3_methods", "outputs", "genotype_environment_results.csv"))
df <- 
  ge_df %>%
  # note NA values occur because there is no genetic variance, we are treating this as a cor of 0 and a sig (count of significant cors) of 0
  # only a very small amount are NAs (less than 0.01%)
  mutate(cor1_rNA = case_when(is.na(cor1_r) ~ 0, TRUE ~ cor1_r)) %>%
  mutate(cor2_rNA = case_when(is.na(cor2_r) ~ 0, TRUE ~ cor2_r)) %>%
  mutate(cor_sigNA = case_when(is.na(cor_sig) ~ 0, TRUE ~ cor_sig)) 

d <-
  df %>%
  filter(adaptive == TRUE) %>%
  # add adaptive loci IDs (1-8, loci are in order; repeat 960 times for 960 simulations)
  mutate(id = rep(1:8, 960)) %>%
  pivot_longer(c(cor1_rNA, cor2_rNA)) %>%
  mutate(cor_r = abs(value)) %>%
  # filter to only include correlations for the adaptive genotypes and their corresponding environment
  filter((id %in% c(1:4) & name == "cor1_rNA") | (id %in% c(5:8) & name == "cor2_rNA")) 

ggplot(d) +
  geom_boxplot(aes(x = name, y = cor_r)) +
  facet_wrap(H~m)
