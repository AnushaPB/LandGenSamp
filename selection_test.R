library(here)
library(dplyr)
library(vcfR)
library(ggplot2)
source(here("general_functions.R"))
phi100_t1000 <- get_data(30, params, "dos")
phi100_t1000df <- get_data(30, params, "gsd")
phi100_t0 <- get_gen("p1_gnxsims/gnx/LGS_data/t0_data/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.vcf")
phi100_t0df <- get_gsd("p1_gnxsims/gnx/LGS_data/t0_data/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.csv")

ggplot(phi100)

phi100_t1000 <- get_data(32, params, "dos")
phi100_t1000df <- get_data(32, params, "gsd")
phi100_t0 <- get_gen("~/Anusha/GitHub/LandGenSamp/p1_gnxsims/parallel/LGS_data/t0_data/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.vcf")
phi100_t0df <- get_gen("~/Anusha/GitHub/LandGenSamp/p1_gnxsims/parallel/LGS_data/t0_data/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.csv")

s1000 <- which(phi100_t1000df$z1 > 0.9)
s0 <- which(phi100_t0df$z1 > 0.9)
p1000 <- colMeans(phi100_t1000[s1000,]/2)
p0 <- colMeans(phi100_t0[s0,]/2)

p1000[1:4] - p0[1:4]


s1000 <- which(phi100_t1000df$z2 > 0.5)
s0 <- which(phi100_t0df$z2 > 0.5)
p1000 <- colMeans(phi100_t1000[s1000,]/2)
p0 <- colMeans(phi100_t0[s0,]/2)

p1000[1:4] - p0[1:4]

df <- data.frame(phi100_t0df, gen0 = phi100_t0[,1], t = 0)

df2 <- data.frame(phi100_t1000df, gen0 = phi100_t1000[,1], t = 1000)

df <- bind_rows(df[1:100,], df2[1:100,])

ggplot(data = df[1:100,], aes(x = x, y = y, col = factor(gen0))) + 
  geom_point() +
  facet_wrap(. ~ t) 


s1000 <- which(phi100_t1000df$x > 0.0 & phi100_t1000df$x < 10.0 & phi100_t1000df$y < 0 & phi100_t1000df$y > -10)
s0 <- which(phi100_t0df$x > 0.0 & phi100_t0df$x < 10.0 & phi100_t0df$y < 0 & phi100_t0df$y > -10)

p1000 <- colMeans(phi100_t1000[s1000,]/2)[1]
p0 <- colMeans(phi100_t0[s0,]/2)[1]
env0 <- mean(phi100_t0df$env1)
env1000 <- mean(phi100_t1000df$env1)

AA1 <- mean(phi100_t0[s0,1] == 2)
Aa1 <- mean(phi100_t0[s0,1] == 1)
aa1 <- mean(phi100_t0[s0,1] == 0)

AA2 <- mean(phi100_t1000[s1000,1] == 2)
Aa2 <- mean(phi100_t1000[s1000,1] == 1)
aa2 <- mean(phi100_t1000[s1000,1] == 0)

sAa = 1 - (Aa2/Aa1)
saa = 1 - (aa2/aa1)


s1000 <- which(phi100_t1000df$x > 0.0 & phi100_t1000df$x < 10.0 & phi100_t1000df$y < 0 & phi100_t1000df$y > -10)
s0 <- which(phi100_t0df$x > 0.0 & phi100_t0df$x < 10.0 & phi100_t0df$y < 0 & phi100_t0df$y > -10)

p1000 <- colMeans(phi100_t1000[s1000,]/2)[1]
p0 <- colMeans(phi100_t0[s0,]/2)[1]
env0 <- mean(phi100_t0df$env1)
env1000 <- mean(phi100_t1000df$env1)

AA1 <- mean(phi100_t0[s0,2] == 2)
Aa1 <- mean(phi100_t0[s0,2] == 1)
aa1 <- mean(phi100_t0[s0,2] == 0)

AA2 <- mean(phi100_t1000[s1000,2] == 2)
Aa2 <- mean(phi100_t1000[s1000,2] == 1)
aa2 <- mean(phi100_t1000[s1000,2] == 0)

sAa = 1 - (Aa2/Aa1)
saa = 1 - (aa2/aa1)

ggplot2(phi100_t1000df, aes(x=x,y=y,col=z1))+geom_point()+scale_color_viridis()


## COMPARSION OF mismatch
phi10 <- get_gsd("p4_analysis/example_data/example_combos/mod-K1_phi10_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv")
phi50 <- get_gsd("p4_analysis/example_data/example_combos/mod-K1_phi50_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv")
phi100 <- get_gsd("p4_analysis/example_data/example_combos/mod-K1_phi100_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv")

phi50_t0 <- get_gsd("~/Anusha/GitHub/LandGenSamp/p1_gnxsims/parallel/LGS_data/t0_data/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.csv")


test <- bind_rows(data.frame(phi50, phi = 0.50), data.frame(phi100, phi = 1), data.frame(phi10, phi = 0.1))
test$mismatch <- rowMeans(data.frame(abs(test$z1 - test$env1), abs(test$z2 - test$env2)))

test3 <- 
  test %>% 
  pivot_longer(c(z1, z2), names_to = "z_name", values_to = "z_value") %>% 
  pivot_longer(c(env1, env2), names_to = "env_name", values_to = "env_value")


ggplot(data = test[test$phi != 0.1,]) +
  #geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", aes(x = env2, y = z2, col = factor(phi), linetype = "z2")) + 
  geom_smooth(method = "lm", aes(x = env1, y = z1, col = factor(phi), linetype = "z1")) + 
  coord_equal() +
  theme_classic()

ggplot(data = test, aes(x = mismatch, fill = factor(phi), col = factor(phi))) +
  geom_density(alpha = 0.5)

ggplot(data = test, aes(x = x, y = y, col = mismatch)) +
  geom_point(alpha = 0.9) +
  facet_grid(~phi) +
  scale_color_viridis(option = "turbo") +
  coord_equal()

test2 <- test %>% pivot_longer(c(z1, z2))
ggplot(data = test2, aes(x = x, y = y, col = value)) +
  geom_point(alpha = 0.9)+
  facet_grid(name~phi) +
  scale_color_viridis(option = "viridis") +
  coord_equal()


phi50_t0 <- get_gsd("~/Anusha/GitHub/LandGenSamp/p1_gnxsims/parallel/LGS_data/t0_data/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.csv")
phi50_t1000 <- get_gsd("~/Anusha/GitHub/LandGenSamp/p1_gnxsims/parallel/LGS_data/mod-K2_phi50_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv")
phi100_t0 <- get_gsd("p1_gnxsims/gnx/LGS_data/t0_data/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.csv")
phi100_t1000 <- get_gsd("p1_gnxsims/gnx/LGS_data/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv")

test <- bind_rows(data.frame(phi50_t0, t = 0, phi = 0.5), data.frame(phi50_t1000, t = 1000, phi = 0.5), data.frame(phi100_t0, t = 0, phi = 1), data.frame(phi100_t1000, t = 1000, phi = 1))
test$mismatch <- rowMeans(data.frame(abs(test$z1 - test$env1), abs(test$z2 - test$env2)))
ggplot(data = test, aes(x = x, y = y, col = mismatch)) +
  geom_point(alpha = 0.9)+
  facet_grid(phi~t) +
  scale_color_viridis(option = "turbo") +
  coord_equal()

test2 <- test %>% pivot_longer(c(z1, z2))
ggplot(data = test2, aes(x = x, y = y, col = value)) +
  geom_point(alpha = 0.9)+
  facet_grid(name~t) +
  scale_color_viridis(option = "viridis") +
  coord_equal()


phi100df <- get_gsd("p1_gnxsims/gnx/LGS_data/t0_data/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-0_spp-spp_0.csv")

phi100 <- read.csv("p1_gnxsims/gnx/LGS_data/stats/mod-K2_phi100_m100_seed1_H50_r60_it-0_spp-spp_0_OTHER_STATS.csv")

phi50 <- read.csv("~/Anusha/GitHub/LandGenSamp/p1_gnxsims/parallel/LGS_data/stats/mod-K2_phi50_m100_seed1_H50_r60_it-0_spp-spp_0_OTHER_STATS.csv")

df <- bind_rows(data.frame(phi100, phi = 100, mismatch = phi100$mean_fit), data.frame(phi50, phi = 50,  mismatch = phi50$mean_fit)) %>% drop_na()

ggplot(df, aes(x = t, y = mismatch,  col = factor(phi))) +
  geom_line() +
  geom_point()

phi100df <- get_gsd("p1_gnxsims/gnx/LGS_data/mod-K2_phi100_m100_seed1_H50_r60_it-0_t-1000_spp-spp_0.csv")




####

files <- list.files("GNX_mod-selection_K2_phi100_m100_seed1_H50_r60/it-0/spp-spp_0/", pattern = ".vcf", full.names = TRUE)

get_genotypes <- function(path){
  
  time <- gsub(".+_t-(\\d+)_.*", "\\1", path)
  
  gsd <- get_gsd(gsub(".vcf", ".csv", path))
  gen <- get_gen(path)
  df <- data.frame(gsd, gen[,1:8], t = time)
  
  set.seed(55)
  s <- sample(nrow(df), 1000)
  df <- df[s, ]
  
  slattice_df <- data.frame(N = nrow(df)*2, matrix(colSums(gen[s,1:8]), nrow = 1), t = time)
  slattice_full <- data.frame(N = nrow(gen)*2, matrix(colSums(gen[,1:8]), nrow = 1), t = time)
  
  return(list(df = df, slattice_df = slattice_df, slattice_full = slattice_full))
}


result <- purrr::map(files, get_genotypes) 
result_dfs <- result %>% list_transpose() %>% map(bind_rows)

  
#The data is a data frame with one row per generation, and two columns named `N` for the total number of chromosomes observed and `N.A` for the total number that carry the allele that is under selection. 
