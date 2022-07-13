source("general_functions.R")

library("here")
library("vcfR")
library("lfmm")
library("raster")
library("automap")

gsd_df <- get_gsd("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/mod-K2_phi50_m25_seed1_H50_r60_it--1_t-1000_spp-spp_0.csv")
gen <- get_gen("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/mod-K2_phi50_m25_seed1_H50_r60_it--1_t-1000_spp-spp_0.vcf")
loci_df <- read.csv("/Users/Anusha/Documents/GitHub/LandGenSamp/p3_methods/data/nnloci_K2_phi50_m25_seed1_H50_r60.csv")


subIDs <- as.character(samples)
gen <- gen[subIDs,]
gsd_df <- gsd_df[subIDs,]


#make grid for kriging
x.range <- c(0,100)
y.range <- c(0,100)
krig_df.grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], len=40),  
                           y=seq(from=y.range[1], to=y.range[2], len=40))
coordinates(krig_df.grd) <- ~x+y
gridded(krig_df.grd) <- TRUE


krig_df = data.frame(x = gsd_df$x,
                     y = gsd_df$y, 
                     prop = gsd_df$env1)

coordinates(krig_df)=~x+y

krig_res <- autoKrige(prop ~ 1, krig_df, krig_df.grd)
krig_raster1 <- raster(krig_res$krige_output)
print(extent(krig_raster1))
extent(krig_raster1) <- c(0,100,0,100)


krig_df = data.frame(x = gsd_df$x,
                     y = gsd_df$y, 
                     prop = gsd_df$env2)

coordinates(krig_df)=~x+y

krig_res <- autoKrige(prop ~ 1, krig_df, krig_df.grd)
krig_raster2 <- raster(krig_res$krige_output)
print(extent(krig_raster2))
extent(krig_raster2) <- c(0,100,0,100)

krig_df = data.frame(x = gsd_df$x,
                     y = gsd_df$y, 
                     prop = gsd_df$env2)

coordinates(krig_df)=~x+y

par(mfrow=c(1,2))
plot(krig_raster1, col=viridis(100))
plot(krig_raster2, col=viridis(100))

tdf <- data.frame(p = gen[1,], gsd_df)
ggplot(tdf)+
  geom_point(aes(x = x, y = y, col = p)) +
  scale_color_viridis(option="viridis")

