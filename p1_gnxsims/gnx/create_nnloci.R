library(here)
nnloci <- data.frame(trait1 = 0:3, trait2 = 4:7)
write.csv(nnloci, here("p1_gnxsims", "gnx", "LGS_data", "nnloci.csv"), row.names = FALSE)
