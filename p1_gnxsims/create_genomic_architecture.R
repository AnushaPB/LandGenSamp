#code to create a genomic architecture file
library(here)
#from gnx docs:
#CSV file with loci as rows and ‘locus’, ‘p’, ‘dom’, ‘r’, ‘trait’, and ‘alpha’ as columns 
#stipulating the locus numbers, starting allele frequencies, dominance values, inter-locus recombination rates, trait names, 
#and effect sizes of all loci; values can be left blank if not applicable).

#total number of loci
nloci <- 100000

#create data frame
df <- data.frame(locus = 1:nloci,
           p = 0.5,
           #codominance = 0 (default)
           dom = 0,
           r = 0,
           trait = NA,
           alpha = NA)

#specify which loci belong to which traits
#in this case there are two traits, each with 4 loci
#their location does not matter because there is no linkage
df$trait[1:4] <- "trait_1"
df$trait[5:8] <- "trait_2"

#write out file
path <- here(dirname(getwd()), "p1_gnxsims", "genomic_architecture.csv")
write.csv(df, path, row.names = FALSE)

