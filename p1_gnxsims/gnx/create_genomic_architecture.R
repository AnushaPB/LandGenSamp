#code to create a genomic architecture file
library(here)
#from gnx docs:
#CSV file with loci as rows and ‘locus’, ‘p’, ‘dom’, ‘r’, ‘trait’, and ‘alpha’ as columns 
#stipulating the locus numbers, starting allele frequencies, dominance values, inter-locus recombination rates, trait names, 
#and effect sizes of all loci; values can be left blank if not applicable).

#total number of loci
nloci <- 10000

#create data frame
df <- data.frame(locus = 0:(nloci-1),
           p = 0.5,
           #codominance = 0 (default)
           dom = 0,
           r = 0.5,
           trait = NA,
           alpha = 0)

#change first r value to 0 (required by gnx - see docs)
df$r[1] <- 0

#specify which loci belong to which traits
#in this case there are two traits, each with 4 loci
#their location does not matter because there is no linkage
df$trait[1:4] <- "trait_1"
df$trait[5:8] <- "trait_2"
#set alpha for trait related loci to cover phenotype spectrum
df$alpha[1:8] <- rep(c(0.25, -0.25), 4)

#IMPORTANT NOTE: because of the python indexing, the loci will be 0:3 for trait 1 and 4:7 for trait 2
#I have accounted for this in the rest of the R code

#write out file
path <- here("p1_gnxsims", "gnx", "genomic_architecture.csv")
write.csv(df, path, row.names = FALSE)

# genomic architecture 2 (4 loci)
df$trait <- NA
df$alpha <- NA
df$trait[1:2] <- "trait_1"
df$trait[3:4] <- "trait_2"
#set alpha for trait related loci to cover phenotype spectrum
df$alpha[1:4] <- rep(c(0.5, -0.5), 2)
path <- here("p1_gnxsims", "gnx", "genomic_architecture2.csv")
write.csv(df, path, row.names = FALSE)

# genomic architecture 3 (one trait)
df$trait <- NA
df$alpha <- NA
df$trait[1:4] <- "trait_1"
#set alpha for trait related loci to cover phenotype spectrum
df$alpha[1:4] <- rep(c(0.25, -0.25), 2)
path <- here("p1_gnxsims", "gnx", "genomic_architecture3.csv")
write.csv(df, path, row.names = FALSE)