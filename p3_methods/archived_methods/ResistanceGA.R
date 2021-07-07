library("here") #paths
library("vcfR") #read VCF files

#note: can't install ResistanceGA using install.packages()
#Instead uncomment these lines to install:
#install.packages("remotes")
#remotes::install_github("wpeterman/ResistanceGA")
library("ResistanceGA") #ResistanceGA

############
#   Data   #
############
#DEFINE NLOCI
nloci <- 1000

#NEEDS TO BE MODIFIED FOR FUTURE REAL DATA
gea_df <- read.csv(here("data", "mod-yosemite_demo_it--1_t-500_spp-Sceloporus graciosus.csv"))
gea_df$tmp <- as.numeric(stringr::str_extract(gea_df$e, '(?<=\\[)[^,]+(?=,)'))
head(gea_df)

vcf <- read.vcfR(here("data","mod-yosemite_demo_it--1_t-500_spp-Sceloporus graciosus.vcf"))
x <- vcfR2genlight(vcf)
gen <- as.matrix(x)

#FOR TESTING USE SAMPLE OF 100 (!!!COMMENT OUT!!!)
set.seed(42)
s <- sample(1:nrow(gea_df),100)
gea_df <- gea_df[s,]
gen <- gen[s,]

#Environmental data
envRast <- raster(here("data","tmp_1980-2010_90x90.tif"))
extent(envRast) <- c(0, 90, 0, 90)
#(?MODIFY?)
#Currently have to flip coordinates since in gnx 
#the coordinates start at 0,0 in the upper left corner 
#whereas in rasters the coordinates start at 0,0 in the lower left corner.
#Alternatively can  set extent to c(0,90,-90,0) and make y coordinates negative
envRast <- flip(envRast, direction='y') 
plot(envRast)

#Environmental data
sdmRast <- raster(here("data","sdm_1980-2010_90x90.tif"))
extent(sdmRast) <- c(0, 90, 0, 90)
#(?MODIFY?)
#Currently have to flip coordinates since in gnx 
#the coordinates start at 0,0 in the upper left corner 
#whereas in rasters the coordinates start at 0,0 in the lower left corner.
#Alternatively can  set extent to c(0,90,-90,0) and make y coordinates negative
sdmRast <- flip(sdmRast, direction='y') 

#make coords object
coords <- cbind(gea_df$x, gea_df$y)
#Check to make sure raster orientation and coordinates from dataframe align
if(cor(extract(envRast,coords), gea_df$tmp) =! 1){"Raster orientation and df coordinates do not align!"}

#CALCULATING PC BASED GENETIC DISTANCE
#perform PCA
pc <- prcomp(gen)
#Calculate PC distance based on first three PCs (?MODIFY?)
pc_dist <- as.matrix(dist(pc$x[,1:64], diag = TRUE, upper = TRUE)) #CHANGE NUMBER OF PCS? (see Shirk et al. 2016:  10.1111/1755-0998.12684)
gdist.response <- lower(pc_dist) #Converts dist matrix into vector of lower half of pairwise matrix for gendist

geo_dist <- as.matrix(dist(gea_df[,c("x","y")], diag = TRUE, upper = TRUE))
geo_dist <- lower(geo_dist)
plot(geo_dist,gdist.response)
###################
#  Resistance GA  #
###################

#ResistanceGA code (see p10 for more details: http://petermanresearch.weebly.com/uploads/2/5/9/2/25926970/resistancega.pdf)
#Setup inputs
#By default, only monomolecular transformations will be assessed for continuous surfaces unless otherwise specified. 
GA.inputs <- GA.prep(ASCII.dir = envRast,
                     Results.dir = here("outputs/"),
                     select.trans="A",
                     maxiter = 10) #!!!FOR TESTING SET MAXITER TO LOW NUMBER CHANGE IN FINAL SCRIPT!!!
gdist.inputs <- gdist.prep(n.Pops = nrow(coords),
                     response = gdist.response,
                     samples = coords)

#run opimization using GA
SS_RESULTS <- SS_optim(gdist.inputs = gdist.inputs,
                       GA.inputs = GA.inputs)

#save solution to PARM object
PARM <- SS_RESULTS$ga[[1]]@solution
results <- SS_RESULTS$ContinuousResults[c(9:11)]

#plot trans
Plot.trans(PARM, 
           envRast, 
           transformation=results$Equation)

#get resistance surface
resRast <- Resistance.tran(transformation = results$Equation, r = envRast, shape = results$shape, max = results$max)
#plot resistance surface
library("wesanderson")
pal <- wes_palette("Zissou1", 100, type = c("continuous"))
palr <- wes_palette("Zissou1", 100, type = c("continuous"))

par(mfrow=c(1,3))
plot(resRast, col=pal, main="Result")
plot(envRast, col=rev(pal), main="Env")
plot(sdmRast, col=pal, main="SDM")

