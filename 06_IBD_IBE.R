##################################################################################
############      ISOLATION BY DISTANCE AND ENVIRONMENT               ############




rm(list=ls())

setwd("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/IBE")

load("IBE_wc.RData")



##################################################################################
############                Installing Packages                       ############


install.packages("adegenet")
install.packages("hierfstat")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::inst1all("SNPRelate")

install.packages("dartR")
install.packages("rgdal")
install.packages("raster")
install.packages("sp")
install.packages("Hmisc")
install.packages("PerformanceAnalytics")
install.packages("StAMPP")
install.packages("vegan")
install.packages("devtools")
devtools::install_github("mirzacengic/climatedata")
install.packages(c("gdalUtils", "httr", "ncdf4", "qpdf", "raster", "RCurl", "RefManageR", "rgdal", "stringr", "sf", "sp", "svMisc", "utils"), dependencies = TRUE)
install.packages("Rcpp")
install.packages("https://gitlab.rrz.uni-hamburg.de/helgejentsch/climdatdownloadr/-/archive/master/climdatdownloadr-master.tar.gz", repos = NULL, type = "source")
install.packages("tiff")


library(adegenet)
library(hierfstat)
library(dartR)
library(SNPRelate)
library(rgdal)
library(raster)
library(sp)
library(Hmisc)
library(PerformanceAnalytics)
library(StAMPP)
library(vegan)
library(devtools)
library(ClimDatDownloadR)
library(tiff)


##################################################################################
############                  Data Formatting                         ############

#### I will convert a VCF file to Genepop using PGDspider. For that I get the pop
#### file for the conversion

# Reading and converting STRUCTURE file into genind
wren.n5nb.genind <- import2genind("populations2.str", onerowperind=FALSE, n.ind=112, n.loc=4409, col.lab=1, col.pop=2, ask=FALSE)

# Reading the names of Populations
poplab <- read.table("C:/Users/USER/Dropbox/Thesis/Molecular_Wrens/Wren_1/R/poplab.txt", header=FALSE, sep="\t")
colnames(poplab) <- c("sample", "pop")

# I wrote the order of populations by latitude from north to south
poplab$OrdLat <- 1:nrow(poplab)

# Read the pop file used in the denovo assembly in Stacks
# This is the same order for samples in the genind object 
# I need this to get the order of samples
stack.ord <- read.table("popmap_wren_n5nb_sh11.txt", col.names = c("Ind", "OrdStack"))

# Adding a column for the order of samples
stack.ord$OrdStack <- 1:nrow(stack.ord)

# Merging the pop and the order files
pop.lab <- merge(poplab, stack.ord, by.x="sample", by.y="Ind")

# Ordering the labels of Populations as in the genind file
pop.lab <- pop.lab[order(pop.lab$OrdStack, decreasing = FALSE),]


# Setting up the populations in the genind data
strata(wren.n5nb.genind) <- pop.lab[,-3]
setPop(wren.n5nb.genind) <- ~pop


##################################################################################
############  Filtering loci that are not in HWE using dartR package  ############

# Converting to genlight object
gl.wren <- gi2gl(wren.n5nb.genind)

# Creates a loc.metrics object to comply to a dartR object
df.loc <- data.frame(RepAvg=runif(nLoc(gl.wren)), CallRate=1)
gl.wren@other$loc.metrics <- df.loc
gl.hwe.fix <- gl.compliance.check(gl.wren)

# Check for HWE
gl.hwe.rep <- gl.report.hwe(gl.wren, plot=TRUE)

# Filtering out loci not at HWE
gl.hwe <- gl.filter.hwe(gl.hwe.fix)

# Check for Linkage desiquilibrium
gl.hwe.ld <- gl.report.ld(gl.hwe)

# Getting pairwise loci with r2>0.1 with strong linked loci
gl.LD <- subset(gl.hwe.ld,R2>0.1)
View(gl.LD)
# Percentage of pairwise loci linked
(length(gl.LD$loc1)/length(gl.hwe.ld$loc1))*100


### Brief explanation of LD outcome

# Values of D close to zero (D=0) indicates LE, random combination of alleles
# Values of D far from zero (D???0) indicates LD, alleles are correlated
# Values of r2 close to zero (r2=0) indicates LE, random combination of alleles
# Values of r2 far from zero (r2???0) indicates LD, alleles are correlated
# To give you a quick and dirty answer, a D� of 0.8 is high disequilibrium. 
# Basically the two SNPs are coinherited roughly 80% of the time. 
# The reason your r2 is low is that this takes account of allele frequency.
# D� and r2 values are widely used but poorly understood. 
# The current "trend" seems to be to take more notice of the r2 value whereas 
# I feel that the D� is more meaningful and easier to understand. 
# The idea of disequilibrium values is that they are a measure of the non-random 
# association of alleles at two or more loci, i.e how often alleles are coinherited. 
# If two loci are not coinherited at all (they are independent) then both the D� 
# and r2 values will be 0.0 irrespective of either allele frequency. 
# As another example if you had two polymorphisms both with a 50% allele frequency and in 
# total disequilibrium then both the D� and r2 values would be 1.0. 
# However the story changes when the allele frequencies are not the same. 
# For example if you had two polymorphisms, one with a 50% allele frequency and the other with 
# a 1% allele frequency that were still in total disequilibrium then the D� 
# value would be 1.0 but the r2 value would only be 0.01. Basically the D� 
# is saying when the rare allele is present it is always inherited with 
# one particular allele of the 50% polymorphism whereas the r2 is saying it is a rare allele 
# so the vast majority of the time the common allele is not found with it 
# (but only because it is rare, not because it is not in disequilibrium). 
# However, even with such a low r2 this SNP adds nothing to an association study 
# because they are in complete disequilibrium.
# a better understanding with another example where the two SNPs are coinherited 
# about half the time. With both SNPs having a 50% allele frequency the D� value would be 
# about 0.5 and the r2 value would be about 0.25 but if one SNP had a 1% allele frequency 
# then the D� value would still be about 0.5 but the r2 value would only be about 0.005





##################################################################################
############               Estimating pairwise Neifst                 ############


# Converting to hierfstat object
genid.hwe.all <- gl2gi(gl.wren)
hier.hwe.all <- genind2hierfstat(genid.hwe.all)

# Pairwise NeiFst
pairfst.all <- pairwise.neifst(hier.hwe.all)

View(pairfst.all)

# Data set without population in Marañon Valley, East Andes
pop(gl.wren)
gl.west <- gl.drop.pop(gl.wren, pop.list=c("Celedin", "Chachapoyas", "Jaen"))
pop(gl.west)
genid.west <- gl2gi(gl.west)
hier.west.all <- genind2hierfstat(genid.west)

# Pairwise NeiFst
pairfst.west <- pairwise.neifst(hier.west.all)

##################################################################################
############                Isolation by Environment                  ############


## GETTING COORDINATES

# Data with coordinates
labels.all <- read.table("Labels.field.csv", header=TRUE, sep=",", na.strings=TRUE, fileEncoding = "UTF-8-BOM")

# Merging with the data with the populations labels
coor.merged <- merge(poplab, labels.all, by.x="sample", by.y="Ind")

# Merge the data set with coordinates and samples order and filtering C. brunneicapillus
metadata.coor <- merge(stack.ord, coor.merged, by.x="Ind", by.y="sample")[,c(1,2,3,11,12)]
metadata.coor <- metadata.coor[order(metadata.coor$OrdStack,  decreasing = FALSE),]

# Making the coordinates numeric
coor <- cbind(as.numeric(as.character(metadata.coor$Long)), as.numeric(as.character(metadata.coor$Lat)))

# Setting the coordinates system as Long/Lat
cord.dec = SpatialPoints(cbind(coor[,1], coor[,2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))

# Converting to UTM, I'm using UTM zone 17S (epsg=32717)
# Most of the territory of Ecuador is in this zone (17S)
#coor.utm <- spTransform(cord.dec, CRS("+init=epsg:32717"))
#colnames(coor.utm@coords) <- c("X", "Y")


## GETTING WORLDCLIM VARIABLES

# Getting the worldclim variables
r <- getData("worldclim",var="bio",res=10)

# Getting the bioclim variables I need: Annual mean prep and Prep seasonality
r <- r[[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)]]
names(r) <- c("AMT","MDR","ISO","TS","MTWM","MTCM","TAR","MTWetQ","MTDQ","MTWarQ","MTCQ","AMP","PWetM","PDM","PS","PWetQ","PDQ","PWarQ","PCQ")

# Extracting values of climate
values <- extract(r, cord.dec)

# Merging the coordinates, Worldclim values, population labels
coor.clim.raw <- cbind.data.frame(coordinates(coor), values, pop.lab[,c(1:2)])

# Estimating the mean per sampling location
clim.value <- aggregate(. ~ pop, coor.clim.raw[,c(-22)], mean)

# Getting Independent Variables for multiple correlations
IV.all <- clim.value[,c(4:22)]
rcorr(as.matrix(IV.all), type = c("pearson"))

# Climate Variables selected: 
# Annual Mean Temperature (AMT)
# Annual Mean Precipitation (AMP)
# Precipitation Seasonality (PS)

# Getting reduced columns for multiple correlations
IV <- clim.value[,c(4,15,18)]
rcorr(as.matrix(IV), type = c("pearson"))

# Plotting multiple correlations
chart.Correlation(IV, method = "pearson", histogram = TRUE, pch = 16)

# Getting matrix of climate variables
Env <- IV

# Getting distances for environmental variables
Env.dist <- vegdist(scale(Env), "euclid")
Env.dist <- log(Env.dist)
Env.dist <- dist(Env.dist)


# Getting the distances for AMT only
Tem.dist <- vegdist(scale(IV[,1]), "euclid")
Tem.dist <- dist(Tem.dist)


# Getting the distances for AMP only
Prec.dist <- vegdist(scale(IV[,2]), "euclid")
Prec.dist <- dist(Prec.dist)


# Getting the distances for PS only
PS.dist <- vegdist(scale(IV[,3]), "euclid")
PS.dist <- dist(PS.dist)



## GETTING CHELSA VARIABLES

Chelsa.Clim.download(parameter = "bio",
                     bio.var =  c(1:19), 
                     version.var = "1.2", 
                     clipping = TRUE, 
                     clip.extent = c(-81.5, -77.6, -7, 1.5), 
                     buffer = 0, 
                     convert.files.to.asc = FALSE, 
                     stacking.data = TRUE, 
                     combine.raw.zip = FALSE,
                     delete.raw.data = FALSE,
                     save.bib.file = TRUE)


# Setting the directory where the clipped raster was saved
wd <- ("C:\\Users\\USER\\Dropbox\\Thesis\\Molecular_Wrens\\Wren_1\\IBE\\bio\\bio_V1.2\\clipped\\")

# Creating a list with the names of all raster files 
list.raster <- list.files(wd, full.names = TRUE)

# Reading and stacking the rasters
stack.ch <- stack(list.raster)

# Extracting the values for the sampling points I have
values.ch <- extract(stack.ch, cord.dec)
values.ch <- values.ch[,20:38]

# Naming the columns
colnames(values.ch) <- c("AMT","MDR","ISO","TS","MTWM","MTCM","TAR","MTWetQ","MTDQ","MTWarQ","MTCQ","AMP","PWetM","PDM","PS","PWetQ","PDQ","PWarQ","PCQ")


# Merging the coordinates, Chelsa values, population labels
coor.clim.raw.ch <- cbind.data.frame(coordinates(coor), values.ch, pop.lab[,c(1:2)])

# Estimating the mean per sampling location
clim.value.ch <- aggregate(. ~ pop, coor.clim.raw.ch[,c(-22)], mean)

# Getting Independent Variables for multiple correlations
IV.ch.all <- clim.value.ch[,c(4:22)]
rcorr(as.matrix(IV.ch.all), type = c("pearson"))

# Plotting multiple correlations
chart.Correlation(IV.ch.all, method = "pearson", histogram = TRUE, pch = 16)

# Getting reduced columns for multiple correlations
IV.ch <- clim.value.ch[,c(4,15,18)]
rcorr(as.matrix(IV.ch), type = c("pearson"))

# Getting reduced columns for multiple correlations
clim.west <- clim.value.ch[-c(4,5,7),c(1:4,15,18)]
View(clim.west)

# Plotting multiple correlations
chart.Correlation(IV.ch, method = "pearson", histogram = TRUE, pch = 16)

IV.west <- clim.west[,4:6]

# Getting matrix of climate variables
Env.ch <- IV.ch
Env.west <- IV.west

# Getting distances for environmental variables
Env.dist.ch <- vegdist(scale(Env.ch), "euclid")
Env.dist.ch <- log(Env.dist.ch)
Env.dist.ch <- dist(Env.dist.ch)


# Getting distances for environmental variables
Env.dist.west <- vegdist(scale(Env.west), "euclid")
Env.dist.west <- log(Env.dist.west)
Env.dist.west <- dist(Env.dist.west)


# Getting the distances for AMT only
Tem.dist.west <- vegdist(scale(IV.west[,1]), "euclid")
Tem.dist.west <- dist(log(Tem.dist.west))


# Getting the distances for AMP only
Prec.dist.west <- vegdist(scale(IV.west[,2]), "euclid")
Prec.dist.west <- dist(log(Prec.dist.west))


# Getting the distances for PS only
PS.dist.west <- vegdist(scale(IV.west[,3]), "euclid")
PS.dist.west <- dist(log(PS.dist.west))


# Getting matrix of geographical variables
Geo <- clim.west[,c(2:3)]

# Getting geographical distances
Geo.dist <- vegdist(scale(Geo), "euclid")
Geo.dist <- log(Geo.dist)
Geo.dist <- dist(Geo.dist)

# Getting genetic distances
Gen.dist <- pairfst.west
Gen.dist <- Gen.dist[order(rownames(Gen.dist)), ]
Gen.dist <- Gen.dist[, order(colnames(Gen.dist))]
Gen.dist <- Gen.dist/(1-Gen.dist)  
Gen.dist <- dist(Gen.dist)


# Checking the dimensions
dim(Env.dist.west)
dim(Geo.dist)
dim(Gen.dist)

View(Gen.dist)

levels <- c("CFPS", "CFPN", 
"CFPS", "CFF", "CFF", "CZB", "CFF", "CZB", 
"CFPN", "CFPN", "CFPN", "CZB", 
"CZB", "CZB", "CFPS", "CFPS")

## RUNNING PARTIAL MANTEL TESTS WITH WORLDCLIM CLIMATE DATA

# Running the Partial Mantel Test between Genetic and Environment controlled by Geography
mantel.partial(Gen.dist, Env.dist.west, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))

# Running the Partial Mantel Test between Genetic and Temperature controlled by Geography
mantel.partial(Gen.dist, Tem.dist, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Partial Mantel Test between Genetic and Precipitation controlled by Geography
mantel.partial(Gen.dist, Prec.dist, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Partial Mantel Test between Genetic and Precipitation Seasonality controlled by Geography
mantel.partial(Gen.dist, PS.dist, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Partial Mantel Test between Genetic and Geography controlled by Environment
mantel.partial(Gen.dist, Geo.dist, Env.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Mantel Test between Genetics and Environment
mantel(Gen.dist, Env.dist, method = "pearson", permutations = 10000, 
       strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Mantel Test between Genetics and Latitude
mantel(Gen.dist, Geo.dist, method = "pearson", permutations = 10000, 
       strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))



## RUNNING PARTIAL MANTEL TESTS WITH CHELSA CLIMATE DATA FOR WEST

# Running the Partial Mantel Test between Genetic and Environment controlled by Geography
mantel.partial(Gen.dist, Env.dist.west, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))

               
# Running the Partial Mantel Test between Genetic and Temperature controlled by Geography
mantel.partial(Gen.dist, Tem.dist.west, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Partial Mantel Test between Genetic and Precipitation controlled by Geography
mantel.partial(Gen.dist, Prec.dist.west, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Partial Mantel Test between Genetic and Precipitation Seasonality controlled by Geography
mantel.partial(Gen.dist, PS.dist.west, Geo.dist, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Partial Mantel Test between Genetic and Geography controlled by Environment
mantel.partial(Gen.dist, Geo.dist, Env.dist.west, method = "pearson", permutations = 10000, 
               strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


# Running the Mantel Test between Genetics and Environment
mantel(Gen.dist, Env.dist.west, method = "pearson", permutations = 10000, 
       strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))


##################################################################################
############                         GDM WEST                         ############

# Changing some columns names
colnames(clim.west)[1:3] <- c("site", "x", "y")

# Selecitng the columns to be used in the GDM
pred.dat <- clim.west[,c("site", "x", "y", "AMT", "AMP", "PS")] 
colnames(pred.dat)
head(pred.dat)

pred.dat.sc <- pred.dat
pred.dat.sc[,4:6] <- scale(pred.dat[,4:6])
head(pred.dat.sc)

# need to add a site ID column
Gen.dist.m <- pairfst.west

# Order matrix Gen.dist.m by row names and column names
Gen.dist.m <- Gen.dist.m[, order(colnames(Gen.dist.m))]
Gen.dist.m <- Gen.dist.m[order(rownames(Gen.dist.m)), ]


site <- rownames(Gen.dist.m)
Gen.dist.m <- cbind(site, Gen.dist.m)
View(Gen.dist.m)



install.packages("gdm")
library(gdm)

# Getting the data in format for GDM using the fuction formatsitepair
gdm.wrens <- formatsitepair(bioData = Gen.dist.m, bioFormat=3, siteColumn="site", 
                            XColumn="x", YColumn="y", dist="euclidian", predData=pred.dat.sc)
head(gdm.wrens)

# Running the GDM
gdm.m1 <- gdm(data=gdm.wrens, geo=TRUE)

# Checking the model
summary(gdm.m1)

# Checking the model
length(gdm.m1$predictors)

# Checking the model
plot(gdm.m1, plot.layout=c(2,3))

# Gettign variables importance
modTest <- gdm.varImp(gdm.wrens, geo=TRUE, predSelect=TRUE, nPerm=10000, cores=10)

# Get percentages of variance explained by each predictor
# How to get the percentage of a percentage

(modTest[[2]]*100)/sum(modTest[[2]])

# Partition variance by component environment and geography
# Make list of variable sets for partitioning
varSet <- vector("list",3)
names(varSet) <- c("temp", "amp", "ps")
varSet$temp <- c("AMT")
varSet$amp <- c("AMP")
varSet$ps <- c("PS")

varSet

# run the function to partition temperature, precipitation, and space (partSpace=TRUE)
scgPart <- gdm.partition.deviance(sitePairTable=gdm.wrens, varSets=varSet, partSpace=FALSE)

?gdm.partition.deviance



##################################################################################
############                         GDM ALL                          ############

# Changing some columns names
colnames(clim.west)[1:3] <- c("site", "x", "y")

# Selecitng the columns to be used in the GDM
pred.dat <- clim.west[,c("site", "x", "y", "AMT", "AMP", "PS")] 
colnames(pred.dat)
head(pred.dat)

pred.dat.sc <- pred.dat
pred.dat.sc[,4:6] <- scale(pred.dat[,4:6])
head(pred.dat.sc)

# need to add a site ID column
Gen.dist.m <- pairfst.west

# Order matrix Gen.dist.m by row names and column names
Gen.dist.m <- Gen.dist.m[, order(colnames(Gen.dist.m))]
Gen.dist.m <- Gen.dist.m[order(rownames(Gen.dist.m)), ]


site <- rownames(Gen.dist.m)
Gen.dist.m <- cbind(site, Gen.dist.m)
View(Gen.dist.m)



install.packages("gdm")
library(gdm)

# Getting the data in format for GDM using the fuction formatsitepair
gdm.wrens <- formatsitepair(bioData = Gen.dist.m, bioFormat=3, siteColumn="site", 
                            XColumn="x", YColumn="y", dist="euclidian", predData=pred.dat.sc)
head(gdm.wrens)

# Running the GDM
gdm.m1 <- gdm(data=gdm.wrens, geo=TRUE)

# Checking the model
summary(gdm.m1)

# Checking the model
length(gdm.m1$predictors)

# Checking the model
plot(gdm.m1, plot.layout=c(2,3))

# Gettign variables importance
modTest <- gdm.varImp(gdm.wrens, geo=TRUE, predSelect=TRUE, nPerm=10000, cores=10)

# Get percentages of variance explained by each predictor
# How to get the percentage of a percentage

(modTest[[2]]*100)/sum(modTest[[2]])

# Partition variance by component environment and geography
# Make list of variable sets for partitioning
varSet <- vector("list",3)
names(varSet) <- c("temp", "amp", "ps")
varSet$temp <- c("AMT")
varSet$amp <- c("AMP")
varSet$ps <- c("PS")

varSet

# run the function to partition temperature, precipitation, and space (partSpace=TRUE)
scgPart <- gdm.partition.deviance(sitePairTable=gdm.wrens, varSets=varSet, partSpace=FALSE)

?gdm.partition.deviance

##################################################################################
#### SAVING ####
save.image("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/IBE/IBE.RData")


##################################################################################