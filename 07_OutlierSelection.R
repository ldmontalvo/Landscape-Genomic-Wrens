################################################################################
################                OUTLIER DETECTION               ################
################################################################################

rm(list=ls())

##### Setting the working directory
setwd("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines")


##### Load the outcomes so far
load("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/Clines.RData")

################################################################################
############                Installing Packages                       ##########


### INSTALL PACKAGES

##### Spatial Analysis Packages
install.packages(c("spdep","raster","rgdal","ClimDatDownloadR", "rgeos", "envirem"), dependencies=TRUE)

##### Statistical Packages
install.packages(c("adespatial","Hmisc","PerformanceAnalytics","coda"), dependencies=TRUE)

##### Genomic and Ecological Analysis Packages
install.packages(c("vcfR"), dependencies=TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LEA")

##### Plotting and data management
install.packages(c("RColorBrewer", "viridis", "reshape2"), dependencies=TRUE)

library(envirem) # Have to be called first to avoid conflict with raster packages
library(spdep)
library(adespatial)
library(raster)
library(viridis)
library(vcfR)
library(LEA)
library(rgdal)
library(Hmisc)
library(PerformanceAnalytics)
library(RColorBrewer)
library(coda)
library(ClimDatDownloadR)
library(rgeos)##### Load packages
library(reshape2)

################################################################################
################                 Data Formatting                ################


##### Calling the vcf and converting to genotype file
wren.vcf.n5nb <- "C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/FastQC_Stacks/denovo_n5nb_sh11/populations.snps.vcf"
vcf2geno(wren.vcf.n5nb, "geno.wren2.txt")

##### Converting Structure to genotype file
struct2geno("populations.structure", ploidy=2, FORMAT=2, extra.row=1, extra.column=2)
## This command give us a file with more than 4409 snps. I'm using here geno.wren2.txt

##### Reading the file produce with vcf2geno
geno.wren <- read.table(text = gsub("", "\t", readLines("geno.wren2.txt")))

##### Transposing the matrix
geno.wren.t <- t(geno.wren)


### GETTING METADATA (coordinates, sampling locations, labels)

##### Reading file with labels
indtable <- read.table("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/R/Genepop/k4popfile.txt", header = TRUE, sep="\t", na.strings=TRUE, fileEncoding = "UTF-8-BOM")

##### Data with coordinates
labels.all <- read.table("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/R/Labels.field.csv", header=TRUE, sep=",", na.strings=TRUE, fileEncoding = "UTF-8-BOM")

##### Read the pop file used in the denovo assembly in Stacks 
##### I need this to get t he order of samples
stack.ord <- read.table("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/R/popmap_wren_n5nb_sh11.txt", col.names = c("Ind", "ord"))

library(dplyr)

##### Getting everything in the same dataframe
meta.data <- merge(stack.ord, labels.all) %>%
  merge(indtable, by.x="Ind", by.y="sample")

##### Selecting the columns we need and ordering
meta.data <- meta.data[,c(1:2,9,10,21)]
meta.data <- meta.data[order(meta.data$ord.x),]


### GETTING COORDINATES IN UTM FORMAT

##### Making the coordinates numeric
coor <- cbind(as.numeric(as.character(meta.data$Long)), as.numeric(as.character(meta.data$Lat)))

##### Setting the coordinates system as Long/Lat
cord.dec = SpatialPoints(cbind(coor[,1], coor[,2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))

##### Converting to UTM, I'm using UTM zone 17s (epsg=32717)
##### Most of the territory is in this zone
coor.utm <- spTransform(cord.dec, CRS("+init=epsg:32717"))
colnames(coor.utm@coords) <- c("X", "Y")

##### Getting the metadata with coordinates in UTM
meta.data <- cbind(meta.data, coor.utm@coords)

##### Getting the final data frame
geno.meta <- cbind(meta.data[,c(1,6:7,5)], geno.wren.t)
View(geno.meta)

### GETTING CHELSA VARIABLES


##### Most of this code came from http://envirem.github.io/ENVIREM_tutorial.html
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

##### Setting the directory where the clipped raster was saved
wd <- ("C:\\Users\\Daniel\\Dropbox\\Thesis\\Molecular_Wrens\\Radseq\\IBE\\bio\\bio_V1.2\\clipped\\")

##### Creating a list with the names of all raster files 
list.raster <- list.files(wd, full.names = TRUE)

##### Reading and stacking the rasters
stack.ch <- stack(list.raster)

##### Extracting the values for the sampling points I have
values.ch <- extract(stack.ch, cord.dec)
values.ch <- values.ch[,20:38]

##### Naming the columns
colnames(values.ch) <- c("AMT","MDR","ISO","TS","MTWM","MTCM","TAR","MTWetQ","MTDQ","MTWarQ","MTCQ","AMP","PWetM","PDM","PS","PWetQ","PDQ","PWarQ","PCQ")

coor.chelsa <- cbind(coor.lag$Localidad, values.ch2)
write.table(coor.chelsa, "coor_chelsa.txt", sep="\t", quote = FALSE)

##### Merging the coordinates, Chelsa values, population labels
coor.clim.raw.ch <- cbind.data.frame(coordinates(coor), values.ch)
colnames(coor.clim.raw.ch)[1:2] <- c("Long", "Lat")

##### Getting Independent Variables for multiple correlations
IV.ch.all <- coor.clim.raw.ch[,c(3:21)]
rcorr(as.matrix(IV.ch.all), type = c("pearson"))

##### Plotting multiple correlations
chart.Correlation(IV.ch.all, method = "pearson", histogram = TRUE, pch = 16)

##### Getting reduced columns for multiple correlations
IV.ch <- coor.clim.raw.ch[,c(3,14,17)]
rcorr(as.matrix(IV.ch), type = c("pearson"))

##### Plotting multiple correlations
chart.Correlation(IV.ch, method = "pearson", histogram = TRUE, pch = 16)

##### Getting matrix of climate variables
AMT <- IV.ch[,1]
AMP <- IV.ch[,2]
PS <- IV.ch[,3]


##### Getting the map of precipitation
prec.raster <- stack.ch[[31]]

##### Making the palette of colors I will use
##### colr <- brewer.pal(11, 'RdYlBu') %>%
##### colorRampPalette()

# Combining meta data and climate variables
climate <- cbind(meta.data, coor.clim.raw.ch)


# Writing the climate date for the samples for other analyses like plumage coloration patterns
write.table(climate, "climate.txt", sep="\t")


### GETTING MONTHLY CHELSA VARIABLES


Chelsa.Clim.download(save.location="../IBE/Chelsa",
                     parameter = c("prec", "temp", "tmax", "tmin"),
                     clipping = TRUE, 
                     clip.extent = c(-81.5, -77.6, -7, 1.5), 
                     buffer = 0, 
                     convert.files.to.asc = FALSE, 
                     stacking.data = TRUE, 
                     combine.raw.zip = FALSE,
                     delete.raw.data = FALSE,
                     save.bib.file = TRUE)


### SOLAR RADIATION RASTERS


##### read in a climatic raster for use as a template
rasterTemplate <- raster(list.raster.all[[1]])

##### calculate monthly solar radiation, defined for the year 2050, output to the current directory
wd.solrad <- ("C:\\Users\\USER\\Dropbox\\Thesis\\Molecular_Wrens\\Wren_1\\IBE\\Chelsa\\solrad\\")

##### Getting the Solar Radiation Raster
ETsolradRasters(rasterTemplate = rasterTemplate, year = 65, outputDir = wd.all, overwrite = TRUE)

##### Verifying is they are properly named
verifyFileStructure(wd.solrad, returnFileNames = FALSE)

##### Setting the directory where the clipped rasters were saved
wd.all <- ("C:\\Users\\USER\\Dropbox\\Thesis\\Molecular_Wrens\\Wren_1\\IBE\\Chelsa\\")

##### Creating a list with the names of all raster files 
list.raster.all <- list.files(wd.all, pattern=c("clipped.tif$|_solrad_"), recursive=TRUE, full.names = TRUE)

##### Reading the first precipitation raster from the list
precRaster <- raster(grep(list.raster.all, pattern =  '_prec_', value=TRUE)[1])

##### Check the file
plot(precRaster)

##### Creating folder for input and outcomes
envirem.in <- "C:/Users/USER/Dropbox/Thesis/Molecular_Wrens/Wren_1/IBE/Chelsa/envirem/input"
dir.create(envirem.in)
envirem.out <- "C:/Users/USER/Dropbox/Thesis/Molecular_Wrens/Wren_1/IBE/Chelsa/envirem/output"
dir.create(envirem.out)

##### Create folder with input
for(i in list.raster.all){file.copy(list.raster.all, envirem.in)}

##### Get all the raster in stack for other analyses
#chelsalist <- list.files(envirem.in, pattern="CHELSA", full.names = TRUE)
#chelsastack <- stack(chelsalist)
#solradfiles <- list.files(envirem.in, pattern = "et_solrad", full.names = TRUE)
#solradstack <- stack(solradfiles)

##### Checking the names
#verifyRasterNames(masterstack=chelsastack, solradstack = solradstack, returnRasters = TRUE)
#verifyFileStructure(envirem.in, returnFileNames = TRUE)
#plot(solradstack[[1]])
#plot(raster(envirem.list[[1]]))
# I found issues with the names of the raster, I had to change the names manually in the folder 

##### Assign names
assignNames(tmin = 'CHELSA_tmin10_##',
            tmax = 'CHELSA_tmax10_##',
            tmean = 'CHELSA_temp10_##',
            precip = 'CHELSA_prec_##',
            solrad = 'et_solrad_##')
#assignNames(reset = TRUE)

##### Check assigned names
envirem::varnames()

##### Generating envirem variables
generateRasters(var="all",
                maindir=envirem.in,
                outputDir=envirem.out,
                nTiles=1)

##### Plotting the results
enviremRasters <- list.files(envirem.out, pattern = '\\.tif$', full.names = TRUE)
enviremRasters <- stack(enviremRasters)

par(mfrow = c(5, 4), mar = c(0.5, 0.5, 2, 0))

for (i in 1:nlayers(enviremRasters)) {
  plot(enviremRasters[[i]], col = inferno(100), box = FALSE, axes = FALSE)
  title(main = names(enviremRasters)[i])
}


##### Extracting the values for the sampling points
values.envirem <- extract(enviremRasters, cord.dec)
cmi.all <- values.envirem[,3]




### GETTING THE DATA IN DIFFERENT DATAFRAME

# Getting coordinate alone
Coord <- data.matrix(geno.meta[,1:2])

# Sampling locations for individuals
Samloc <- as.data.frame(geno.meta[,3])

# Genotypes alone
Loci <- data.matrix(geno.meta[,4:ncol(geno.meta)])


################################################################################
################     Moran spectral outlier detection (MSOD)    ################

#### From: https://popgen.nescent.org/2016-12-13_MEM_outlier.html


##### Neighbor definition
nb <- graph2nb(gabrielneigh(Coord), sym=TRUE)

##### Spatial weights matrix
listW <- nb2listw(nb, style="W")

##### Add longlat=T for lat/long coordinates
disttri <- nbdists(nb, Coord)

##### Use inverse distance weights
fdist <- lapply(disttri, function(x) (x+0.000001)^(-1))

##### Revised spatial weights matrix
listW <- nb2listw(nb, glist=fdist, style="W") 

##### Eigen analysis
tmp <- scores.listw(listW, MEM.autocor = "all") 

##### MEM eigenvectors and eigenvalues
mem <- list(vectors = as.matrix(tmp), values = attr(tmp, "values")) 

##### Rescale eigenvalues to Moran's I
mem$values <- mem$values / abs(sum(mem$values))     

##### Add the individuals and Gabriel graph
# plot(nb, coords=Coord, col=1, pch=16, cex=0.8, add=T)  

##### Plot the selection surface
prec.raster.utm <- projectRaster(prec.raster, crs=crs("+init=epsg:32717"))
par(mar=c(1, 1, 1, 1))
plot(prec.raster.utm, axes=F, legend=F, box=F,                 
     col=c("gray70","gray98")) 

##### Add the individuals and Gabriel graph
plot(nb, coords=Coord, col=1, pch=16, cex=0.8, add=T)  



### CORRELATIONS BETWEEN THE LOCI AND MEM AXES


##### Calculate R.YV, which contains for each locus the vector of its correlations with all MEM axes. 

##### R.YV = Correlation with MEM axes
R.YV <- cor(Loci, mem$vectors, use="pairwise.complete.obs")     

##### S = Average power spectrum
S    <- apply(R.YV^2, 2, mean)                              


### PLOT LOT POWER SPECTRA


##### First locus
barplot((R.YV^2)[1,], ylim=c(0, 0.12))  

##### Second locus
barplot((R.YV^2)[2,], ylim=c(0, 0.12))  

##### Average power spectrum (all 100 loci)
barplot(S, ylim=c(0, 0.12))             

##### Cutoffs (can be modified!)
##### Three cutoff
cutoffs <- abs(qnorm(c(0.05, 0.01, 0.001)/2))  


### CALCULATE Z-SCORE FOR POWER SPECTRA


##### Subtract average power spectrum from each locus.
Dev <- sweep(R.YV^2, 2, S, "/") - 1 

##### Set positive deviations to zero.
Dev[Dev > 0] <- 0        

##### Sum of negative deviations
Dev <- apply(Dev, 1, sum)     

##### Standardize
z <- scale(Dev)                                


### PLOT Z-SCORE FOR POWER SPECTRA

##### Plot the z-scores
plot(z, ylim=c(-7,5))   

##### Add lines for the three cutoffs
for(h in 1:length(cutoffs))                    
{
  lines(c(0,4409), rep(cutoffs[h],2), lty=h)
  lines(c(0,4409), rep(-cutoffs[h],2), lty=h)
}


### MEM CORRELATION FOR ENV AND COORDINATES (as spurious predictors)


##### Just the middle cutoff of 0.05
cutoff.msod <- cutoffs[2]                              

##### Candidate loci at this cutoff
Candidates.msod <- c(1:4409)[abs(z)>cutoff.msod]        

##### Set a cutoff & number of permutations for MSR

# cutoff.msr <- 0.05    # Set a less stringent cutoff

##### Set number of permutations for MSR test (may choose e.g. 499 or 999)
nPerm <- 1000

R.XV.AMT <- cor(scale(AMT), mem$vectors)
R.XV.AMP <- cor(scale(AMP), mem$vectors)
R.XV.PS <- cor(scale(PS), mem$vectors)
R.XV.xcoord <- cor(Coord[,1], mem$vectors)
R.XV.ycoord <- cor(Coord[,2], mem$vectors)


### FUNCTION TO PERFORM MSR TEST


get.pvalue.msr <- function(r.XV=R.XV, r.YV=R.YV, nPerm=1000)
{
  R.XV.rand <- matrix(r.XV, nPerm, ncol(r.XV), byrow=TRUE) 
  R.XV.rand <- R.XV.rand * sample(c(-1,1), length(R.XV.rand), replace=TRUE)
  Cor.obs <- abs(as.vector(r.YV %*% t(r.XV)))
  Cor.rand <- abs(r.YV %*% t(R.XV.rand))
  P.values.MSR <- apply((cbind(Cor.obs,Cor.rand) >= Cor.obs), 1, mean)
  P.values.MSR
}


##### MSR test for candidate outlier loci detected by MSOD

b.AMT <- get.pvalue.msr(r.XV=R.XV.AMT, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.AMP <- get.pvalue.msr(r.XV=R.XV.AMP, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.PS <- get.pvalue.msr(r.XV=R.XV.PS, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.X <- get.pvalue.msr(r.XV=R.XV.xcoord, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)
b.Y <- get.pvalue.msr(r.XV=R.XV.ycoord, r.YV=R.YV[Candidates.msod,], nPerm=nPerm)

cutoff.msr <- 0.01

print(paste("Loci significantly associated with AMT: ", names(b.AMT)[b.AMT < cutoff.msr]))
print(paste("Loci significantly associated with AMP: ", names(b.AMP)[b.AMP < cutoff.msr]))
print(paste("Loci significantly associated with PS: ", names(b.PS)[b.PS < cutoff.msr]))
print(paste("Loci significantly associated with X: ", names(b.X)[b.X < cutoff.msr]))
print(paste("Loci significantly associated with Y: ", names(b.Y)[b.Y < cutoff.msr]))


msod.01 <- c(names(b.AMP)[b.AMP < cutoff.msr],
             names(b.PS)[b.PS < cutoff.msr],
             names(b.Y)[b.Y < cutoff.msr])

sig.snp <- as.data.frame(table(msod.01))
sig.snp[,1]

### Get chromosomes names from vcf
library(vcfR)
vcf <- read.vcfR(wren.vcf.n5nb, verbose = FALSE )

## Getting the names of chromosomes and other data
chrm <- as.data.frame(getFIX(vcf))
head(chrm)
dim(chrm)

## Adding columns of continous number from 1 to 4409
chrm$loci <- 1:4409

## Select row in chrm that correspond to the loci in sig.snp[,1]
chrm.sig <- chrm[chrm$loci %in% sig.snp[,1],]


### PLOT THE SPATIAL DISTRIBUTIONOF ALLELES


allele_colors <- viridis(3)          # Set allele colors

##### Making the palette of colors I will use
prec.col <- brewer.pal(11, 'RdYlBu') %>%
  colorRampPalette()

##### Allele correlated with Precipitation
par(mfrow=c(1, 2))
plot(prec.raster.utm, axes=F, legend=F, box=F, col=prec.col(20)) 
plot(nb, coords=Coord, add=T) 
points(Coord, col = allele_colors[Loci[,123] + 1], pch=20)
title("L123 - True Positive Detection")

##### Allele correlated with Latitude
plot(prec.raster.utm, axes=F, legend=F, box=F, col=prec.col(20))
plot(nb, coords=Coord, add=T)
points(Coord, col = allele_colors[Loci[,946] + 1], pch=20)
title("L946 - True Positive Detection")

##### Color codes
# 2 = yellow
# 1 = green
# 0 = purple



################################################################################
################                 Bayescenv Analysis             ################
 

#### RUN BAYESCENV WITH  AMP ####

#### Preparing Environmental variables

##### Getting the data I need
bayes.amp.all <- cbind(meta.data[,c(1,5)], coor.clim.raw.ch[,c(1,2,14)])
colnames(bayes.amp.all)
View(coor.clim.raw.ch)

##### BAYESCENV get the order of sampling location for order of appearance
ord.pop <- as.data.frame(cbind(as.factor(unique(bayes.amp.all$Pop)), seq(1,16,1)))
colnames(ord.pop) <- c("alphab", "appear")
ord.pop <- ord.pop[order(ord.pop$alphab, decreasing=FALSE),]

##### Getting the average of AMP per sampling locations
bayes.IV <- bayes.amp.all %>%
  group_by(Pop) %>%
  summarise(AMP_mean = mean(AMP, na.rm = TRUE))

##### Getting the order as in BAYESCENV
bayes.IV <- cbind(bayes.IV, ord.pop$appear)
colnames(bayes.IV) <- c('samloc', "AMP_mean", "appear")
bayes.IV <- bayes.IV[order(bayes.IV$appear, decreasing=FALSE),]

##### Scaling the variables to get standarized differences with the mean
bayes.IV$amp_sc <- scale(bayes.IV$AMP_mean)

# This is the same numbers for the second run of structure. 
# Same sampling locations as input of bayescenv

AMPp <- as.vector(bayes.IV$amp_sc)
AMPp <- as.numeric(format(round(AMPp, 2), nsmall=2))
AMPp <- t(scale(AMPp))

##### Write the vector to a .txt file for bayescenv 
write.table(AMPp, "./runs/AMPp.txt", row.names = FALSE, col.names = FALSE)

##### This AMP will be the input for bayescenv-1.1

### CHECKING CONVERGENCE AND Q-VALUES

##### I ran bayescenv in Ubuntu. After that I just upload the outcome files

##### After running bayescenv-1.1 check for convergence

##### Read the .sel file where the values for the iterations are
chain <- read.table("./runs/run.amp.res1.sel", header = TRUE)

##### Adapt thin to its actual value
chain <- mcmc(chain, thin = 10)

##### To test for convergence
heidel.diag(chain)

##### To compute effective sample size
effectiveSize(chain)    

##### To look for auto-correlation
autocorr.diag(chain) 

##### To plot the "trace" and the posterior distribution
plot(chain)

#View(chain)  
#summary(chain)

##### checking output
run.fst <- read.table("./runs/run.amp.res1_fst.txt", header=TRUE)
View(run.fst)

##### Check the lowest prob of getting a right model and g parameter
min(run.fst$PEP_g)


#### RUN BAYESCENV WITH RESIDUALS OF AMT VS LATITUDE FOR INDEPENDENT EFFECTS OF TEMPERATURE  #####

#### Preparing Environmental Variables (Temperature)


##### Run a lineal regression between AMT and Latitude
amt.ind <- lm(scale(AMT) ~ scale(coor.clim.raw.ch$Lat))

##### Getting the residuals for temperature
res.amt.lat <- residuals(amt.ind)

##### Adding the residuals to the dataframe
bayes.amp.all <- cbind(bayes.amp.all, res.amt.lat)

##### Getting the average of AMT absolute residual per sampling locations
amt.res <- as.data.frame(summarize(abs(bayes.amp.all$res.amt.lat), by=bayes.amp.all$Pop, FUN=mean))[,2]

colnames(bayes.amp.all)

##### Getting dataframe with residuals and orders
amt.res <- cbind(ord.pop, amt.res)

##### Ordering according to sampling locations appearance (ordering coming from Stack) as in BYAESCENV
amt.res <- amt.res[order(amt.res$appear, decreasing=FALSE),]

##### Putting in one single row as required
AMTp.res <- t(amt.res$amt.res) 

##### Write the vector to a .txt file for bayescenv 
write.table(AMTp.res, "./runs/AMTp.res.txt", row.names = FALSE, col.names = FALSE)


#### Checking Convergence and q-values (Temperature)

##### I ran bayescenv in Ubuntu. Next, I uploaded the outcome files

##### After running bayescenv-1.1 check for convergence

##### Read the .sel file where the values for the iterations are
chain1 <- read.table("./runs/run.amt.res1.sel" , header = TRUE)
chain1 <- chain1[-c(1)]

##### Adapt thin to its actual value
chain1 <- mcmc(chain1, thin = 10)

##### To plot the "trace" and the posterior distribution
plot(chain1)

#View(chain)  
summary(chain1)

##### To look for auto-correlation
autocorr.diag(chain1) 

##### To compute effective sample size
effectiveSize(chain1)    

##### Testing convergence with Geweke's convergence diagnostic
geweke.diag(chain1, frac1=0.1, frac2=0.5)

##### To test for convergence with Heidelberg and Welch's convergence diagnostic
heidel.diag(chain1, eps = 0.1, pvalue = 0.01)

##### checking output
run1.fst <- read.table("./runs/run.amt.res1_fst.txt", header=TRUE)
View(run1.fst)

##### Check the lowest prob of getting a right model and g parameter
min(run1.fst$PEP_g)



#### RUN BAYESCENV WITH RESIDUALS OF AMP VS LATITUDE FOR INDEPENDENT EFFECTS OF PRECIPITATION  ####


#### Preparing Environmental variables


##### Run a lineal regression between AMP and Latitude
amp.ind <- lm(scale(AMP) ~ scale(coor.clim.raw.ch$Lat))

##### Getting the residuals for precipitation
res.amp.lat <- residuals(amp.ind)

##### Adding the residuals to the dataframe
bayes.amp.all <- cbind(bayes.amp.all, res.amp.lat)

##### Getting the average of AMP absolute residual per sampling locations
amp.res <- as.data.frame(summarize(abs(bayes.amp.all$res.amp.lat), by=bayes.amp.all$Pop, FUN=mean))[,2]

##### Getting dataframe with residuals and orders
amp.res <- cbind(ord.pop, amp.res)

##### Ordering according to sampling locations appearance (ordering coming from Stack) as in BYAESCENV
amp.res <- amp.res[order(amp.res$appear, decreasing=FALSE),]

##### Putting in one single row as required
AMPp.res <- t(amp.res$amp.res) 

##### Write the vector to a .txt file for bayescenv 
write.table(AMPp.res, "./runs/AMPp.res.txt", row.names = FALSE, col.names = FALSE)


#### Checking convergence and q-values

##### I ran bayescenv in Ubuntu. After that I just upload the outcome files

##### After running bayescenv-1.1 check for convergence

##### Read the .sel file where the values for the iterations are
chain2 <- read.table("./runs/run.amp.res1.sel" , header = TRUE)
chain2 <- chain2[-c(1)]

##### Adapt thin to its actual value
chain2 <- mcmc(chain2, thin = 10)

##### To plot the "trace" and the posterior distribution
plot(chain2)

#View(chain)  
summary(chain2)

##### To look for auto-correlation
autocorr.diag(chain2) 

##### To compute effective sample size
effectiveSize(chain2)    

##### Testing convergence with Geweke's convergence diagnostic
geweke.diag(chain2, frac1=0.1, frac2=0.5)

##### To test for convergence with Heidelberg and Welch's convergence diagnostic
heidel.diag(chain2, eps = 0.1, pvalue = 0.01)

##### checking output
run2.fst <- read.table("./runs/run.amp.res1_fst.txt", header=TRUE)
View(run2.fst)

##### Check the lowest prob of getting a right model and g parameter
min(run2.fst$PEP_g)



#### RUN BAYESCENV WITH RESIDUALS OF PS VS LATITUDE FOR INDEPENDENT EFFECTS OF PREC SEASONALITY  ####


#### Preparing Environmental variables


##### Run a lineal regression between PS and Latitude
ps.ind <- lm(scale(PS) ~ scale(coor.clim.raw.ch$Lat))

##### Getting the residuals for precipitation seasonality
res.ps.lat <- residuals(ps.ind)

##### Adding the residuals to the dataframe
bayes.amp.all <- cbind(bayes.amp.all, res.ps.lat)

##### Getting the average of PS absolute residual per sampling locations
ps.res <- as.data.frame(summarize(abs(bayes.amp.all$res.ps.lat), by=bayes.amp.all$Pop, FUN=mean))[,2]

##### Getting dataframe with residuals and orders
ps.res <- cbind(ord.pop, ps.res)

##### Ordering according to sampling locations appearance (ordering coming from Stack) as in BYAESCENV
ps.res <- ps.res[order(ps.res$appear, decreasing=FALSE),]

##### Putting in one single row as required
PSp.res <- t(ps.res$ps.res)

##### Write the vector to a .txt file for bayescenv 
write.table(PSp.res, "./runs/PSp.res.txt", row.names = FALSE, col.names = FALSE)


#### Checking convergence and q-values

##### I ran bayescenv in Ubuntu. After that I just upload the outcome files

##### After running bayescenv-1.1 check for convergence

##### Read the .sel file where the values for the iterations are
chain3 <- read.table("./runs/run.ps.res1.sel" , header = TRUE)
chain3 <- chain3[-c(1)]

##### Adapt thin to its actual value
chain3 <- mcmc(chain3, thin = 10)

##### To plot the "trace" and the posterior distribution
plot(chain3)

#View(chain)  
summary(chain3)

##### To look for auto-correlation
autocorr.diag(chain3) 

##### To compute effective sample size
effectiveSize(chain3)    

##### Testing convergence with Geweke's convergence diagnostic
geweke.diag(chain3, frac1=0.1, frac2=0.5)

##### To test for convergence with Heidelberg and Welch's convergence diagnostic
heidel.diag(chain3, eps = 0.1, pvalue = 0.01)

##### checking output
run3.fst <- read.table("./runs/run.ps.res1_fst.txt", header=TRUE)
View(run3.fst)

##### Check the lowest prob of getting a right model and g parameter
min(run3.fst$PEP_g)



#### RUN BAYESCENV WITH RESIDUALS OF LATITUDE VS AMP FOR INDEPENDENT EFFECTS OF LATITUDE  ####


#### Preparing geographic variables


##### Run a lineal regression between Latitude and precipitation
Lat.ind <- lm(scale(coor.clim.raw.ch$Lat) ~ scale(AMP))

##### Getting the residuals for Latitude
res.lat.amp <- residuals(Lat.ind)

##### Adding the residuals to the dataframe
bayes.amp.all <- cbind(bayes.amp.all, res.lat.amp)

##### Getting the average of Latitude absolute residual per sampling locations
lat.res <- as.data.frame(summarize(abs(bayes.amp.all$res.lat.amp), by=bayes.amp.all$Pop, FUN=mean))[,2]

##### Getting dataframe with residuals and orders
lat.res <- cbind(ord.pop, lat.res)

##### Ordering according to sampling locations appearance (ordering coming from Stack) as in BYAESCENV
lat.res <- lat.res[order(lat.res$appear, decreasing=FALSE),]

##### Putting in one single row as required
Latp.res <- t(lat.res$lat.res)

##### Write the vector to a .txt file for bayescenv 
write.table(Latp.res, "./runs/Latp.res.txt", row.names = FALSE, col.names = FALSE)


#### Checking convergence and q-values


##### I ran bayescenv in Ubuntu. After that I just upload the outcome files

##### After running bayescenv-1.1 check for convergence

##### Read the .sel file where the values for the iterations are
chain4 <- read.table("./runs/run4.sel" , header = TRUE)
chain4 <- chain4[-c(1)]

##### Adapt thin to its actual value
chain4 <- mcmc(chain4, thin = 10)

##### To plot the "trace" and the posterior distribution
plot(chain4)

#View(chain)  
summary(chain4)

##### To look for auto-correlation
autocorr.diag(chain4) 

##### To compute effective sample size
effectiveSize(chain4)    

##### Testing convergence with Geweke's convergence diagnostic
geweke.diag(chain4, frac1=0.1, frac2=0.5)

##### To test for convergence with Heidelberg and Welch's convergence diagnostic
heidel.diag(chain4, eps = 0.1, pvalue = 0.01)

##### checking output
run4.fst <- read.table("./runs/run.lat.res1_fst.txt", header=TRUE)
View(run4.fst)

##### Check the lowest prob of getting a right model and g parameter
min(run4.fst$PEP_g)

#### RUN BAYESCENV WITH RESIDUALS OF CMI VS LATITUDE FOR INDEPENDENT EFFECTS OF CMI  ####


#### Preparing Environmental variables


##### Run a lineal regression between ClimaticMositureIndex (cmi) and latitude
cmi.ind <- lm(scale(cmi.all) ~ scale(coor.clim.raw.ch$Lat))

##### Getting the residuals for cmi
res.cmi.lat <- residuals(cmi.ind)

##### Adding the residuals to the dataframe
bayes.amp.all <- cbind(bayes.amp.all, res.cmi.lat)

##### Getting the average of cmi absolute residual per sampling locations
cmi.res <- as.data.frame(summarize(abs(bayes.amp.all$res.cmi.lat), by=bayes.amp.all$Pop, FUN=mean))[,2]

##### Getting dataframe with residuals and orders
cmi.res <- cbind(ord.pop, cmi.res)

##### Ordering according to sampling locations appearance (ordering coming from Stack) as in BYAESCENV
cmi.res <- cmi.res[order(cmi.res$appear, decreasing=FALSE),]

##### Putting in one single row as required
CMIp.res <- t(cmi.res$cmi.res)

##### Write the vector to a .txt file for bayescenv 
write.table(CMIp.res, "./runs/CMIp.res.txt", row.names = FALSE, col.names = FALSE)


#### Checking convergence and q-values


##### I ran bayescenv in Ubuntu. After that I just upload the outcome files

##### After running bayescenv-1.1 check for convergence

##### Read the .sel file where the values for the iterations are
chain5 <- read.table("./runs/run5.sel" , header = TRUE)
chain5 <- chain5[-c(1)]

##### Adapt thin to its actual value
chain5 <- mcmc(chain5, thin = 10)

##### To plot the "trace" and the posterior distribution
plot(chain5)

#View(chain)  
summary(chain5)

##### To look for auto-correlation
autocorr.diag(chain5) 

##### To compute effective sample size
effectiveSize(chain5)    

##### Testing convergence with Geweke's convergence diagnostic
geweke.diag(chain5, frac1=0.1, frac2=0.5)

##### To test for convergence with Heidelberg and Welch's convergence diagnostic
heidel.diag(chain5, eps = 0.1, pvalue = 0.01)

##### checking output
run5.fst <- read.table("./runs/run5_fst.txt", header=TRUE)
View(run5.fst)

##### Check the lowest prob of getting a right model and g parameter
min(run5.fst$PEP_g)


#### TESTING  ####


##### Read the .sel file where the values for the iterations are
chain.test <- read.table("./test/test1.sel" , header = TRUE)
chain.test <- chain.test[-c(1)]

##### Adapt thin to its actual value
chain.test <- mcmc(chain.test, thin = 5)

##### To plot the "trace" and the posterior distribution
plot(chain.test)

#View(chain)  
summary(chain.test)

##### To look for auto-correlation
autocorr.diag(chain.test) 

##### To compute effective sample size
effectiveSize(chain.test)    

##### Testing convergence with Geweke's convergence diagnostic
geweke.diag(chain.test, frac1=0.1, frac2=0.5)

##### To test for convergence with Heidelberg and Welch's convergence diagnostic
heidel.diag(chain.test, eps = 0.1, pvalue = 0.01)

##### checking output
test.fst <- read.table("./test/test1_fst.txt", header=TRUE)
View(test.fst)

##### Check the lowest prob of getting a right model and g parameter
min(test.fst$PEP_g)




################################################################################
################        Comparison MSOS vs BAYESCENV            ################



##### Getting the SNPs outliers found by MSOD
msr.val <- melt(do.call(cbind, list(AMT=b.AMT, AMP=b.AMP, PS=b.PS, lat=b.Y)))
colnames(msr.val) <- c("snp", "var", "pvalue")
msr.val <- msr.val[msr.val$pvalue <= 0.01,]

##### Number of snp under selection with MSOD
length(unique(msr.val$snp))

##### Count the number of SNPs per variable
as.data.frame(table(msr.val$var))

##### Count the number of times the SNP are present across variables
snp.msr.count <- as.data.frame(table(msr.val$snp))[order(as.data.frame(table(msr.val$snp))$Freq, decreasing=TRUE),]

##### Check what variables are associated to this SNP
msr.val[msr.val$snp == 123,] 

##### Number of SNPs per frequency 
table(snp.msr.count$Freq)



##### Getting the SNP outliers found by BAYESCENV
bayesc.val <- do.call(rbind, list(AMT=run1.fst, AMP=run2.fst, PS=run3.fst, Lat=run4.fst, CMI=run5.fst))
bayes.out <- bayesc.val[bayesc.val$qval_g <= 0.01,]

##### Name of the SNP outliers
snp.bayes <- do.call(rbind.data.frame, strsplit(as.vector(row.names(bayes.out)), "[.]"))
colnames(snp.bayes) <- c("var", "snp")

length(unique(snp.bayes$snp))

##### Count the number of times the SNP are present across variables
snp.bay.count <- as.data.frame(table(snp.bayes$snp))[order(as.data.frame(table(snp.bayes$snp))$Freq, decreasing=TRUE),]


##### Check how many SNPs are shared by techniques
merge(msr.val, snp.bayes, by="snp")

## Get the SNP outlier for both techniques : For MSR the ones that associate to 4 variables)
snp.out <- c(as.vector(snp.msr.count[snp.msr.count$Freq == 4,1]), as.vector(snp.bayes$snp))
length(snp.out)

##### Get the SNP outlier for both techniques 
snp.out <- unique(c(as.vector(snp.msr.count$Var1), as.vector(snp.bay.count$Var1)))
length(x = snp.out)




################################################################################
################           BLASTING SNPs OUTLIERS               ################
## Select row in chrm that correspond to the loci in sig.snp[,1]
## This is useful to BLAST the sequences in populations.loci.fa from the same run of STACKS
## The first column Chrm is the same number in populations.loci.fa for loci name
chrm.sig <- chrm[chrm$loci %in% snp.out,]
chrm.sig$CHROM <- paste(">CLocus_", chrm.sig$CHROM, sep="")
View(chrm.sig)


##### Reding the whole sequences of fragments that contains the SNP outliers
loci.fa.path <- "C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/FastQC_Stacks/denovo_n5nb_sh11/populations.loci.fa"

##### Select loci in loci.fa.patch fasta file that match the SNP outliers using the columns   vcsa CHROM that matches the loci name
loci.fa <- read.table(loci.fa.path, header=FALSE, sep = "\n")

##### Get index of the loci that match the SNP outliers and the line folling the match
loci.index <- which(loci.fa$V1 %in% chrm.sig$CHROM)
loci.index <- c(loci.index, loci.index+1)

##### Order the index
loci.index <- sort(loci.index)

##### Get the fasta sequences of the SNP outliers
loci.fa.sig <- loci.fa[loci.index,]

##### Write the fasta file with the sequences of the SNP outliers
write.table(loci.fa.sig, file = "C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/sig.snps.fa", quote = FALSE, row.names = FALSE, col.names = FALSE)

View(loci.index)
length(loci.fa)

##### Write the fasta file with the sequences of the SNP outliers
write.table(loci.fa, file = "C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/Clines.fa", quote = FALSE, row.names = FALSE, col.names = FALSE)


#################################################################################################
################       DEPRECATED: Get the predicted genes from BLAST            ################

##### Read results of BLAST
blast <- read.table("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/GO/blastn_wren_taegut_jun28_online_9TMRJJ93016-Alignment.txt", header=FALSE, sep = "\n")
dim(blast)

library(stringr)

## Extract the number of chromosomes in blast, which is the number after the word chromosome in the lines that start with ">Taeniopygia guttata isolate Blue55 chromosome"
chrm.blast <- str_extract_all(blast$V1, ">Taeniopygia guttata isolate Blue55 chromosome (.*?)\\s", simplify = TRUE)

## Extract everything after the string chromosome
chrm.blast <- str_extract_all(chrm.blast, "chromosome (.*?)\\s", simplify = TRUE)

## Get only the strings in chrm.blast
chrm.blast <- as.data.frame(chrm.blast[grepl("chromosome", chrm.blast)])

## remove commas and spaces in chrm.blast
chrm.blast <- as.data.frame(gsub(",", "", chrm.blast$chrm.blast))
colnames(chrm.blast) <- c("V1")

## Remove strings in crhm.blast
chrm.blast <- as.data.frame(gsub("chromosome ", "", chrm.blast$V1))
colnames(range.blast) <- c("chrm")


## Substract everything, text and numbers after the string Range 1: in blast
range.blast <- as.data.frame(blast[grepl("Range 1:", blast$V1),])
colnames(range.blast) <- c("V1")

## Remove strings in range.blast
range.blast <- as.data.frame(gsub("Range 1: ", "", range.blast$V1))
colnames(range.blast) <- c("V1")

## Split range.blast into a data.frame of two columns by the string " to "
range.blast <- as.data.frame(str_split_fixed(range.blast$V1, " to ", 2))

## Combine chr.blast and range.blast
chrm.pos.blast <- cbind(chrm.blast, range.blast)

## Columns names in chrm.pos.blast
colnames(chrm.pos.blast) <- c("chrm", "start", "end")

## Remove space in chrm.pos.blast$chrm
chrm.pos.blast$chrm <- gsub(" ", "", chrm.pos.blast$chrm)

## Making them regions
chrm.pos.reg <- paste(chrm.pos.blast$chrm, chrm.pos.blast$start, chrm.pos.blast$end, sep = ":")

View(chrm.pos.reg)

################################################################################
################                 Getting GO_id                  ################


## Read csv of Hit table from BLASTN
path.blastn <- "C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/GO/blastn_wren_taegut_jun28_online_9T2HKPE4013-Alignment-HitTable.csv"
blastn <- read.csv(path.blastn, header=TRUE, sep = ",")

View(blastn)

## Installing Biomart
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)

## Connect to the desired BioMart database. The biomart argument should be given a valid name
ensembl <- biomaRt::useEnsembl(biomart = "genes")

## Check the dataset names
datasets <- biomaRt::listDatasets(ensembl)
View(datasets)

## Another way to check the dataset name more specifically
biomaRt::searchDatasets(mart = ensembl, pattern = "tguttata")

## Select the dataset
mart <- biomaRt::useEnsembl(biomart = "ensembl",
                  dataset = "tguttata_gene_ensembl")


gen.id <- biomaRt::getBM(attributes = c("entrezgene_id", "ensembl_gene_id", 
                              "external_gene_name", 'go_id', 'hgnc_symbol', 
                              "chromosome_name", "start_position","end_position"),
                        filters = c('chromosomal_region'),
                        values = chrm.pos.reg,
                        mart = mart)

                        View(chrm.pos.reg)

## number of genes associated with our 17 SNP outliers
length(unique(gen.id$go_id))

## Checking results
View(head(gen.id))
dim(gen.id)
colnames(gen.id)
length(unique(gen.id$hgnc_symbol))
View(gen.id)

## Check if ensembl_gene_id has NA
sum(is.na(gen.id$ensembl_gene_id))


################################################################################

BiocManager::install("GO.db")
library(GO.db)

## Checking
columns(GO.db)

## Select uniques genes in gen.id
gene.name <- unique(gen.id$go_id)

## GEt the Gene description and category of gene ontology
go.out <- AnnotationDbi::select(GO.db, keys=gene.name, keytype="GOID", columns=c("TERM","ONTOLOGY") )

View(go.out)

################################################################################


## Get the Galllus gallus data base
BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)

## Enrichment Analysis
wrengo <- clusterProfiler::enrichGO(gene          = gene.name,
                                    keyType = "SYMBOL",
                                    OrgDb         = "org.Gg.eg.db",
                                    ont           = "ALL",
                                    pAdjustMethod = "bonferroni",
                                    pvalueCutoff  = 0.05,
                                    #qvalueCutoff  = 0.1,
                                    readable      = TRUE)


## Check results
head(wrengo)
class(wrengo)
View(wrengo)

## Converting to data.frame
wrengo.dt <- as.data.frame(wrengo)

## Checking reuslts in the dataframe
View(wrengo.dt)
length(unique(wrengo.dt$Description))
length(unique(wrengo.dt$ID))

## Filter wrengo.dt where Ontology is BP
wrengo.bp <- wrengo.dt[wrengo.dt$ONTOLOGY == "BP",]
View(wrengo.bp)
length(unique(wrengo.bp$Description))
length(unique(wrengo.bp$ID))

## Select row in wrengo.bp with the lowest p.adjust value
View(wrengo.bp[wrengo.bp$p.adjust == min(wrengo.bp$p.adjust),])

## Filter wrengo.dt where Ontology is MF
wrengo.mf <- wrengo.dt[wrengo.dt$ONTOLOGY == "MF",]
View(wrengo.mf)
length(unique(wrengo.mf$Description))
length(unique(wrengo.mf$ID))

## Select row in wrengo.mf with the lowest p.adjust value
View(wrengo.mf[wrengo.mf$p.adjust == min(wrengo.mf$p.adjust),])


## Filter wrengo.dt where Ontology is CC
wrengo.cc <- wrengo.dt[wrengo.dt$ONTOLOGY == "CC",]
View(wrengo.cc)
length(unique(wrengo.cc$Description))
length(unique(wrengo.cc$ID))

## Select row in wrengo.cc with the lowest p.adjust value
View(wrengo.cc[wrengo.cc$p.adjust == min(wrengo.cc$p.adjust),])


## Plotting the results
BiocManager::install("enrichplot")
library(enrichplot)

## dotplot
dotplot(wrengo, showCategory=20)

??enrichDAVID

################################################################################


## Installing biomartr
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomartr")

library(biomartr)

## Getting the gene ontology
GO_tbl <- getGO(organism="Taeniopygia guttata", filter="external_gene_name", 
                genes=unique(gen.id$external_gene_name), 
                verbose=TRUE)

dim(GO_tbl)
head(GO_tbl)
unique(GO_tbl$goslim_goa_description)


################################################################################


#### SAVING  ####
save.image("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/Clines.RData")
