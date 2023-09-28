##################################################################################
############           Estimation of Diversity Indices                ############ 


rm(list=ls())


# Setting the working directory
setwd("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/R/Genepop")


# Load the outcomes so far
load("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/R/Genepop/Fst_diversity.RData")


##################################################################################
############                  Installing Packages                     ############
hist(rnorm(200))

# Installing necessary packages
install.packages("adegenet")
install.packages("genepop")
install.packages("hierfstat")
install.packages("poppr")
install.packages("devtools", repos="https://cran.r-project.org/web/packages/devtools/index.html", type="source")
install.packages("BiocManager")
install.packages("popkin")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("doBy")
install.packages("doParallel")
library("devtools")
install_github("zhengxwen/SNPRelate",force = TRUE)
install_github("green-striped-gecko/dartR")

library(adegenet)
library(genepop)
library(hierfstat)
library(poppr)
library(devtools)
library(BiocManager)
library(popkin)
library(ggplot2)
library(reshape2)
library(doBy)
library(SNPRelate)
library(dartR)
library(dplyr)
library(parallel)
library(doParallel)



no_cores<- detectCores(logical=TRUE)
cl <- makeCluster(no_cores)
registerDoParallel(cl)



##################################################################################
############                  Data Formatting                         ############



# Reading and converting STRUCTURE file into genind
wren.n5nb.genind <- import2genind("populations2.str", onerowperind=FALSE, n.ind=112, n.loc=4409, 
                                  col.lab=1, col.pop=2, NA.char="0", ask=FALSE)

#### I converted a VCF file to Genepop using PGDspider. 
#### Fist, I got the meta data of sampling locations

# Meta data with sampling locations and observational
# classification of species, subspecies and hybrids
indtable <- read.table("wren_pop.txt", header = TRUE, sep="\t")

# Latitude as numeric
indtable$Ord_Lat <- as.numeric(indtable$Ord_Lat)

# Order of samples as in genomic data
# Order used in assembly in Stack
stack.ord <- read.table("popmap_wren_n5nb_sh11.txt", col.names = c("Ind", "Ord_stack"))

# Adding a column for the order of samples
stack.ord$Ord_stack <- 1:nrow(stack.ord)

# Metadata of sampling location and Stack order
pop.lab <- merge(indtable, stack.ord, by.x="Lev6_ind", by.y="Ind")

# Giving the same order as genomic data
pop.lab <- pop.lab[order(pop.lab$Ord_stack, decreasing = FALSE),]

# Adding code to breeding group for Las Golondrinas
pop.lab$Lev5_grp[pop.lab$Lev6_ind == "1836"] <- "CZ080"
pop.lab$Lev5_grp[pop.lab$Lev6_ind == "1835"] <- "CZ080"

# Read the text file with the structure assignment probabilities
struc <- read.table("k4popfile.txt", header = TRUE, sep="\t")

# Merging the structure assignments and population metadata
struc.pop <- merge(pop.lab, struc, by.x="Lev6_ind", by.y="sample")

# Write metadata for PGDspider for later analyses
write.table(pop.lab, "poplab.txt", row.names = FALSE, col.names = FALSE)

# Filtering the genind object with metadata
# To take out C. brunneicapillus samples
wren.n5nb.genind <- wren.n5nb.genind[indNames(wren.n5nb.genind) %in% pop.lab$Lev6_ind]

# Setting up the populations in the genind data
strata(wren.n5nb.genind) <- struc.pop[,c(1:6,32)]
setPop(wren.n5nb.genind) <- ~top1name

# Getting a vector with the sites in order of Latitude
site.lat.order2 <- unique(indtable[order(indtable$Ord_Lat),]$Lev4_pop)[1:16]


##################################################################################
############                Data for Structure2                       ############

# this is a data formatting to put the identity of location to a structure format
# file so i can run a model with LOCPRIOR = 1. 

# Read the text file with different designation of populations
data.raw <- read.table("populations_n5nbx.txt", header = TRUE)

# Select the info of population
pop.factor <- pop.lab[,c('Lev6_ind', "Lev4_pop", "Ord_stack")]

# Adding a id for popualtion in numeric format
pop.factor$pop <- as.numeric(as.factor(pop.lab$Lev4_pop))

# Getting the file with all information
data.pop.str2 <- merge (pop.factor, data.raw, by="Lev6_ind")

# Ordering according to the Stack order
data.pop.str2 <- data.pop.str2[order(data.pop.str2$Ord_stack.x, decreasing=FALSE),]

# Removing the columns that are not needed
data.pop.str2 <- data.pop.str2[,-c(2,3,5,6)]  

write.table(data.pop.str2, "population_n5nb_str2.txt", row.names = FALSE, col.names = TRUE)


##################################################################################
############      Testing Hardy-Weinberg Equilibrium (Genepop)        ############

# Got Genepop format in PDGSpider

# Performing the test for HWE
hwe.wn5nb <- genepop::test_HW("wn5nb.genepopPGDspider.txt", which = "Proba", verbose=interactive())
# Results storaged in the working directory as "wn5nb.genepopPGDspider.txt.P"
# If all P-values across loci in each population are >0.05, 
# do not reject the null hypothesis that all loci are in HWE.

# Taking a look of the results
View(read.table("wn5nb.genepopPGDspider.txt.P", header=FALSE, sep="\t"))

# I modified the outcome "filename.txt.P" to have just sites and prob by loci
# Reading the file outcome from genepop
HWE.genepop <- read.table("wn5nb.genepop_HWE.txt", header=TRUE, sep="\t")

# Getting loci with less than P-val 0.05 and failed to be in HWE
HWE.genepop1 <- subset(HWE.genepop, P.val>0.05)


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
gl.hwe.ld <- gl.report.ld(gl.hwe, ncores=18)

# Getting pairwise loci with r2>0.1 with strong linked loci
gl.LD <- subset(gl.hwe.ld,R2>0.5)
View(gl.LD)

# Percentage of pairwise loci linked
(length(gl.LD$loc1)/length(gl.hwe.ld$loc1))*100


### Brief explanation of LD outcome

# Values of D close to zero (D=0) indicates LE, random combination of alleles
# Values of D far from zero (D???0) indicates LD, alleles are correlated
# Values of r2 close to zero (r2=0) indicates LE, random combination of alleles
# Values of r2 far from zero (r2???0) indicates LD, alleles are correlated
# To give you a quick and dirty answer, a D? of 0.8 is high disequilibrium. 
# Basically the two SNPs are coinherited roughly 80% of the time. 
# The reason your r2 is low is that this takes account of allele frequency.
# D? and r2 values are widely used but poorly understood. 
# The current "trend" seems to be to take more notice of the r2 value whereas 
# I feel that the D? is more meaningful and easier to understand. 
# The idea of disequilibrium values is that they are a measure of the non-random 
# association of alleles at two or more loci, i.e how often alleles are coinherited. 
# If two loci are not coinherited at all (they are independent) then both the D? 
# and r2 values will be 0.0 irrespective of either allele frequency. 
# As another example if you had two polymorphisms both with a 50% allele frequency and in 
# total disequilibrium then both the D? and r2 values would be 1.0. 
# However the story changes when the allele frequencies are not the same. 
# For example if you had two polymorphisms, one with a 50% allele frequency and the other with 
# a 1% allele frequency that were still in total disequilibrium then the D? 
# value would be 1.0 but the r2 value would only be 0.01. Basically the D? 
# is saying when the rare allele is present it is always inherited with 
# one particular allele of the 50% polymorphism whereas the r2 is saying it is a rare allele 
# so the vast majority of the time the common allele is not found with it 
# (but only because it is rare, not because it is not in disequilibrium). 
# However, even with such a low r2 this SNP adds nothing to an association study 
# because they are in complete disequilibrium.
# a better understanding with another example where the two SNPs are coinherited 
# about half the time. With both SNPs having a 50% allele frequency the D? value would be 
# about 0.5 and the r2 value would be about 0.25 but if one SNP had a 1% allele frequency 
# then the D? value would still be about 0.5 but the r2 value would only be about 0.005






##################################################################################
############       Basic genetic statistics NO bootstrapping          ############

# Estimating General Genetics Statistics 
gl.basics <- gl.basic.stats(gl.hwe)

# Changing format for plotting
gl.Fis <- as.data.frame(melt(gl.basics$Fis))
colnames(gl.Fis)<- c('loci', 'pop', 'Fis')

# Making graph for outcome from gl.basics.stats
# Fis per locus and Geopop, CI show the variation of Fis across locus with geopop
plot.fis.all <- ggplot(gl.Fis, aes(x=factor(pop, level=site.lat.order2), y=Fis)) +
                geom_boxplot() +
                ggtitle("Fis with all individuals")+
                theme(plot.title = element_text(hjust = 0.5))+
                geom_boxplot(notch = TRUE, fill = "lightgray")+
                stat_summary(fun = mean, geom = "point",
                shape = 18, size = 2.5, color = "#FC4E07") +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plotting Heterozygosity
# Observed
# Mean of Het across locus within geopop 
Ho.all <- as.data.frame(colMeans(gl.basics$Ho, na.rm=TRUE))
Ho.all$pop <- rownames(Ho.all)
colnames(Ho.all)[1] <- "value"

# Gene Diversity
Hs.all <- as.data.frame(colMeans(gl.basics$Hs, na.rm=TRUE))
Hs.all$pop <- rownames(Hs.all)
colnames(Hs.all)[1] <- "value"

# Get both together
Het.all <- rbind(Ho.all, Hs.all)
Het.all$stat <- c(rep("Ho",16), rep("Hs", 16))


plot.het.all <- ggplot(Het.all, aes(fill=stat, x=factor(pop, level=site.lat.order2), y=value))+
                geom_bar(position="dodge", stat="identity")+
                ggtitle("Heterozygosity and Gene Diversity with All Individuals")+
                theme(plot.title = element_text(hjust = 0.5))+
                theme_bw()+
                theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
                xlab("Populations")+ylab("Value")+
                scale_fill_grey()


# Plotting the observed vs the expected heterozygosity

# Converting back to Genind object
gi.hwe <- gl2gi(gl.hwe)

pop(gi.hwe)

# Subsetting genind object per sampling location
gi.hwe.pop <- list()
for(i in 1:length(unique(pop(gi.hwe)))){
  gi.hwe.pop[[i]] <- popsub(gi.hwe, sublist=unique(pop(gi.hwe))[i])}


# Getting the summaries of genind per pop.
summary.hwe.all <- list()
for(i in 1:length(gi.hwe.pop)){
  summary.hwe.all[[i]] <- adegenet::summary(gi.hwe.pop[[i]])
}

# Making scatterplots of Hexp vs Hobs per population
par(mfrow = c(4, 4), mar=c(3,3,3,2)) 
for(i in 1:length(summary.hwe.all)){
  plot(summary.hwe.all[[i]]$Hobs, summary.hwe.all[[i]]$Hexp, pch=20, cex=3, xlim=c(0,1), ylim=c(0,1))
  abline(0,1,lty=2)}

# Getting Ho in a table
Ho.mean.sd <- data.frame(pop=as.character(unique(pop(gi.hwe))), 
                         mean=rep(0, 16), sd=rep(0, 16))
for(i in 1:length(summary.hwe.all)){
  Ho.mean.sd$mean[i] <- mean(summary.hwe.all[[i]]$Hobs)
  Ho.mean.sd$sd[i] <- sd(summary.hwe.all[[i]]$Hobs)
}

# Getting He in a table
He.mean.sd <- data.frame(pop=as.character(unique(pop(gi.hwe))), 
                         mean=rep(0, 16), sd=rep(0, 16))
for(i in 1:length(summary.hwe.all)){
  He.mean.sd$mean[i] <- mean(summary.hwe.all[[i]]$Hexp)
  He.mean.sd$sd[i] <- sd(summary.hwe.all[[i]]$Hexp)
}



##################################################################################
############        Basic genetic statistics one individual           ############


# Performing basic genetic statistics for one subsample 
# With one individual per breeding group and 4 ind per pop


# Looking at the number of individuals per pop
table(pop(gi.hwe))
# The site Manglareschurute has less than 4 individuals 

# Checking number of breeding groups per geopop
grp.pop <- pop.lab %>%
  group_by(Lev4_pop) %>%
  summarise(n_distinct(Lev5_grp))

# Geopop with more than 3 breeding groups
grp.pop4 <- subset(grp.pop, `n_distinct(Lev5_grp)`>3)$Lev4_pop

# Subsetting pop.lab by the pop with more than 4 group
pop.lab.4 <- pop.lab[pop.lab$Lev4_pop %in% grp.pop4, ]

# Getting the subsample of individuals per breeding group
pop.lab.subsample <- do.call(rbind, lapply(split(pop.lab.4, pop.lab.4$Lev5_grp), 
                                           function(x) x[sample(nrow(x), 1), ]))

# Subsetting the data with population with four or more individuals
gi.hwe.4 <- popsub(gi.hwe, sublist=grp.pop4)

# Converting genind to hier object
hier.hwe <- genind2hierfstat(gi.hwe.4)

# Reduce the hier data to get the Fis 
hier.sample <- as.data.frame(hier.hwe[rownames(hier.hwe) %in% pop.lab.subsample$Lev6_ind,])

# Getting the basic genetic stats
hier.sample.basic <- basic.stats(hier.sample)

# Getting the Fis per loci per pop
pp.sample.Fis <- hier.sample.basic$Fis

# Making a boxplot of Fis per population with mean per loci
par(mfrow=c(1,1), mar=c(8,4,1,1))
boxplot(pp.sample.Fis, col=colr, las=3, ylab="Fis")
mtext(text="Populations",
      side=1, line=6)
points(colMeans(pp.sample.Fis, na.rm=TRUE), col="black") # showing also the mean


# Converting to dosage data
dos.sample <- as.data.frame(fstat2dos(hier.sample[,-c(1)], diploid=TRUE))

# Estimating Fis per pop
fs.dos.sample <- fs.dosage(dos.sample, hier.sample[,1])

# Estimating bootstrapping for Fis
Fis.sample.boot <- boot.ppfis(dat=hier.sample,nboot=100000,quant=c(0.025,0.975),diploid=TRUE,dig=4)

# Getting the mean of Fis per pop
sample.fis.mean <- fs.dos.sample$Fs[1,1:ncol(fs.dos.sample$Fs)-1]

sample.fis.ci <- Fis.sample.boot$fis.ci
sample.fis.ci$pop <- colnames(pp.sample.Fis)
sample.fis.ci$Fis <- sample.fis.mean
site.lat.order3 <- site.lat.order2[site.lat.order2 %in% colnames(pp.sample.Fis)]

# Plotting mean and CI of Fis per population
ggplot(sample.fis.ci, aes(x=factor(pop, level=site.lat.order2), Fis)) +
  geom_point() +
  geom_errorbar(aes(ymin = ll, ymax = hl))


# Population specific Fsts and bootstrap CIs
Beta.fst.sample <- betas(hier.sample, nboot=100000, lim=c(0.05, 0.95), diploid=TRUE)


fst.ci.sample <- as.data.frame(t(rbind(Beta.fst.sample$betaiovl, Beta.fst.sample$ci)))
fst.ci.sample$Pop <- row.names(fst.ci.sample)
colnames(fst.ci.sample) <- c("Fst", "ll", "hl", "Pop")

ggplot(fst.ci.sample, aes(x=factor(Pop, level=site.lat.order3), y= Fst)) +
  geom_point() +
  geom_errorbar(aes(ymin = ll, ymax = hl))


# When calculated within one population, Fst is a measure of average pairwise distances 
# between pairs of individuals (haplotypes) in terms of allele frequencies. 
# For between population comparisons Fst is more intuitive and basically is 
# a "difference" between allele frequencies of two populations.




##################################################################################
############               Ho and Hs Bootstrapping                    ############ 


# Performing Ho and Hs with bootstrapping across individuals in the same
# breeding group with one individual per group with at least 4 individuals

set.seed(200)

# Number of bootstrapping
nsamples <- 50

# Converting genind to hier object
hier.hwe.all <- genind2hierfstat(gi.hwe)

# Information of individuals, sampling site, breeding group
# col 1: individual identifier
# col 2: species identifier
# col 3: sampling site
# col 4: breeding group identifier
pop.boot <- pop.lab.4[,-c(3,4,7,8)]

# Statistics we want to get
b.stat <- c("Ho", "Hs", "Fis")


################################################################################

# Function to get N replicates (bootstrapping) and get the mean across replicates
# of the selected genetic basic statistics from Hierfstat::basic.stats
# arg data has to be a hierfstat object with a group identifier (usually population)
# as first column. Row names should be unique (usually the id of the samples)
# Rest of columns are the loci.

boot.basic <- function(gen.data, nboots, stat, pop.meta, ord){
  
  require(hierfstat)
  require(stringr)
  require(parallel)
  require(doParallel)
  require(ggplot2)
  
  # Setting parallelization
  no_cores <- detectCores(logical=TRUE)
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  # Getting lists with sub-samples with one ind per group
  pop.list <- list()
  for(i in 1:nboots){
    pop.list[[i]] <- do.call(rbind, lapply(split(pop.meta, pop.meta[,4]), 
                                           function(x) x[sample(nrow(x), 1, replace = TRUE),]))}

  # Getting the population
  pop.label <- ord
  
  # Subsetting the hier object by the samples in pop.list
  hier.list <- list()
  for(i in 1:length(pop.list)){
    hier.list[[i]] <- as.data.frame(gen.data[rownames(gen.data) %in% pop.list[[i]][,1],])}
  
  # Changing the levels attributes, it may happens that the dataframe preserve the levels attribute of the original
  for(i in 1:length(hier.list)){
    hier.list[[i]][,1] <- factor(hier.list[[i]][,1])}
  
  # Exporting data
  parallel::clusterEvalQ(cl, c(library(hierfstat), library(stringr), library(ggplot2)))
  parallel::clusterExport(cl, c("hier.list"), envir = environment())
  
  # Applying Hierfstat::basic.stats to the replicates
  par.fun <- function(x){hierfstat::basic.stats(x)}  
  system.time(outcome <-  c(parLapply(cl, hier.list, fun=par.fun)))
  
  # Functions to estimate lower and upper bounds across replicates
  lower <- function(x){
    me <- qnorm(.95)*(sd(x, na.rm=TRUE)/sqrt(nboots))
    ll <- mean(x, na.rm=TRUE) - me
    return(ll)
  }
  
  upper <- function(x){
    me <- qnorm(.95)*(sd(x, na.rm=TRUE)/sqrt(nboots))
    hl <- mean(x, na.rm=TRUE) + me
    return(hl)
  }

  # Creating the list for the outcome
  Avg <- list()
  ll <- list()
  hl <- list()
  
  # Getting mean per group across replicates
  for(i in 1:length(stat)){
    L <- sapply(outcome, "[", stat[i]) 
    L <- array(unlist(L), dim=c(nrow(L[[1]]),
                                ncol(L[[1]]),
                                length(L)))
    
    # Setting the names in array
    dimnames(L)[[2]] <- colnames(outcome[[1]]$Ho)
    
    # Getting mean and CI across replicates
    Avg[[i]] <- apply(L, c(1,2), mean, na.rm=TRUE)
    ll[[i]] <- apply(L, c(1,2), lower)
    hl[[i]] <- apply(L, c(1,2), upper)
    #colnames(Avg[[i]]) <- pop.label
    #colnames(ll[[i]]) <- pop.label
    #colnames(hl[[i]]) <- pop.label
  }
  
  # Naming each element of the list
  names(Avg) <- stat
  names(ll) <- stat
  names(hl) <- stat
  
  # Getting all together
  G <- list(Avg, ll, hl)
  names(G) <- c("mean", "ll", "hl")
  
  # Makings a data frame
  loci <- t(as.data.frame(do.call(cbind, as.data.frame(G))))
  
  # Split the names to convert in a dataframe
  names <- as.data.frame(t(as.data.frame(str_split(rownames(loci), "[.]", n=3))))
  names[,3] <- gsub(".", " ", names[,3], fixed=TRUE)
  rownames(names) <- c(1:nrow(names))
  
  # Combine names and statistics in the same dataframe
  Stats <- as.data.frame(cbind(names, loci))
  
  # Row names
  rownames(Stats) <- 1:nrow(Stats)
  
  # Names of columns one to three
  colnames(Stats)[1:3] <- c("Stat", "Fstat", "Pop")
  
  stopCluster(cl)  #stop the cluster
  
  Ho.mean <- t(Stats[with(Stats, Stat== "mean" & Fstat=="Ho"),][,-c(1:2)])
  colnames(Ho.mean) <- Ho.mean[1,]
  Ho.mean <- sapply(as.data.frame(Ho.mean[-1,]), as.numeric)

  Hs.mean <- t(Stats[with(Stats, Stat== "mean" & Fstat=="Hs"),][,-c(1:2)])
  colnames(Hs.mean) <- Hs.mean[1,]
  Hs.mean <- sapply(as.data.frame(Hs.mean[-1,]), as.numeric)
  
  # Making scatterplots of Hs vs Ho per population from the basic.stats function of Hierfstat
  par(mfrow = c(2, 4), mar=c(5,5,5,2)) 
  scatterplot.4ind <- for(i in 1:ncol(Ho.mean)){
                        plot(Hs.mean[,i], Ho.mean[,i], pch=20, cex=2, xlim=c(0,1), ylim=c(0,1),
                        main=colnames(Ho.mean)[i], xlab="Hs", ylab="Ho")
                        abline(0,1,lty=2)}
  
  print(scatterplot.4ind)
  
  # Plotting Heterozygosity

  # Observed
  Ho.sample <- as.data.frame(colMeans(Ho.mean, na.rm=TRUE))
  Ho.sample$pop <- rownames(Ho.sample)
  colnames(Ho.sample)[1] <- "value"

  # Expected
  Hs.sample <- as.data.frame(colMeans(Hs.mean, na.rm=TRUE))
  Hs.sample$pop <- rownames(Hs.sample)
  colnames(Hs.sample)[1] <- "value"

  # Get both together
  Het.sample <- rbind(Ho.sample, Hs.sample)
  Het.sample$stat <- c(rep("Ho",nrow(Ho.sample)), rep("Hs", nrow(Hs.sample)))

  plot.het.4ind <-  ggplot(Het.sample, aes(fill=stat, x=factor(pop, level=ord), y=value))+
                    geom_bar(position="dodge", stat="identity")+
                    ggtitle("Heterozygosity with Four Individuals per Population")+
                    theme(plot.title = element_text(hjust = 0.5))+
                    theme_bw()+
                    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
                    xlab("Sampling Sites")+ylab("Heterozygosity")+
                    scale_fill_grey()
  
  print(plot.het.4ind)
  

  #### Plotting Heterozygosity Function
  
  Plot.het.pop <- function(x, orden, stat.het){ggplot(x, aes(x=factor(pop, level=orden), y=stat.het)) +
              geom_boxplot() +
              stat_summary(fun=mean, geom="point", shape=2, size=2)+
              theme_bw() + xlab("Sampling Sites") + 
                         ylab(colnames(x)[3])}

  #### Plotting Ho
  
  # Changing format for plotting
  Ho.mean.melted <- as.data.frame(melt(Ho.mean))
  colnames(Ho.mean.melted)<- c('loci', 'pop', 'Ho')

  # Plotting Ho per geopop and standard error across loci
  print(Plot.het.pop(Ho.mean.melted, orden=ord, stat.het=Ho.mean.melted$Ho))
  
  # Check the mean and SD per population
  Ho.sum <- summaryBy(Ho ~ pop, data=Ho.mean.melted,
            FUN=function(x) {c(m=mean(x, na.rm=TRUE), s=sd(x, na.rm=TRUE))})
  
  #### Plotting Hs

  # Changing format for plotting
  Hs.mean.melted <- as.data.frame(melt(Hs.mean))
  colnames(Hs.mean.melted)<- c('loci', 'pop', 'Hs')

  # Plotting Ho per geopop and standard error across loci
  print(Plot.het.pop(Hs.mean.melted, orden=ord, stat.het=Hs.mean.melted$Hs))
  
  # Check the mean and SD per population
  Hs.sum <- summaryBy(Hs ~ pop, data=Hs.mean.melted,
          FUN=function(x) {c(m=mean(x, na.rm=TRUE), s=sd(x, na.rm=TRUE))})

 return(list(Stats, Ho.sum, Hs.sum))
  
}

################################################################################


# Data simulation 
dat <- sim.genot(size=5,nbloc=10,nbal=2, nbpop=4)
pop <- rep(c("A","B","C","D"), 5)
dat$Pop <- as.factor(pop)
dat[7,10] <- "NA"
dat[,c(2:ncol(dat))] <- lapply(dat[,c(2:ncol(dat))], as.integer)

# Pop sim data
pop.test <- data.frame(1:20, rep(1:2, 10), pop, rep(1:10,2))

# Order we want the boxplots
ord.test <- c("D", "B", "A", "C")
Ord.7 <- c("Las Golondrinas","PVM","Patricia Pilar","Chone","Calceta","Machalilla","Zapotillo")

##### RUNNING THE FUNCTION: simulated data
basic.boot.test <- boot.basic(gen.data=dat, nboots=nsamples, stat=b.stat, 
                              pop.meta=pop.test, ord=ord.test)


##### RUNNING THE FUNCTION: real data
basic.boot <- boot.basic(gen.data=hier.hwe.all, nboots=nsamples, stat=b.stat, 
                              pop.meta=pop.boot, ord=Ord.7)


################################################################################
# NO BOOTSTRAPPING WITHIN BREEDING GROUP, WITH HWE FILTERING


#### Function for Plotting Heterozygosity and Fis

# x: A list, the outcome of hierfstats::basic.stats function.
# ord: A character vector with the order we want the sampling sites
# frame: A vector with two elements following the par(mfrow=c(nrows,ncols)) format. 
# frame: nrows x ncols >= n sampling sites

plot.basic.hier.stats <- function(x, ord, frame){
  
  require(ggplot2)
  
  ##### FUNCTIONS

  #### Scatterplot
  
  # Making scatterplots of Hs vs Ho per population from the basic.stats function of Hierfstat
  scatterplot <- function(x, orden){
    ggplot(x, aes(x=Ho, y=Hs)) +
      geom_point()+
      facet_wrap(~Var2)+
      theme_bw()+
      theme(strip.background =element_rect(fill="white"))+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
      geom_abline(intercept = 0)}
  

  #### Barplot

  # Bar plot with Hs and Ho together 
  barplot <-  function(x, orden){ggplot(x, aes(fill=stat, x=factor(pop, level=ord), y=value))+
                              geom_bar(position="dodge", stat="identity")+
                              ggtitle("Heterozygosity per Population")+
                              theme(plot.title = element_text(hjust = 0.5))+
                              theme_bw()+
                              theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
                              xlab("Sampling Sites")+ylab("Heterozygosity")+
                              scale_fill_grey()}

  #### Boxplot

  # Plotting boxplot for each type of heterozygosity 
  boxplot <- function(x, orden, stat.het){ggplot(x, aes(x=factor(pop, level=orden), y=stat.het)) +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=2, size=2)+
    theme_bw() + xlab("Sampling Sites") + 
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
    ylab(colnames(x)[3])}


  ##### DATA MANAGEMENT
  
  ##### Scatterplot
  
  sca.het <- cbind(melt(x$Ho),melt(x$Hs)[,3])
  colnames(sca.het)[c(3,4)] <- c("Ho", "Hs")
  sca.het$Var2 <- factor(sca.het$Var2, levels = site.lat.order2)
  
  ##### Barplot

  # Observed
  Ho.sample <- as.data.frame(colMeans(x$Ho, na.rm=TRUE))
  Ho.sample$pop <- rownames(Ho.sample)
  colnames(Ho.sample)[1] <- "value"

  # Expected
  Hs.sample <- as.data.frame(colMeans(x$Hs, na.rm=TRUE))
  Hs.sample$pop <- rownames(Hs.sample)
  colnames(Hs.sample)[1] <- "value"

  # Get both together
  Het.sample <- rbind(Ho.sample, Hs.sample)
  Het.sample$stat <- c(rep("Ho",nrow(Ho.sample)), rep("Hs", nrow(Hs.sample)))

  ##### Boxplot

  # Changing format for plotting
  Ho.mean.melted.all <- as.data.frame(melt(x$Ho))
  colnames(Ho.mean.melted.all)<- c('loci', 'pop', 'Ho')

  # Changing format for plotting
  Hs.mean.melted.all <- as.data.frame(melt(x$Hs))
  colnames(Hs.mean.melted.all)<- c('loci', 'pop', 'Hs')

  # Changing format for plotting
  Fis.mean.melted.all <- as.data.frame(melt(x$Fis))
  colnames(Fis.mean.melted.all)<- c('loci', 'pop', 'Fis')

  
  #### PLOTTING

  # Scatterplot
  par(mfrow = c(frame[1], frame[2]), mar=c(5,5,5,2)) 
  print(scatterplot(x=sca.het, orden=orc))

  # Barplot Ho vs Hs
  print(barplot(x=Het.sample, orden=ord))

  # Individual statistics

  # Ho
  print(boxplot(Ho.mean.melted.all, orden=ord, stat.het=Ho.mean.melted.all$Ho))

  # Hs
  print(boxplot(Hs.mean.melted.all, orden=ord, stat.het=Hs.mean.melted.all$Hs))

  # Fis
  print(boxplot(Fis.mean.melted.all, orden=ord, stat.het=Fis.mean.melted.all$Fis))

}


##### RUNNING THE FUNCTION
plot.basic.hier.stats(x=gl.basics, ord=site.lat.order2, frame=c(4,4))

################################################################################
# NO BOOTSTRAPPING WITHIN BREEDING GROUP, NO HWE FILTERING

# Converting genind object to hierfstat
hier.all <- genind2hierfstat(wren.n5nb.genind, pop=pop(wren.n5nb.genind))

# Getting the F-statistics
hier.outcome <- hierfstat::basic.stats(hier.all)

##### RUNNING THE FUNCTION
plot.basic.hier.stats(x=hier.outcome, ord=site.lat.order2, frame=c(4,4))


##################################################################################
########             Plot Hobs vs Hexp (Adegenet summaries)           ############

# Subsampling of the complete genind
# each of the elements of the list
# is a genind object with 40 individuals
# a single individual per sampling location
gi.hwe.4ind <- list()
for(j in 1:length(pop.list)){
  gi.hwe.4ind[[j]] <- gi.hwe[i=c(pop.list[[j]]$Lev6_ind)]
}

# Subsetting genind object per sampling location
# each element of the list is also a list with
# one genind object per sampling locations
# He and Hs will be estimated within location
gi.hwe.pop4ind <- list()
for(i in 1:length(gi.hwe.4ind)){
 gi.pop <- list()
 for(j in 1:length(unique(pop(gi.hwe.4ind[[i]])))){
    gi.pop[[j]] <- popsub(gi.hwe.4ind[[i]], sublist=unique(pop(gi.hwe.4ind[[i]]))[j])}
 gi.hwe.pop4ind[[i]] <- gi.pop
}

# Getting summaries for all the replicates and populations
# Get basic statistics for every replicate "i" and every sampling location "j"
gi.sum.pop4ind <- list()
for(i in 1:length(gi.hwe.pop4ind)){
  gi.pop <- list()
  for(j in 1:length(gi.hwe.pop4ind[[i]])){
    gi.pop[[j]] <- adegenet::summary(gi.hwe.pop4ind[[i]][[j]])}
  gi.sum.pop4ind[[i]] <- gi.pop}


#### GETTING HEXP

# Converting nested list to an array with just Hexp
gi.array.hexp <- list()
for(i in 1:length(gi.sum.pop4ind)){
  test <-sapply(gi.sum.pop4ind[[i]], "[", "Hexp")
  test <-sapply(test, "length<-", max(lengths(test)))
  gi.array.hexp[[i]] <- test}
gi.array.hexp <- array(unlist(gi.array.hexp), dim=c(4409,7,50))
Hexp.mean <- apply(gi.array.hexp, c(1,2), mean, na.rm=TRUE)# Taking the mean of Hexp across repetitions

# Putting the names in order
names.pop <- rep(NA,7)
for(i in 1:7){
  names.pop[i] <- names(gi.sum.pop4ind[[1]][[i]]$pop.n.all)}
colnames(Hexp.mean) <- names.pop

#### Plotting Hexp

# Changing format for plotting
Hexp.mean.melted <- as.data.frame(melt(Hexp.mean))
colnames(Hexp.mean.melted)<- c('loci', 'pop', 'value')
Hexp.mean.melted$stat <- rep("Hexp", nrow(Hexp.mean.melted))

# Plotting Hexp per geopop and standard error across loci
ggplot(Hexp.mean.melted, aes(x=factor(pop, level=site.lat.order3), y=value)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=2, size=2)+
  theme_bw() + xlab("Geographical Populations") + ylab("Expected Heterozygosity")

# Check the mean and SE per population
sum.Hexp <-summaryBy(value ~ pop, data=Hexp.mean.melted,
          FUN=function(x) {c(m=mean(x, na.rm=TRUE), s=sd(x, na.rm=TRUE)/sqrt(4409))})


#### GETTING HOBS

# Converting nested list to an array with just Hobs
gi.array.hobs <- list()
for(i in 1:length(gi.sum.pop4ind)){
  test <-sapply(gi.sum.pop4ind[[i]], "[", "Hobs")
  test <-sapply(test, "length<-", max(lengths(test)))
  gi.array.hobs[[i]] <- test}
gi.array.hobs <- array(unlist(gi.array.hobs), dim=c(4409,7,50))
Hobs.mean <- apply(gi.array.hobs, c(1,2), mean, na.rm=TRUE)# Taking the mean of Hobs across repetitions

# Putting the names in order
names.pop <- rep(NA,7)
for(i in 1:7){
  names.pop[i] <- names(gi.sum.pop4ind[[1]][[i]]$pop.n.all)}
colnames(Hobs.mean) <- names.pop


#### Plotting Hobs

# Changing format for plotting
Hobs.mean.melted <- as.data.frame(melt(Hobs.mean))
colnames(Hobs.mean.melted)<- c('loci', 'pop', 'value')
Hobs.mean.melted$stat <- rep("Hobs", nrow(Hobs.mean.melted))

# Plotting Ho per geopop and standard error across loci
ggplot(Hobs.mean.melted, aes(x=factor(pop, level=site.lat.order3), y=value)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="point", shape=2, size=2)+
  theme_bw() + xlab("Geographical Populations") + ylab("Observed Heterozygosity")

# Check the mean and SE per population
sum.Hobs <- summaryBy(value ~ pop, data=Hobs.mean.melted,
          FUN=function(x) {c(m=mean(x, na.rm=TRUE), s=sd(x, na.rm=TRUE)/sqrt(4409))})

# Get both together
Het.all.sum <- rbind(Hobs.mean.melted, Hexp.mean.melted)

# Plotting the mean of Hexp and Hobs
plot.het.sum <- as.data.frame(summaryBy(value ~ pop+stat, data=Het.all.sum,
                              FUN=function(x) {c(m=mean(x, na.rm=TRUE), s=sd(x, na.rm=TRUE))}))

plot.het.4ind.sum <- ggplot(plot.het.sum, aes(fill=stat, x=factor(pop, level=site.lat.order3), y=value.m))+
  geom_bar(position="dodge", stat="identity")+
  ggtitle("Heterozygosity with Four Individuals per Population")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  xlab("Populations")+ylab("Heterozygosity")+
  scale_fill_grey()



##################################################################################
########        Estimating Fst, Fis and Confidence intervals          ############


#### GETTING Fst AND Fis PER SAMPLING LOCATION

# Converting the element in the list to dosage format
dos.list <- list()
for(i in 1:length(hier.list)){
  dos.list[[i]] <- as.data.frame(fstat2dos(hier.list[[i]][,-c(1)], diploid=TRUE))}

# Getting Fst and Fis per populations
pop.all <- as.vector(hier.list[[1]][,1])
clusterExport(cl, "pop.all")
dos.fst.fis.fun <- function(x){
  hierfstat::fs.dosage(x, pop=pop.all)}
system.time(dos.fst.fis <- c(parLapply(cl, dos.list, fun=dos.fst.fis.fun)))

#### GETTING CI PER SAMPLING LOCATION

# Getting the confidence intervals for Fis
dos.fis.ci.fun <- function(x){
  hierfstat::boot.ppfis(dat=x,nboot=100000,quant=c(0.025,0.975),diploid=TRUE,dig=4)}
system.time(dos.fis.ci <- c(parLapply(cl, hier.list, fun=dos.fis.ci.fun)))

# Getting the confidence intervals for Fst
dos.fst.ci.fun <- function(x){
  hierfstat::betas(x, nboot=100000, lim=c(0.025, 0.975), diploid=TRUE)}
system.time(dos.fst.ci <- c(parLapply(cl, hier.list, fun=dos.fst.ci.fun)))

# Getting the co-ancestry coefficients
dos.ac.ci.fun <- function(x){
  hierfstat::betas(x, nboot=100000, lim=c(0.025, 0.975), diploid=TRUE, betaijT=TRUE)}
system.time(dos.ac.ci <- c(parLapply(cl, hier.list, fun=dos.ac.ci.fun)))


#### GETTING MEANS OF Fst AND Fis across the 50 replicates

# Getting F-stat in a list
fis.fst.list <- list()
for(i in 1:length(dos.fst.fis)){
  fis.fst.list[[i]] <- dos.fst.fis[[i]]$Fs}

# Converting the list with Fst to array
fis.fst.array <- array(unlist(fis.fst.list), dim=c(nrow(fis.fst.list[[1]]), 
                                               ncol(fis.fst.list[[1]]), 
                                               length(fis.fst.list)))

# Taking the mean across repetitions
fis.fst.mean <- apply(fis.fst.array, c(1,2), mean, na.rm=TRUE)

# Putting names to columns and rows
colnames(fis.fst.mean) <- colnames(fis.fst.list[[1]])
rownames(fis.fst.mean) <- rownames(fis.fst.list[[1]])


#### GETTING MEANS OF Fis ACROSS LOCI and 50 replicates 

# Getting Fis in a list
fis.list <- list()
for(i in 1:length(dos.fst.fis)){
  fis.list[[i]] <- dos.fst.fis[[i]]$Fi}

# Getting the Fis list into array
for(i in 1:length(fis.list)){
  fis.list[[i]] <- sapply(fis.list[[i]], "[", 1:max(lengths(fis.list[[i]])))}
fis.array <- array(unlist(fis.list), dim=c(9,7,50))

# Taking the mean across repetitions
fis.mean <- apply(fis.array, c(2), median, na.rm=TRUE)



###########################################
#   Getting the ci of Fis in a list       
fis.ci.list <- list()
for(i in 1:length(dos.fis.ci)){
  fis.ci.list[[i]] <- dos.fis.ci[[i]]$fis.ci}

# Converting the ci.list to array
fis.ci.array <- array(unlist(fis.ci.list), dim=c(nrow(fis.ci.list[[1]]), 
                                             ncol(fis.ci.list[[1]]), 
                                             length(fis.ci.list)))

# Taking the mean across repetition
fis.ci.mean <- apply(fis.ci.array, c(1,2), mean, na.rm=TRUE)

# Putting names to columns and rows
colnames(fis.ci.mean) <- colnames(fis.ci.list[[1]])

# Getting a dataframe with mean Fis and ci
fis.ci <- cbind(fis.ci.mean, fis.fst.mean[1,-8])
#fis.ci <- cbind(fis.ci.mean, fis.mean)

# Putting names to columns
colnames(fis.ci)[3] <- "Fis"

# Converting array to dataframe
fis.ci <- as.data.frame(fis.ci)

# Adding the name of pop
fis.ci$pop <- colnames(fis.fst.mean)[1:7]

# Plotting the outcome
plot.fis.ci.4ind <- ggplot(fis.ci, aes(x=factor(pop, level=site.lat.order2), y= Fis)) +
                    geom_point() +
                    geom_errorbar(aes(ymin = ll, ymax = hl))


###########################################
# Getting the ci of Fst in a list
fst.ci.list <- list()
for(i in 1:length(dos.fst.ci)){
  fst.ci.list[[i]] <- dos.fst.ci[[i]]$ci
}

# Converting the ci.list to array
fst.ci.array <- array(unlist(fst.ci.list), dim=c(nrow(fst.ci.list[[1]]), ncol(fst.ci.list[[1]]), length(fst.ci.list)))

# Taking the mean across repetition
fst.ci.mean <- t(apply(fst.ci.array , c(1,2), mean, na.rm=TRUE))

# Getting the fst in a list
fst.list <- list()
for(i in 1:length(dos.fst.ci)){
  fst.list[[i]] <- as.data.frame(dos.fst.ci[[i]]$betaiovl)
}

# Converting the fst.list to array
fst.array <- array(unlist(fst.list), dim=c(nrow(fst.list[[1]]), ncol(fst.list[[1]]), length(fst.list)))

# Taking the mean across repetition
fst.mean <- apply(fst.array , c(1,2), mean, na.rm=TRUE)

# Getting a dataframe with mean Fis and ci
fst.ci <- cbind(fst.ci.mean, fst.mean)

# Putting names to columns
colnames(fst.ci) <- c("ll", "hl", "Fst")

# Converting aray to dataframe
fst.ci <- as.data.frame(fst.ci)

# Adding the name of pop
fst.ci$pop <- colnames(fs.mean)[1:7]

# Plotting the outcome
plot.fst.ci.4ind <- ggplot(fst.ci, aes(x=factor(pop, level=site.lat.order3), y= Fst)) +
                    geom_point() +
                    geom_errorbar(aes(ymin = ll, ymax = hl))


###########################################
# Getting a table with Hobs, Hexp and Fis
stat.all <- merge(summ.Hobs, summ.Hexp, by="pop")
stat.all <- merge(stat.all, fis.ci, by="pop")
stat.all <- stat.all[,c(1,2,3,4,5,8,6,7)]# change the order of columns
colnames(stat.all) <- c("Population", "Ho (mean)", "Ho (SE)", "He (mean)", 'He (SE)', "Fis", "Fis LowCI", "Fis HighCI")


write.table(stat.all, "stat.all.txt")

##################################################################################
############        HET accounting for kinship (BestHet)              ############


# Converting genind object to hierfstat
hier.all <- genind2hierfstat(wren.n5nb.genind, pop=pop(wren.n5nb.genind))

# Replacing NA with zeros
hier.all[is.na(hier.all)] <- 0

# Converting to dosage format
dos.all <- as.data.frame(fstat2dos(hier.all[,-1], diploid=TRUE))

# Dropping the columns with zeros (NAs)
dos.all <- select(dos.all, -contains(".0"))

# Write dosage data
write.table(dos.all, "dos_all.txt", row.names = FALSE, col.names = FALSE)

# Getting the kinship coefficient
kin.coef <- beta.dosage(dos.all)

# Rescaling the kinship matrix
kin.res <- rescale_popkin(kin.coef, min_kinship = min(kin.coef))

# Look at the histogram
hist(kin.res)

# Write kinship coefficients
write.table(kin.res, "kin_coef1.txt", row.names = FALSE, col.names = FALSE)

# Adding a group breeding-code for the museum samples
{
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F002"] <- "CF080" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F026"] <- "CF080" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F010"] <- "CF081" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F025"] <- "CF081" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F008"] <- "CF082" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F016"] <- "CF082" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F022"] <- "CF082" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F004"] <- "CF083" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F013"] <- "CF083" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F014"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F011"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F015"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F009"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F001"] <- "CF085" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F020"] <- "CF085" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F023"] <- "CF085" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F019"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F021"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F024"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F003"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F006"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F007"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F017"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F012"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F005"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "1836"] <- "CZ080" # LasGolondrinas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "1835"] <- "CZ080" # LasGolondrinas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "93"] <- "CF088" # Cazaderos
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "94"] <- "CF088" # Cazaderos
}

# Creating a matrix to match individuals in the same breeding group (1)
eq.grp <- matrix('NA', nrow=112, ncol=112)
grp <- pop.lab$Lev5_grp
for(i in 1:length(grp)){
  for(j in 1:length(grp)){
    ifelse(grp[i]==grp[j], eq.grp[i,j]<-1, eq.grp[i,j]<-0)
  }
}
diag(eq.grp) <- 0

# Creating a matrix with simplified kinship coefficient
kin <- matrix('NA', nrow=112, ncol=112)
kin <- ifelse(eq.grp == 1 & kin.res>0.6, 0.25, 0)
diag(kin) <- 0.5

# Checking the number of related individuals
sum(kin)/0.25-(112/0.5)

# Write kinship coefficients
write.table(kin, "kin_coef_simp.txt", row.names = FALSE, col.names = FALSE)



###########################################
##### Working Kinship Coefficient 
kin.coef.w <- read.table("kin_coef1.txt", header = TRUE, sep="\t")

# making it long format
kin.coef.l <- melt(kin.coef.w)

# Write kinship coefficients
write.table(kin.coef.l, "kin.coef.l.txt", row.names = FALSE, col.names = FALSE)



###########################################
##### Per population

# Splitting genlight object by sampling location
gi.all.pop <- list()
for(i in 1:length(unique(pop(wren.n5nb.genind)))){
  gi.all.pop[[i]] <- popsub(wren.n5nb.genind, sublist=unique(pop(wren.n5nb.genind))[i])}


# Converting to hierfstat
hier.all <- list()
for(i in 1:length(gi.all.pop)){
  hier.all[[i]] <- genind2hierfstat(gi.all.pop[[i]], pop=pop(gi.all.pop[[i]]))}

# Dropping last row of NA
for(i in 1:length(hier.all)){
  hier.all[[i]] <- hier.all[[i]][-nrow(hier.all[[i]]),]
  # Replacing NA with zeros
  hier.all[[i]][is.na(hier.all[[i]])] <- 0
}

# Converting to Dosage format
dos.all <- list()
for(i in 1:length(hier.all)){
  dos.all[[i]] <- as.data.frame(fstat2dos(hier.all[[i]][,-1], diploid=TRUE))}

# Dropping the columns with zeros (NAs)
for(i in 1:length(dos.all)){
  dos.all[[i]] <- select(dos.all[[i]], -contains(".0"))}

# Write dosage data
write.table(dos.all[[1]], "dos_golondrinas.txt", row.names = FALSE, col.names = FALSE)

# Getting the kinship coefficient
kin.coef <- list()
for(i in 1:length(dos.all)){
  kin.coef[[i]] <- beta.dosage(dos.all[[i]])}

# Rescaling the kinship matrix
kin.res <- list()
for(i in 1:length(kin.coef)){
  kin.res[[i]] <- rescale_popkin(kin.coef[[i]], min_kinship=min(kin.coef[[i]]))}

# Look at the histogram
hist(kin.res[[1]])

# Adding a group breeding-code for the museum samples
{
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F002"] <- "CF080" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F026"] <- "CF080" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F010"] <- "CF081" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F025"] <- "CF081" # Celedin
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F008"] <- "CF082" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F016"] <- "CF082" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F022"] <- "CF082" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F004"] <- "CF083" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F013"] <- "CF083" # Chachapoyas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F014"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F011"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F015"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F009"] <- "CF084" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F001"] <- "CF085" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F020"] <- "CF085" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F023"] <- "CF085" # Jaen
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F019"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F021"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F024"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F003"] <- "CF086" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F006"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F007"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F017"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F012"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "F005"] <- "CF087" # Sullana
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "1836"] <- "CZ080" # LasGolondrinas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "1835"] <- "CZ080" # LasGolondrinas
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "93"] <- "CF088" # Cazaderos
  pop.lab$Lev5_grp[pop.lab$Lev6_ind == "94"] <- "CF088" # Cazaderos
}

# Creating a matrix to match individuals in the same breeding group (1)
grp <- list()
eq.grp <- list()
for(k in 1:length(kin.res)){
  grp[[k]] <- rownames(hier.all[[k]])
  grp[[k]] <- pop.lab$Lev5_grp[match(grp[[k]], pop.lab$Lev6_ind)]
  eq.grp[[k]] <- matrix(NA, nrow=length(grp[[k]]), ncol=length(grp[[k]]))
  for(i in 1:length(grp[[k]])){
    for(j in 1:length(grp[[k]])){
      ifelse(grp[[k]][i]==grp[[k]][j], eq.grp[[k]][i,j]<-1, eq.grp[[k]][i,j]<-0)}}
  diag(eq.grp[[k]]) <- 0}

# Creating a matrix with simplified kinship coefficient
kin <- list()
for(i in 1:length(kin.res)){
  kin[[i]] <- matrix(NA, nrow=nrow(kin.res[[i]]), ncol=ncol(kin.res[[i]]))
  kin[[i]] <- ifelse(eq.grp[[i]]==1 & kin.res[[i]]>0.4, 0.25, 0)
  diag(kin[[i]]) <- 0.5}

# Checking the number of related individuals
sum(kin[[1]])/0.25-(nrow(kin[[1]])/0.5)

# Write kinship coefficients
write.table(kin[[1]], "kin_golondrinas.txt", row.names = FALSE, col.names = FALSE)


##################################################################################
#######   Hier AMOVA Morphological Subspp, Pop, and Breeding Groups   ############


# Taking out populations with less than 4 individuals
levels <- pop.lab[which(pop.lab$Lev4_pop !='Manglareschurute'),]

# Taking out populations without  breeding group information
levels <- levels[!is.na(levels$Lev5_grp),]

# converting genind to hier object
wren.n5nb.hier.15pop <- genind2hierfstat(wren.n5nb.genind)

# Subsetting the hier object by the samples with got in levels
wren.n5nb.hier.15pop <- as.data.frame(wren.n5nb.hier.15pop[rownames(wren.n5nb.hier.15pop) %in% levels$Lev6_ind,])

# Performing the hierarchical AMOVA
hiervar <- varcomp.glob(levels=levels[,c("Lev2_ssp","Lev4_pop", "Lev5_grp")],
                        loci=wren.n5nb.hier.15pop[,-1])

### Getting a p-value of differentiation among groups

# Testing breeding groups within populations 
testbg.bwpop <- test.between.within(data=wren.n5nb.hier.15pop[,-1], within=levels$Lev4_pop, 
                               test.lev=levels$Lev5_grp, rand.unit=levels$Lev6_ind, nperm=1000)

# Testing populations within morphological subspecies 
testpop.bwssp <- test.between.within(data=wren.n5nb.hier.15pop[,-1], within=levels$Lev2_ssp, 
                                    test.lev=levels$Lev4_pop, rand.unit=levels$Lev6_ind, nperm=1000)



### AMOVA dropping the transitions populations

# Dropping transition populations 
wren.n5nb.genind.11pop <- popsub(wren.n5nb.genind, exclude=c("PatriciaPilar", "Chone", "Calceta", "Montecristi"))

# Check the population's names
popNames(wren.n5nb.genind.11pop)

# converting genind to hier object
wren.n5nb.hier.11pop <- genind2hierfstat(wren.n5nb.genind.11pop)

# Taking out populations with less than 4 individuals
levels.11pop <- pop.lab[which(pop.lab$Lev4_pop !="PatriciaPilar" &  pop.lab$Lev4_pop !="Chone" & pop.lab$Lev4_pop !="Calceta" &
                                pop.lab$Lev4_pop !="Montecristi" & pop.lab$Lev4_pop !="Manglareschurute"),]

# Taking out populations without  breeding group information
levels.11pop <- levels.11pop[!is.na(levels.11pop$Lev5_grp),]

# Subsetting the hier object by the samples with got in levels
wren.n5nb.hier.11pop <- as.data.frame(wren.n5nb.hier.11pop[rownames(wren.n5nb.hier.11pop) %in% levels.11pop$Lev6_ind,])


# Performing the hierarchical AMOVA
hiervar.11pop <- varcomp.glob(levels=levels.11pop[,c("Lev2_ssp","Lev4_pop", "Lev5_grp")],
                        loci=wren.n5nb.hier.11pop[,-1])

### Getting a p-value of differentiation among groups

# Testing breeding groups within populations 
testbg.bwpop.11pop <- test.between.within(data=wren.n5nb.hier.11pop[,-1], within=levels.11pop$Lev4_pop, 
                                    test.lev=levels.11pop$Lev5_grp, rand.unit=levels.11pop$Lev6_ind, nperm=1000)

# Testing populations within morphological subspecies 
testpop.bwssp.11pop <- test.between.within(data=wren.n5nb.hier.11pop[,-1], within=levels.11pop$Lev2_ssp, 
                                     test.lev=levels.11pop$Lev4_pop, rand.unit=levels.11pop$Lev6_ind, nperm=1000)

# Estimate the pairwise Fst
wren.n5nb.hier <- genind2hierfstat(wren.n5nb.genind)
wn5nb.pwFst <- pairwise.neifst(wren.n5nb.hier, diploid=TRUE)


par(mar=c(11,4,4,4)) # Setting the margins
temp <- as.matrix(wn5nb.pwFst)
diag(temp) <- NA
boxplot(temp, col=funky(nPop(wren.n5nb.genind)), las=3, ylab="Pairwise Fst")
mtext(text="Populations",
      side=1, line=10)




##################################################################################
############         Ho, Hs & Fis using relatedness index             ############

##################################################################################
# Data Management

setPop(wren.n5nb.genind) <- ~top1name

# Converting to genlight object
gl.wren <- gi2gl(wren.n5nb.genind)

# Creates a loc.metrics object to comply to a dartR object
df.loc <- data.frame(RepAvg=runif(nLoc(gl.wren)), CallRate=1)
gl.wren@other$loc.metrics <- df.loc
gl.wren.fix <- gl.compliance.check(gl.wren)


################################################################################
###################       Heterozygosity without filter              ###########


# Estimating Heretozygosity
gl.het <- gl.report.heterozygosity(gl.wren.fix)

# Estimating Diversity
gl.div <- gl.report.diversity(gl.wren.fix)


################################################################################
############ Getting no related individuals with kinship coefficient ###########

# Converting genind to hier object
hier.hwe.all2 <- genind2hierfstat(wren.n5nb.genind, pop=pop(wren.n5nb.genind))

# Replacing NA with zeros
hier.hwe.all2[is.na(hier.hwe.all2)] <- 0

# Convert to dosage format
hier.dos2 <- as.data.frame(fstat2dos(hier.hwe.all2[,-1]))

# Dropping the columns with zeros
hier.dos2 <- select(hier.dos2, -contains(".0"))

# Set row names
rownames(hier.dos2) <- rownames(hier.hwe.all2)

# Get groups with more than one individuals
grp.ind <- as.data.frame(table(pop.lab$Lev5_grp))
grp.ind <- grp.ind[which(grp.ind$Freq != 1),]

# Getting the pop data with groups with more than one individual
pop.grp.2 <- struc.pop[which(struc.pop$Lev5_grp %in% grp.ind$Var1), ]

# split pop.lab by breeding groups
kin.in.list <- split(pop.grp.2, f=pop.grp.2$Lev5_grp)

# Splitting the dosage data by groups
grp.list <- list()
  for(i in 1:length(kin.in.list)){
      grp.list[[i]] <- hier.dos2[which(rownames(hier.dos2) %in% kin.in.list[[i]]$Lev6_ind),]}

# Estimating the kinship coefficient
grp.kin.1 <- list()
for(i in 1:length(kin.in.list)){
  grp.kin.1[[i]] <- beta.dosage(grp.list[[i]])
}



# 
for(i in 1:length(grp.kin.1)){
  grp.kin.1[[i]][upper.tri(grp.kin.1[[i]])] <- 0
  grp.kin.1[[i]] <- as.matrix(grp.kin.1[[i]])
  diag(grp.kin.1[[i]]) <- 1
  grp.kin.1[[i]] <- as.data.frame(as.table(grp.kin.1[[i]]))
  grp.kin.1[[i]] <- grp.kin.1[[i]][grp.kin.1[[i]]$Freq != 0, ]
}


# Converting in dataframe
grp.kin <- do.call(rbind.data.frame, grp.kin.1)

# Adding ID of the breeding groups
grp.kin <- merge(grp.kin, pop.lab[,c(1,6)], by.x='Var2', by.y='Lev6_ind')
grp.kin <- merge(grp.kin, pop.lab[,c(1,6)], by.x='Var1', by.y='Lev6_ind')
grp.kin.2 <- grp.kin[which(grp.kin$Freq != 1),]
dim(grp.kin)
write.table(grp.kin.2, file="kinship_coefficient.txt", sep="\t", row.names=FALSE)

# Samples that have positive kinship coef
id.nok <- grp.kin.2[which(grp.kin.2$Freq >= 0),]
 
# Selection of independent samples 
# I could use either of these vectors, but no samples from both
# vector 1
vec1 <- c('41B', '42', '47B', '45', '48', '55', '54', '64', '65B', '73', '75B', 
          'F008', 'F011', 'F020', 'F023', 'F021', 'F003',  'F006', '40', 'F017', 'F005')

# vector 2
vec2 <- c('43', '46', '49', '56', '66', '72', '77B', 'F016', 'F014', 'F001', 
          'F019', 'F007', '39', 'F012')

# Checking if the vectors have elements in common
# Intersection among these should be zero
#intersect(vec1, vec2)

# Independent samples
id.ok <- setdiff(unique(pop.lab$Lev6_ind), vec1)

length(id.ok)

# Select hierfstat object with independent samples
hier.wren.ind <- wren.n5nb.genind[indNames(wren.n5nb.genind) %in% id.ok]

##################################################################################
### FSt 

# Converting genind object to hierfstat
hier.wren <- genind2hierfstat(hier.wren.ind)

# Estimating general genetics statatistics
library(hierfstat)
basic.hier <- hierfstat::basic.stats(hier.wren)

# Estimating Fst
fst.wc <- pairwise.WCfst(hier.wren)
fst.nei <- pairwise.neifst(hier.wren)
fst.betas <- pairwise.betas(hier.wren)

# Write table
write.table(fst.nei, "fst_nei.csv", row.names = TRUE, sep=",", quote=FALSE, col.names =TRUE)

##################################################################################
## Heterozygosity

## Filtering gl.wren.fix with independent samples based on Kinship coefficients
gl.wren.ind <- gl.keep.ind(gl.wren.fix, ind.list=id.ok, recalc = TRUE)

# Estimating Heretozygosity per individuals
gl.het.ind <- gl.report.heterozygosity(gl.wren.ind, method="ind")
gl.het.ind$He <- 2*gl.het.ind$f.hom.ref*gl.het.ind$f.hom.alt
gl.het.ind2 <- merge(struc.pop[,c(1,32)], gl.het.ind, by.x="Lev6_ind", by.y="ind.name")
colnames(gl.het.ind2)[6] <- "Ho"
gl.het.ind2 <- rbind(gl.het.ind2[,c(1,2,3)], gl.het.ind2[,c(1,2,6)])
colnames(gl.het.ind2) <- c("Sample", "GC", "Value")
gl.het.ind2$Statistic <- c(rep("Observed", 91), rep("Expected", 91))

library(dplyr)

gl.het.sum <- gl.het.ind2 %>%
	dplyr::group_by(GC, Statistic) %>%
	dplyr::summarize(
		y0=quantile(Value, 0.025),
		y25=quantile(Value, 0.25),
		y50 = median(Value), 
    y75 = quantile(Value, 0.75), 
   	y100 = quantile(Value, 0.975))

write.csv(gl.het.sum, file="het_gc.csv", row.names=FALSE)

jpeg(file="het_gc.jpeg", width=10000, height=5000, units="px", res=300, quality=100)
ggplot(gl.het.ind2, aes(x=factor(GC, level=spp.ord), fill=Statistic), 
 	position=position_dodge(width=0.8))+
  geom_boxplot(aes(y=Value), position=position_dodge(width=0.8), outlier.shape=NA, alpha=0.2) +
  geom_point(aes(y=Value, color=Statistic), position=position_jitterdodge(dodge.width=0.8), size=3) +
  scale_fill_manual(values=c("blue2", "red2")) +
  scale_color_manual(values=c("blue2", "red2")) + 
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=1,face="bold")) +
  theme(axis.title.x = eloement_text(vjust=-1)) +
  theme(axis.title.y = element_text(vjust=2)) +
  theme(text = element_text(size = 40)) +
  xlab("Subspecies") + ylab("Heterozygosity")
dev.off()


# Estimating Heretozygosity per individuals
gl.het.pop <- gl.report.heterozygosity(gl.wren.ind, method="pop")










################################################################################
###################                     AMOVA                        ###########

# Setting the strata hierarchies

# Nesting the C. fasicatus and pallescens Subspecies in species
struc.pop$subsp <- struc.pop$top1name
struc.pop$subsp <- gsub("C. f. pallescens North", "C. f. pallescens", struc.pop$subsp) 
struc.pop$subsp <- gsub("C. f. pallescens South", "C. f. pallescens", struc.pop$subsp)

# Nesting the C. f . pallescens genetic cluster in species
struc.pop$sp <- struc.pop$subsp
struc.pop$sp <- gsub("C. f. pallescens", "C. fasciatus", struc.pop$sp) 
struc.pop$sp <- gsub("C. f. fasciatus", "C. fasciatus", struc.pop$sp) 


# partiotion of variance with hierfstat
varcomp.hier <- varcomp.glob(levels=struc.pop[struc.pop$Lev6_ind %in% id.ok, c(32,1)],
loci=as.matrix(hier.wren[,-1]))


data<- read.table("examplehier.txt", header=TRUE)

strata(gl.wren.ind) <- struc.pop[struc.pop$Lev6_ind %in% id.ok, c(1,32,38,39)]
#strata(gl.wren.fix) <- struc.pop[,c(1,32,38,39)]
setPop(gl.wren.ind) <- ~top1name

# Performing an AMOVA
gl.amv <- gl.amova(gl.wren.ind, permutations = 10000)

# Getting the FSt
gl.fst <- gl.fst.pop(gl.wren.fix, nboots=10000, nclusters=10)

# Getting the Fst values in a table
fst.gc <- gl.fst$Bootstraps[,c(1:2,10003:10006)]

# Write table
write.table(fst.gc, "fst_gc.csv", row.names = FALSE, sep=",", quote=FALSE, col.names =TRUE)


# Getting the Fst values in a table
amv.gc <- gl.amv$Bootstraps[,c(1:2,10003:10006)]


################################################################################
###################       Plotting the Heterozygosity                ###########

# Getting margin of error of heterozygosity
Ho.me <- 1.96*(gl.het.po$HoSD/sqrt(gl.het.po$nInd))
He.me <- 1.96*(gl.het.po$HeSD/sqrt(gl.het.po$nInd))

# Adding the margin of error
gl.het.po$Ho.me <- Ho.me
gl.het.po$He.me <- He.me

# Adding the lower and upper boundaries
gl.het.po$Ho.lb <- gl.het.po$Ho - gl.het.po$Ho.me
gl.het.po$Ho.ub <- gl.het.po$Ho + gl.het.po$Ho.me
gl.het.po$He.lb <- gl.het.po$He - gl.het.po$He.me
gl.het.po$He.ub <- gl.het.po$He + gl.het.po$He.me

# Converting to long format
gl.het <- gl.het.po[,c(1,8,19,12,20)] # Removing FIS
colnames(gl.het)[4:5]  <- colnames(gl.het)[2:3]
het.long <- rbind(gl.het[,c(1,2,3)], gl.het[,c(1,4,5)])
het.long$het <- c(rep('Observed', 4), rep('Expected', 4))
colnames(het.long) <- c("GC", "mean", 'me', 'Statistic')

# Write Heterozygosity
write.table(gl.het.po, "he_gc.csv", row.names = FALSE, sep=",", quote=FALSE, col.names =TRUE)

jpeg(file="het_site.jpeg", width=10000, height=5000, units="px", res=300, quality=100)
ggplot(het.long, aes(x=factor(GC, level=spp.ord), y=mean, color=Statistic)) +
  geom_point(position=position_dodge(width=0.4), size=8) +
  geom_errorbar(aes(ymin = mean-me , ymax = mean+me), width=0.4, size=3,
                position=position_dodge(width=0.4)) +
  scale_color_manual(values=c("blue2", "red2")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=1,face="bold")) +
  theme(axis.title.x = element_text(vjust=-1)) +
  theme(axis.title.y = element_text(vjust=2)) +
  theme(text = element_text(size = 40)) +
  xlab("Subspecies") + ylab("Heterozygosity")
dev.off()


##################################################################################




##################################################################################

# Select randomly 12 rows per category in the column top1name from the dataframe struc.pop
# and save the result in a new dataframe called struc.pop2

# Creating a number of resampling
n_resample <- 10
rs.list <- list()
gl.12.temp <- list()
gl.fst.12.temp <- list()
fst.12 <- list()
pvalue.12 <- list()
nboot.12 <- list()
het.12 <- list()
nboot <- 10000
a=nboot+3
b=nboot+6
nsamp=12

library(dplyr)
for(i in 1:n_resample){
 rs.list[[i]] <- struc.pop %>% 
 group_by(top1name) %>% 
 sample_n(size=nsamp, replace = TRUE)

 # Filtering gl.wren.ind using the new data set struc.pop
 gl.12.temp[[i]] <- gl.keep.ind(gl.wren.ind, ind.list=rs.list[[i]]$Lev6_ind, recalc = TRUE)

 # Get fst among populations per each resampling
 gl.fst.12.temp[[i]] <- gl.fst.pop(gl.12.temp[[i]], nboots=nboot, nclusters=10)

 fst.12[[i]] <- as.data.frame(gl.fst.12.temp[[i]]$Fsts)
 pvalue.12[[i]] <- as.data.frame(gl.fst.12.temp[[i]]$Pvalues)
 nboot.12[[i]] <- as.data.frame(gl.fst.12.temp[[i]]$Bootstraps)

 #Getting heterozygosity per each resampling
 het.12[[i]] <- gl.report.heterozygosity(gl.12.temp[[i]], method="pop")
}


# Convert fst.12 to an array
L <- array(unlist(fst.12), dim = c(4,4,n_resample))
P <- array(unlist(pvalue.12), dim = c(4,4,n_resample))
N <- list()
for(i in 1:n_resample){
N[[i]] <- as.matrix(nboot.12[[i]][,a:b])
}
N <- apply(simplify2array(N), c(1,2), mean, na.rm=TRUE)
# Remove first column of het.12
H <- lapply(het.12, function(x) as.matrix(x[,-1]))

gl.fst.12.temp[[1]]$Fsts
L[[1,,3]]
fst.12[[1]]
dim(het.12[[1]])
length(het.12)
class(H[[1]])

# Getting average across the resamplings
gl.fst.12 <- apply(L, c(1,2), mean, na.rm=TRUE)
gl.pvalue.12 <- apply(P, c(1,2), mean, na.rm=TRUE)
gl.nboot.12 <- apply(simplify2array(N), c(1,2), mean, na.rm=TRUE)
gl.het.12 <- apply(simplify2array(H), c(1,2), mean, na.rm=TRUE)

# Putting colnames and rownames
colnames(gl.fst.12) <- colnames(fst.12[[1]])
rownames(gl.fst.12) <- rownames(fst.12[[1]])
colnames(gl.pvalue.12) <- colnames(fst.12[[1]])
rownames(gl.pvalue.12) <- rownames(fst.12[[1]])
colnames(gl.nboot.12) <- colnames(nboot.12[[1]][,a:b])

gl.fst.12
gl.pvalue.12
gl.nboot.12
gl.het.12


##################################################################################
# Samples selection for Whole Genome Sequencing Project
# This uses the kin coefficients estimated above but
# all the code inside this section is not part of this project.

# Reading data with original collection codes
sample_lib <- read.csv("SamplesReseq.csv", header = TRUE)
sample_lib$No <- sample_lib$?..No
sample_lib <- sample_lib[,-1]
View(sample_lib)

# Get independent samples
sample_sel <- sample_lib[sample_lib$No %in% id.ok, ]
sample_sel$Museum.code <- substr(sample_sel$Museum.code,3,7)
View(sample_sel)

# Reading and selecting independent samples from FMNH
fmnh <- read.csv("Genus_Campy.csv", header = TRUE)
fmnh$Museum.code <- fmnh$?..Museum.code
fmnh <- fmnh[,-1]
fmnh.sel <- fmnh[fmnh$Museum.code %in% sample_sel$Museum.code, ]
View(fmnh.sel)

# Write the samples in txt
write.table(fmnh.sel, "fmnh.sel.txt", sep="\t", row.names = FALSE, col.names = TRUE)



#################################################################################
# Getting results

spp.ord <- c("C. z. brevirostris", "C. f. pallescens North", "C. f. pallescens South", "C. f. fasciatus")

# Running function
hier.hwe.basic <- basic.stats(hier.hwe.ind)

# Running the plot function
plot.basic.hier.stats(x=hier.hwe.basic, ord=site.lat.order2, frame=c(4,4))


################################################################################
# Function to get a table
mean.ci <- function(x, vec, num){
  
  m.l <- list(ho <- hier.hwe.basic$Ho,
              hs <- hier.hwe.basic$Hs,
              Fis <- hier.hwe.basic$Fis)
  
 med.ci <- matrix(NA, ncol(m.l[[1]]), 9)
 
fun.ci <- function(idx){ 
  ci.lw <- rep(NA, ncol(m.l[[1]]))
  ci.up <- rep(NA, ncol(m.l[[1]]))
  m <- rep(NA, ncol(m.l[[1]]))
  sd <- rep(NA, ncol(m.l[[1]]))
  se <- rep(NA, ncol(m.l[[1]]))
  for(j in 1:ncol(m.l[[idx]])){
  m[j] <- mean(m.l[[idx]][,j], na.rm=TRUE)
  sd[j] <- sd(m.l[[idx]][,j], na.rm=TRUE)
  se <- sd/sqrt(length(m.l[[idx]][,j])) 
  ci.lw[j] <- m[j]-1.96*se[j]
  ci.up[j] <- m[j]+1.96*se[j]}
  ci <- cbind(ci.lw, ci.up)
  return(ci)
}

 temp1 <- fun.ci(1)
 temp2 <- fun.ci(2) 
 temp3 <- fun.ci(3)
 
for(j in 1:ncol(m.l[[1]])){
    
  med.ci[j,1] <- mean(m.l[[1]][,j], na.rm=TRUE)
  med.ci[j,4] <- mean(m.l[[2]][,j], na.rm=TRUE)
  med.ci[j,7] <- mean(m.l[[3]][,j], na.rm=TRUE)
  
  med.ci[j,2] <- temp1[j,1]
  med.ci[j,3] <- temp1[j,2]
  
  med.ci[j,5] <- temp2[j,1]
  med.ci[j,6] <- temp2[j,2]
  
  med.ci[j,8] <- temp3[j,1]
  med.ci[j,9] <- temp3[j,2]
 }

   rownames(med.ci) <- colnames(m.l[[1]])
 
   med.ci <- cbind(rownames(med.ci), med.ci)
   
   med.ci <- med.ci[match(vec, med.ci[,1]), ] 
   
   colnames(med.ci) <- c("n", "Mean_Ho", "Lw_Ci_Ho", "Up_Ci_Ho",
                       "Mean_Hs", "Lw_Ci_Hs", "Up_Ci_Hs",
                       "Mean_Fis", "Lw_Ci_Fis", "Up_Ci_Fis")
   
   med.ci <- as.data.frame(med.ci)
   
   for(i in 1:ncol(med.ci)){
     med.ci[,i]<-as.numeric(med.ci[,i])
   }
   
   med.ci <- round(med.ci, digits=3)
   
   num <- num[match(vec, num[,1]), ] 
   
   med.ci[,1] <- num[,2]
   
   med.ci <- cbind(rownames(med.ci), med.ci)
   
   colnames(med.ci)[1] <- 'Pop'
   
   return(med.ci)

}
################################################################################


# Running the function
basic.table <- mean.ci(hier.hwe.basic, vec=site.lat.order2, num=pop.num)

write.table(basic.table, "Fis_H.txt", quote=FALSE, sep="\t", col.names=TRUE)

# Check the table
View(basic.table)

# Formula for checking. I can estimate individual values to check those in the DF
mean((hier.hwe.basic$Hs[,4])) + qt(0.975, df=4409-1)*sd(hier.hwe.basic$Hs[,4])/sqrt(4409)


# Plotting the results
ggplot(basic.table, aes(x=factor(Pop, level=site.lat.order2), y=Mean_Ho)) +       
  geom_point() +
  geom_errorbar(aes(ymin = Lw_Ci_Ho , ymax = Up_Ci_Ho)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("Sampling Sites") + ylab("Observed Heterozygosity")


ggplot(basic.table, aes(x=factor(Pop, level=site.lat.order2), y=Mean_Hs)) +       
  geom_point() +
  geom_errorbar(aes(ymin = Lw_Ci_Hs , ymax = Up_Ci_Hs)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  xlab("Sampling Sites") + ylab("Expected Heterozygosity")


ggplot(basic.table, aes(x=factor(Pop, level=spp.ord), y=Mean_Fis)) +       
  geom_point() +
  geom_errorbar(aes(ymin = Lw_Ci_Fis , ymax = Up_Ci_Fis)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("Sampling Sites") + ylab("Inbreeding Coefficient")

### Both heterozygosity in the same plot

# Converting to long format
basic.het <- basic.table[,-c(2,9:11)] # Removing FIS
colnames(basic.het)[5:7]  <- colnames(basic.het)[2:4]
het.long <- rbind(basic.het[,c(1,2,3,4)], basic.het[,c(1,5,6,7)])
het.long$het <- c(rep('Observed', 16), rep('Expected', 16))
colnames(het.long)[2:5] <- c("mean", 'lw.ci', 'up.ci', 'Statistic')


# Plotting the results

library(ggplot2) 
library(grid)
library(RColorBrewer)


make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
    matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"), 
    interpolate = TRUE
  )
}

g <- make_gradient(deg = 180, n = 500, cols = brewer.pal(9, "RdBu"))


jpeg(file="het_site.jpeg", width=8000, height=5000, units="px", res=300, quality=100)


ggplot(het.long, aes(x=factor(Pop, level=site.lat.order2), y=mean, color=Statistic)) +
    geom_point() +
  geom_errorbar(aes(ymin = lw.ci , ymax = up.ci), width=0.2, size=1) +
  scale_color_manual(values=c("gray49", "black")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=1)) +
  theme(text = element_text(size = 15)) +
  xlab("Subspecies") + ylab("Heterozygosity")



ggplot(het.long2, aes(x=factor(Pop, level=spp.ord), y=mean, color=Statistic)) +
  geom_point() +
  geom_errorbar(aes(ymin = lw.ci , ymax = up.ci), width=0.1, size=1) +
  scale_color_manual(values=c("gray49", "black")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust=1)) +
  theme(text = element_text(size = 15)) +
  xlab("Subspecies") + ylab("Heterozygosity")

    
dev.off()




annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
geom_segment(aes(x = 3, y = 0.11, xend = 1, yend = 0.11, size=5),
            arrow = arrow(length = unit(1, "cm")), colour="black", show.legend = F) +
annotate(geom="text", size=15, x=2, y=0.1, label="North", color="black") +
geom_rect(aes(xmin = 13.5, xmax = 16.5, ymin = 0.06, ymax = 0.15), 
          alpha=0.01, colour="black", fill="#EE2C2C", size=1, show.legend = F)





install.packages("Rmisc")
library(Rmisc)

mean.ci.het <- as.data.frame(t(do.call(cbind,lapply(hier.hwe.basic[3:5], function(x)apply(na.omit(x),2,CI)))))
pop <- as.vector(rep(row.names(mean.ci.het)[1:16], 3))
mean.ci.het$Pop <- rep(row.names(mean.ci.het)[1:16], 3)
mean.ci.het$Stat <- rep(c("Ho","Hs","Fis"), each=16)
View(mean.ci.het)


jpeg(file="het_site.jpeg", width=10000, height=5000, units="px", res=300, quality=100)

ggplot(mean.ci.het[mean.ci.het$Stat != "Fis",], aes(x=factor(Pop, level=site.lat.order2), y=mean, color=Stat)) +       
  geom_point() +
  geom_errorbar(aes(ymin = lower , ymax = upper)) +
  scale_color_manual(values=c("red", "blue")) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  theme(text = element_text(size = 30)) +
  xlab("Sampling Sites") + ylab("Heterozygosity")

dev.off()



##################################################################################

save.image("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/R/Genepop/Fst_diversity.RData")
