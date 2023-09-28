################################################################################
################            ANALYSIS WITH INTROGRESS            ################


rm(list=ls())

# Set the working directory
# Setting the working directory
setwd("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Introgress")

# Load the outcomes so far
load("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Introgress/Introgress.RData")

memory.limit(size=100000)


################################################################################
############                Installing Packages                       ##########

# Installing packages
install.packages("introgress")
install.packages("ggplot2")
install.packages("rgdal")
install.packages("devtools")
devtools::install_github("btmartin721/ClineHelpR")


# Calling packages
library(introgress)
library(ggplot2)
library(rgdal)
library("ClineHelpR")

################################################################################
################                 Data Formatting                ################

# Reading the STRUCTURE with info with 16 pop from the struc2 model
struct2 <- read.table("populations_n5nb_str2.txt", header=FALSE, sep="\t")

# Preparing admix.gen data object
# transpose the structure file
struct2 <- struct2[, c(2,1,c(3:4411))]

# Transposing the database
struc2.t <- t(struct2)
struc2.t <- struc2.t[,-1]
struc2.t[struc2.t == 0] <- NA

#### We need to classified individuals of parental populations
# Read the clumpp outfile for K=2
k2 <- read.table("K2.outfile")
k2 <- k2[,c(6,7)]

# Read the pop file used in the denovo assembly in Stacks
# This is the same order for samples
stack.ord <- read.table("popmap_wren_n5nb_sh11.txt", col.names = c("Ind", "ord"))
stack.ord$ord <- seq(1,112, 1)

# Putting names of individuals in k2 matrix
k2$Ind <- stack.ord[,1]
colnames(k2) <- c("CF", "CZ", "Ind")

# Reading the data where I have the right Sampling Location Names
poplab <- read.table("C:/Users/USER/Dropbox/Thesis/Molecular_Wrens/Wren_1/R/poplab.txt", header=FALSE, sep="\t")
colnames(poplab) <- c("Ind", "pop")



#### GETTING PARENTALS

#### P1 C zonatus
# Getting most northern sampling locations for parental P1
pop.p1 <- poplab[poplab$pop %in% c("Las Golondrinas", "PVM", "Pedro Carbo"),]
  
# Getting a vector with individuals with Q>0.9 of parental for C zonatus (P1)
P1 <- k2[k2$CZ >= 0.9,3]
P1 <- P1[P1 %in% pop.p1$Ind]

# Subsetting struct2 to get parental P1 matrix
P1.m <- struct2[struct2$V1 %in% P1,]
P1.m <- P1.m[,-c(1:2)]
P1.m <- rbind(struct2[1,-c(1:2)], P1.m)
P1.m <- t(P1.m)
P1.m <- P1.m[,-1]


#### P2 C fasciatus
# Getting most southern sampling locations for parental P2
pop.p2 <- poplab[poplab$pop %in% c("Jaen", "Chachapoyas", "Celedin"),]

# Getting a vector with individuals with Q>0.9 of parental for C fasciatus (P2)
P2 <- k2[k2$CF >= 0.9,3]
P2 <- P2[P2 %in% pop.p2$Ind]

# Subsetting struct2 to get parental P2 matrix
P2.m <- struct2[struct2$V1 %in% P2,]
P2.m <- P2.m[,-c(1:2)]
P2.m <- rbind(struct2[1,-c(1:2)], P2.m)
P2.m <- t(P2.m)
P2.m <- P2.m[,-1]

# Getting the loci.data matrix
loci.d <- as.data.frame(t(struct2[1,-c(1:2)]))
loci.d$C <- rep("C", 4409)
loci.d <- as.matrix(loci.d)
colnames(loci.d) <- c("locus", "type")


################################################################################
################               Running Introgress               ################

# Preparing data
intro.data <- prepare.data(admix.gen=struc2.t,
                         loci.data=loci.d,
                         parental1=P1.m,
                         parental2=P2.m,
                         pop.id=TRUE,
                         ind.id=TRUE,
                         fixed=FALSE,
                         sep.columns=TRUE)

# Getting the hybrid index
hi.index <- est.h(introgress.data=intro.data, loci.data=loci.d)

# Dropping off the sampling locations for C brunneicapillus
poplab2 <- subset(poplab, !(poplab$pop%in%c("Texas", "Arizona", "Nevada")))

# Adding the sampling location in the hybrid index
hi.index.pop <- cbind(hi.index, poplab2$pop)

# Columns names for Hybrid Index
colnames(hi.index.pop)[4] <- c('pop')

# Write the hyrbrid index in a table
write.csv(hi.index.pop, "hi.index.csv", sep=",", row.names = FALSE, col.names = TRUE)

#### Plotting the Hybrid Index

# Reading sampling locations with latitude order
pop.lat <- read.table("wren_pop.txt", header=TRUE)

# Dropping off the sampling locations for C brunneicapillus
pop.lat <- subset(pop.lat, !(pop.lat$Lev4_pop%in%c("Texas", "Arizona", "Nevada")))

# Ordering the sampling locations by latitude
pop.lat <- pop.lat[order(pop.lat$Ord_Lat, decreasing=FALSE),]

# Getting the sampling locations in the latitudinal order
site.lat.order <- unique(pop.lat$Lev4_pop)

# Plotting the Hybrid index per sampling location 
plot.hi.index <- ggplot(hi.index.pop, aes(x=factor(pop, level=site.lat.order), y=h)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title="Hybrid Index by Sampling Locations", y="Hybrid Index", x="Sampling Locations")
 



#### Plotting Hybrid Index by Individual

# Data with everything including labels and coordinates per individuals
labels.all <- read.table("Labels.field.csv", header=TRUE, sep=",", na.strings=TRUE, fileEncoding = "UTF-8-BOM")
coor <- merge(stack.ord, labels.all, by="Ind")[,c(1,2,9,10)]
coor <- coor[order(coor$ord,  decreasing = FALSE),]

# Making the coordinates numeric
coor2 <- cbind(as.numeric(as.character(coor$Long)), as.numeric(as.character(coor$Lat)))

# Setting the coordinates system as Long/Lat
cord.dec = SpatialPoints(cbind(coor2[,1], coor2[,2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))

# Converting to UTM, I'm using UTM zone 17s (epsg=32717)
# Most of the territory is in this zone
coor.utm <- spTransform(cord.dec, CRS("+init=epsg:32717"))
colnames(coor.utm@coords) <- c("X", "Y")

# Adding the coordinates to the Hybrid index
hi.index.pop <- cbind(hi.index.pop, coor.utm@coords, stack.ord)

# Ordering individuals by Latitude
hi.index.pop <- hi.index.pop[order(hi.index.pop$Y, decreasing = TRUE),]


## Reading table with Structure ancestry proportions
struc <- read.table("k4popfile.txt", header = TRUE, sep="\t")

site <- c("Pedro Carbo", "Las Golondrinas", "PVM", "Patricia Pilar", "Chone", "Calceta", "Montecristi", "Machalilla", 
          "Manglares Churute", "Arenillas", "Cazaderos", "Zapotillo", "Sullana", "Jaen", "Chachapoyas", "Celedin")

spp.ord <- c("C. z. brevirostris", "C. f. pallescens North", "C. f. pallescens South", "C. f. fasciatus")

# Merging everything together
hi.index.pop2 <- merge(struc, hi.index.pop, by.x="sample", by.y="Ind")


# creating a personalized Blue&Red gradient pallete for 16 colors
colr.3 <- c("Pedro Carbo"="#000080", "Las Golondrinas"="#0041C2", PVM="#1569C7",  "Patricia Pilar"="#357EC7", 
            Chone="#3090C7", Calceta="#79BAEC", Montecristi="#7BCCB5", Machalilla="#4CC417",
            "Manglares Churute"="#B2C248", Arenillas="#EDE275", Cazaderos="#FFD801", Sullana="#FBB117", 
            Zapotillo="#FF6C00",Jaen="#FF6100", Chachapoyas="#FF3600", Celedin="#FF0000")


# Choosing a new pallete with four colors
colr.2 <- c("C. f. fasciatus"="#FF0000","C. z. brevirostris"="#157DEC","C. f. pallescens North"="#4CC417", "C. f. pallescens South"="#FFD801")



#### Plotting the Hybrid index per sampling location 
plot.hi.index.ind2 <- ggplot() +
  geom_point(hi.index.pop2, mapping=aes(x=reorder(sample,-Y), y=h, color=factor(top1name, levels = spp.ord)), size=5) +
  geom_errorbar(hi.index.pop2, mapping=aes(x=reorder(sample,-Y), ymin=lower, ymax=upper, color=factor(top1name, levels = spp.ord)), size=2, width=2)+
  scale_color_manual("Genetic Cluster K=4", values = colr.2) +
  
  #geom_text(aes(x=103, label="C. f. fasciatus", y=0.90), colour="#F62217",  angle=0, size=12, face='bold')+
  #geom_text(aes(x=11, label="C. z. brevirostris", y=0.05), colour="#000080",  angle=0, size=12, face='bold')+
  theme_bw()+
  labs(y="Hybrid Index", x="Individuals")+
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
        
        legend.title=element_text(size=35), 
        legend.text=element_text(size=33),
        axis.title.x = element_text(size=40),
        axis.title.y = element_text(angle=90, size=40),
        axis.text.x = element_text(angle = 90, size=18, vjust = 0.5, hjust=1),
        axis.text.y = element_text(angle = 0, size=18, vjust = 0.5, hjust=1),
        legend.position = c(.2,.8))

# Saving the plot with right proportions
jpeg(file="HI_ind.jpg", width=8000, height=5000, res=300, quality=100)
plot.hi.index.ind2 
dev.off()



ggsave("HI_ind.jpg", plot.hi.index.ind2, scale=1, dpi=300, limitsize=FALSE)




#### Genomic Clines

# Preparing loci.data 
locus.data <- as.data.frame(rownames(intro.data$Count.matrix))
locus.data$type <- "C"
colnames(locus.data)<- c('locus', "type")
head(locus.data)


# Estimating Genomic Clines
clines.out <- genomic.clines(introgress.data = intro.data,
                             hi.index = hi.index.pop,
                             loci.data = locus.data,
                             sig.test = TRUE,
                             method = "permutation",
                             classification = TRUE,
                             n.reps=1000)


composite.clines(clines.out)



save.image("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Introgress/Introgress.RData")
rm(list=ls())
