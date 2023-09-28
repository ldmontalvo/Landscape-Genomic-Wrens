
####################        DAPC and PCA for Wrens     #########################


rm(list=ls())

setwd("/home/ldmontalvo/PCA")
setwd("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/DAPC_PCA/HPC")


load("dapc_wren_n5nb.RData")


#install.packages("vcfR")
#install.packages("poppr")
#install.packages("ape")
#install.packages("RColorBrewer")
#install.packages("adegenet")
#install.packages("ggplot2")
#install.packages("ggforce")
#install.packages("rgdal")
#install.packages("ade4")
#install.packages("splancs")
#install.packages("spdep")
#install.packages("akima")
#install.packages("adespatial")


library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(adegenet)
library(ggplot2)
library(ggforce)
library(rgdal)
library(ade4)
library(splancs)
library(spdep)
library(akima)
library(adespatial)



#####  DATA MANAGEMENT  ######


# Calling the vcf file
wren.vcf.n5nb <- read.vcfR("./populations.snps.vcf")

# Read the pop file used in the denovo assembly in Stacks 
# I need this to get the order of samples
stack.ord <- read.table("./popmap_wren_n5nb_sh11.txt", col.names = c("Ind", "ord"))

# Adding a column for the order of samples
stack.ord$ord <- 1:nrow(stack.ord)

# Data with everything including labels and coordinates with info of morpho subspp
labels.all <- read.table("../R/Labels.field.csv", header=TRUE, sep=",", na.strings=TRUE, fileEncoding = "UTF-8-BOM")

# Merge the data set with coordinates and samples order
metadata.coor <- merge(stack.ord, labels.all, by="Ind")[,c(1,2,9,10)]
metadata.coor <- metadata.caoor[order(metadata.coor$ord,  decreasing = FALSE),]

## GETTING COORDINATES

# Getting coordinates and ordering them by samples
coor <- metadata.coor[order(metadata.coor$ord,  decreasing = FALSE), -c(1:2)]

# Making the coordinates numeric
coor2 <- cbind(as.numeric(as.character(coor$Long)), as.numeric(as.character(coor$Lat)))

# Reading the file with different levels of grouping including populations
wren_pop <- read.table("../R/wren_pop.txt", header=TRUE)

# Filtering the samples of C. brunnecapillus
wren_pop <- merge(wren_pop, stack.ord, by.x="Lev6_ind", by.y="Ind")

mia# Confirming if the names in vcf and pop data match
all(sort(colnames(wren.vcf.n5nb@gt)[-1]) == sort(wren_pop$Lev6_ind))

# We need to order the pop data Labels.field in the same order as the VCF file
# So we can assign names of populations
# GEt names in the order in the VCF
names <- colnames(wren.vcf.n5nb@gt)[-1]

# Order Labels.field in the same order as in VCF
Labels.field.n5nb <- wren_pop[match(names, wren_pop$Lev6_ind),]

#convert vcf to genlight to use with poppr or adegenet
gl.wren.n5b <- vcfR2genlight(wren.vcf.n5nb)

#Specifying the ploid of wrens as diploid
ploidy(gl.wren.n5b) <- 2

#Adding the pop label in the genlight object using Scientific Name from pop.data
pop(gl.wren.n5b) <- Labels.field.n5nb$Lev4_pop



#####  VISUALIZING THE PCA  ##### 

# Visualizing eigenvalues
wren.pca <- glPca(gl.wren.n5b, nf=2)
barplot(100*wren.pca$eig/sum(wren.pca$eig), col=heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percentage of variance Explained", line =2)
title(xlab = "Eigenvalues", line=1)

# Visualizing the PCA
wren.pca.scores <- as.data.frame(wren.pca$scores)
wren.pca.scores$pop <- pop(gl.wren.n5b)


# Making the palette of colors I will use
colr <- palette(brewer.pal(11, 'RdYlBu'))

# Adding one more color cause we need 12. RdYlBu just have 11
colr <- c(colr, "#0000A0")


# Plotting the PCA
#set.seed(9)
#ggplot(wren.pca.scores, aes(x=PC1, y=PC2)) +
#       geom_mark_ellipse(aes(color  = pop)) +  
#        geom_point(aes(color=pop), size=3) +
#        theme_bw() +
#        theme(legend.title = element_text(size=22)) +
#        theme(legend.text = element_text(size=20))


# Version with names in the ellipses
set.seed(9)
ggplot(wren.pca.scores, aes(x=PC1, y=PC2)) +
 geom_mark_ellipse(aes(color  = pop,
                               label  = pop),
                   expand = unit(0.5,"mm"),
                   label.buffer = unit(-5, 'mm')) +  
 geom_point(aes(color=pop)) +
 theme_bw() +
 theme(legend.title = element_text(size=22)) +
 theme(legend.text = element_text(size=20))  


# Different with labels in each point.
#set.seed(9)
#p <- ggplot(wren.pca.scores, aes(x=PC1, y=PC2, colour=pop))
#p <- p + geom_point(size=2)
#p <- p + stat_ellipse(level = 0.95, size = 1)
#p <- p + scale_color_manual(values = cols) 
#p <- p + geom_hline(yintercept = 0) 
#p <- p + geom_vline(xintercept = 0)
#p <- p + theme_bw()
#p <- p + theme(legend.position="none")
#p <- p + theme(legend.title = element_text(size=22))
#p <- p + geom_label(label=wren.pca.scores$pop)
#p




#####   sPCA (Spatial PCA)   #####

# Reading and converting STRUCTURE file into Genind object
wren.n5nb.genind <- import2genind("./populations2.str", onerowperind=FALSE, n.ind=112, n.loc=4409, col.lab=1, col.pop=2, ask=FALSE)


# Install rgdal to transform coordinates

# Setting the coordinates system as Long/Lat
cord.dec = SpatialPoints(cbind(coor2[,1], coor2[,2]), proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"))

# Converting to UTM, I'm using UTM zone 17s (epsg=32717)
# Most of the territory is in this zone
coor.utm <- spTransform(cord.dec, CRS("+init=epsg:32717"))
colnames(coor.utm@coords) <- c("X", "Y")

# Running the sPCA
cn.t5 <- chooseCN(xy=coor.utm@coords, ask=FALSE, type=5, d1=1, d2=2000)
cn.t6 <- chooseCN(xy=coor.utm@coords, type=6, k=10)
cn.t6.50 <- chooseCN(xy=coor.utm@coords, type=6, k=10)


# Running spatial PCA
wren_pca.t5 <- spca(obj=wren.n5nb.genind, cn=cn.t5)
wren_pca <- spca(obj=wren.n5nb.genind, cn=cn.t6)
wren_pca6.50 <- spca(obj=wren.n5nb.genind, cn=cn.t6.50)


# Plot the Spatial CA
#screeplot(wren_pca)
spca.t5 <- plot(wren_pca.t5)
spca.t6.50 <- plot(wren_pca6.50)

# Plotting PCA socres vs geography
pca.xy <- colorplot(wren_pca.t5, cex=3, main="sPCA for C. zonatus and C. fasciatus")
pca.xy6.50 <- colorplot(wren_pca6.50, cex=3, main="sPCA for C. zonatus and C. fasciatus")


##### Making Personalized sPCA plot  #####

# Check the components
head(wren_pca$li)
dim(wren_pca$li)

# Reading the names of Populations
poplab <- read.table("C:/Users/USER/Dropbox/Thesis/Molecular_Wrens/Wren_1/R/poplab.txt", header=FALSE, sep="\t")
colnames(poplab) <- c("sample", "pop")

# Get the components in a separate dataframe
pca.scores <- wren_pca.t5$li
pca.scores$sample <- rownames(pca.scores)
pca.scores <- merge(pca.scores, poplab, by="sample")

# Checking the pca scores
head(pca.scores)

# Writing the pca.scores for others analyses like the plumage coloration patterns
write.table(pca.scores[,-c(22,23)], "pca.scores.txt", sep="\t")

# creating a personalized Blue&Red gradient pallete for 16 colors
colr.3 <- c(PedroCarbo="#000080",LasGolondrinas="#0041C2", PVM="#1569C7", PatriciaPilar="#357EC7", 
            Chone="#3090C7", Calceta="#79BAEC", Montecristi="#7BCCB5", Machalilla="#4CC417",
            Manglareschurute="#B2C248", Arenillas="#EDE275", Cazaderos="#FFD801", Sullana="#FBB117", 
            Zapotillo="#FF6C00",Jaen="#FF6100", Chachapoyas="#FF3600", Celedin="#FF0000")


site.lat.ord <- c("PedroCarbo", "LasGolondrinas", "PVM", "PatriciaPilar", "Chone", "Calceta", "Montecristi", "Machalilla", 
                  "Manglareschurute", "Arenillas", "Cazaderos", "Zapotillo", "Sullana", "Jaen", "Chachapoyas", "Celedin")


# Version with names in the ellipses
set.seed(9)
ggplot(pca.scores, aes(x=pca.scores$`Axis 1`, y=pca.scores$`Axis 2`)) +
        geom_point(aes(color=factor(pop, level=site.lat.ord)), size=4) +
        geom_mark_ellipse(aes(color  = factor(pop, level=site.lat.ord),
                              label  = pop),
                          label.fontsize=8,
                          label.buffer = unit(-5, 'mm'),
                          expand = unit(0.2,"mm"),
                          label.fill=NA,
                          con.cap=0,
                          con.size=0.5) +  
        scale_color_manual(values = colr.3, name="Sampling \n Locations") +
        labs(y='PC 2', x = "PC 1") +
        theme_bw() +
        theme(legend.title = element_text(size=14)) +
        theme(legend.text = element_text(size=12)) 


#####   Multispati (With Adespatial)   #####


### Creating a listw object

# Creating the coordinates object
#coor.pca <- coor.utm@coords
#coor.pca <- as.data.frame(cbind(sample=stack.ord$Ind, coor.pca))
#coor.pca$X <- as.numeric(coor.pca$X)
#coor.pca$Y <- as.numeric(coor.pca$Y)
#coordinates(coor.pca) <- ~ X + Y

# Creating the Neighbor object
knea <- knearneigh(coordinates(coor.pca),k=10, longlat=FALSE)
neib <- knn2nb(knea)

# Creating the listw object
rn <- row.names(coor.pca)
all.linked <- max(unlist(nbdists(neib, coor.pca)))
listw.m2 <- dnearneigh(coor.pca, 10, all.linked, row.names=rn)
summary(listw.m2, coor.pca)

### Creating the dudi object, this means make the PCA
### I used samples (112) as rows and alleles on columns (12945)
dd1 <- dudi.pca(wren.n5nb.genind@tab, scannf = FALSE)

# Running the Multispati function
ms.pca1 <- multispati(dd1, listw.m2, scannf=TRUE, nfposi=2, nfnega=2)



#####   Other Analyses   #####


### Global test
Gtest <- global.rtest(wren.n5nb.genind@tab, wren_pca$lw, nperm=9999)
Gtest
plot(Gtest)


Gtest6.50 <- global.rtest(wren.n5nb.genind@tab, wren_pca6.50$lw, nperm=9999)
Gtest6.50
plot(Gtest6.50)


### Local test
Ltest <- local.rtest(wren.n5nb.genind@tab, wren_pca$lw, nperm=9999)
Ltest
plot(Ltest)


Ltest6.50 <- local.rtest(wren.n5nb.genind@tab, wren_pca6.50$lw, nperm=9999)
Ltest6.50
plot(Ltest6.50)


head(wren_pca6.50$c1)[,1:3]



### Interpolate Component Scores

x <- coor.utm@coords[,1]
y <- coor.utm@coords[,2]

# Interpolation using k=10
int.pca <- interp(x, y, wren_pca$li[,1], duplicate="mean")
image(int.pca, col=azur(100))
points(x, y)

mypal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
        title("sPCA - interpolated map of individual scores")
        points(x,y)
}
filled.contour(int.pca, color.pal=mypal, nlev=50,
               key.title=title("lagged \nscore 1"), plot.title=annot())


# Interpolation using k=50
int.pca6.50 <- interp(x, y, wren_pca6.50$li[,1], duplicate="mean")
image(int.pca6.50, col=azur(100))
points(x, y)

filled.contour(int.pca6.50, color.pal=mypal, nlev=50,
               key.title=title("lagged \nscore 1"), plot.title=annot())







#### SAVING ####
save.image("dapc_wren_n5nb.RData")

