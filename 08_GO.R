################################################################################
################                  Gene Ontology                 ################
################################################################################

rm(list=ls())

##### Setting the working directory
setwd("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines")


##### Load the outcomes so far
load("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/go.RData")


###############################################################################
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
colnames(chrm.blast) <- c("chrm")


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
keytypes(GO.db)

## Select uniques genes in gen.id
gene.name <- unique(gen.id$go_id)

## GEt the Gene description and category of gene ontology
go.out <- AnnotationDbi::select(GO.db, keys=gene.name, keytype="GOID",
                         columns=c("DEFINITION", "TERM","ONTOLOGY"))

## Remove NA in go.out
go.out <- na.omit(go.out)
dim(go.out)

View(go.out)

## Merge gen.id and go.out by GO_id
gen.go <- merge(gen.id, go.out, by.x = "go_id", by.y = "GOID")

View(gen.go)

BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)

## Checking
columns(org.Gg.eg.db)
keytypes(org.Gg.eg.db)

## Select uniques genes
gene.sym <- unique(gen.id$external_gene_name)

## Remove empty ""
gene.sym <- gene.sym[gene.sym != ""]

## GEt the Gene description and category of gene ontology
go.gen.name <- AnnotationDbi::select(org.Gg.eg.db, keys=gene.sym, keytype="SYMBOL",
                         columns=c("GO", "GENENAME"))

View(go.gen.name)

## ADD GENENAME TO gen.go by GO_id
gen.go2 <- merge(gen.go, go.gen.name, by.x = "external_gene_name", by.y = "SYMBOL")
View(gen.go2)
colnames(gen.go2)

## Selecting columns external_gene_name, go_id, DEFINITION, TERM, ONTOLOGY.x, GENENAME
gen.go3 <- gen.go2[,c(1,2,6:12,15)]
View(gen.go3)
colnames(gen.go3)

length(unique(gen.go3$go_id))
unique(gen.go3$GENENAME)
unique(gen.go3$external_gene_name)


dim(gen.go3)
## REmove duplicated rows in all columns in gen.go3
gen.go3 <- gen.go3[!duplicated(gen.go3$TERM),]
dim(gen.go3)

## Sumarize the number of GO term per ONTOLOGY.x
table(gen.go3$ONTOLOGY.x)

## Cound number of GO term per each gene
table(gen.go3$external_gene_name)

## Write table to csv
write.csv(gen.go3, "C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/GO/gen.go3.csv", row.names = FALSE)



################################################################################


#### SAVING  ####
save.image("C:/Users/Daniel/Dropbox/Thesis/Molecular_Wrens/Radseq/Clines/go.RData")
