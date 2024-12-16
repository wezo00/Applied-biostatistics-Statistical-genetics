# install bioconductor packages
#source ("http://bioconductor.org/biocLite.R")
#list.of.packages <- c("snpStats", "SNPRelate","rtracklayer", "biomaRt")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) biocLite(new.packages)


#list.of.packages <- c('plyr', 'LDheatmap', 'doParallel', 'ggplot2', 'coin' ,'igraph', 'devtools', 'downloader')
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# GenABEL has moved to CRAN archive. The below command for local installation from CRAN archive.
#install.packages("https://cran.r-project.org/src/contrib/Archive/GenABEL/GenABEL_1.7-6.tar.gz", 
#                type = "source", repos = NULL)

#library(devtools)
#install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")

### STEP 1

library(snpStats)

# Read in PLINK files
geno <- read.plink("GWAStutorial.bed", "GWAStutorial.bim", "GWAStutorial.fam", na.strings = ("-9"))

# Obtain the SnpMatrix object (genotypes) table from geno list
# Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
print(genotype) 

#Obtain the SNP information from geno list
genoBim <- geno$map
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))
# Remove raw file to free up memory
rm(geno)

# Read in clinical file
clinical <- read.csv("GWAStutorial_clinical.csv",colClasses=c("character", "factor", "factor", rep("numeric", 4)))
rownames(clinical) <- clinical$FamID
print(head(clinical))
# Subset genotype for subject data
genotype <- genotype[clinical$FamID, ]
print(genotype)  # Tutorial: All 1401 subjects contain both clinical and genotype data

### STEP 2

# Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print(head(snpsum.col))

# Setting thresholds
call <- 0.95
minor <- 0.01

# Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),
    "SNPs will be removed due to low MAF or call rate.\n") #203287 SNPs will be removed
# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]
print(genotype)     
#save(genotype, snpsum.col, genoBim, clinical, file=working.data.fname(2))
library(snpStats)
library(SNPRelate)               # LD pruning, relatedness, PCA
library(plyr)

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

head(snpsum.row)

# Setting thresholds
sampcall <- 0.95    # Sample call rate cut-off
hetcutoff <- 0.1    # Inbreeding coefficient cut-off

sampleuse <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall & abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE    # remove NA's as well
cat(nrow(genotype)-sum(sampleuse), 
    "subjects will be removed due to low sample call rate or inbreeding coefficient.\n") #0 subjects removed

# Subset genotype and clinical data for subjects who pass call rate and heterozygosity crtieria
genotype <- genotype[sampleuse,]
clinical<- clinical[ rownames(genotype), ]

### STEP 3

# Checking for Relatedness

ld.thresh <- 0.2    # LD cut-off
kin.thresh <- 0.1   # Kinship cut-off

# Create gds file, required for SNPRelate functions
snpgdsBED2GDS("GWAStutorial.bed", "GWAStutorial.bim", "GWAStutorial.fam", "GWAStutorial.gds")

genofile <- openfn.gds("GWAStutorial.gds", readonly = FALSE)

# Automatically added "-1" sample suffixes are removed
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

#Prune SNPs for IBD analysis
set.seed(1000)
geno.sample.ids <- rownames(genotype)

### PROBLEM
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh, sample.id = geno.sample.ids, snp.id = colnames(genotype))
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")  # Tutorial: expect 72812 SNPs


