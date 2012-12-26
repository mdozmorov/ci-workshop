source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
biocLite() # Install or update Bioconductor

# Preliminary steps: Image analysis and calculation of expression value
setwd("your/working/directory/with/CEL/files") # Set your working directory
# biocLite("affy") # Install affy package, if absent
biocLite("affy") # Install affy package
library(affy) # Load affy package
affy.data <- ReadAffy(celfile.path="../data/Su_CELs/")
write.table(exprs(affy.data),"Su_raw_matrix.txt", sep="\t") # Save unnormalized data
eset.mas5 <- mas5(affy.data) # Summarize with MAS5 algorithm
eset.mas5 # A summary of AffyBatch object
head(exprs(eset.mas5)) # Peek on the expression matrix
pData(eset.mas5) # Check experimental phenoData
boxplot(log2(exprs(eset.mas5))) # Boxplot of log2-transformed data
image(affy.data[1]) # Check image of the first CEL file
exprSet.nologs = exprs(eset.mas5) # Expression matrix
colnames(exprSet.nologs) # List the column (chip) names
# Rename the column names if we want 
colnames(exprSet.nologs) = c("brain.1", "brain.2", 
                             "fetal.brain.1", "fetal.brain.2",
                             "fetal.liver.1", "fetal.liver.2", 
                             "liver.1", "liver.2")
exprSet = log(exprSet.nologs, 2) # Log2-transformed expression
eset.rma <- justRMA(celfile.path="../data/Su_CELs/") # RMA summarization of the CEL files
biocLite("gcrma") # Install GCRMA package, if absent
library(gcrma) # Load GCRMA package
eset.gcrma <- justGCRMA(celfile.path="../data/Su_CELs/") # GCRMA summarization of the CEL files
eset.dChip = expresso(affy.data, normalize.method="invariantset", 
                      bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong") # dCHIP summarization
write.table(exprSet, file="Su_mas5_matrix.txt", quote=F, sep="\t")
# Run the Affy A/P call algorithm on the CEL files we processed above
data.mas5calls = mas5calls(affy.data)
# Get the actual A/P calls
data.mas5calls.calls = exprs(data.mas5calls)
# Print the calls as a matrix
write.table(data.mas5calls.calls, file="Su_mas5calls.txt", quote=F, sep="\t")

# Part I. Normalization of expression data
exprSetRaw = read.delim("Su_raw_matrix.txt",sep="\t") # Set separator as tab
trmean.col.1 = mean(exprSetRaw[,1], trim=0.02) # Trimmed mean in the first column
trmean = apply(exprSetRaw, 2, mean, trim=0.02) # Trimmed mean in all columns
trmean
sd = apply(exprSetRaw, 2, sd)
sd
median = apply(exprSetRaw, 2, median)
median
mean.of.trmeans = mean(trmean)
exprSet.trmean = exprSetRaw / trmean * mean.of.trmeans
boxplot(log2(exprSetRaw)) # Data before
boxplot(log2(exprSet.trmean)) # And after trimmed mean normalization
write.table(exprSet.trmean, file="Su_mas5_trmean_norm.txt", quote=F, sep="\t")
exprSet = exprSet.trmean
biocLite("limma")
library(limma)
exprSet.quantile = normalizeQuantiles(exprSet)
boxplot(log2(exprSet.quantile))
brain.fetalbrain.2color = read.maimages("../data/brain.fetalbrain.2color.data.txt", 
                                        columns=list(G="brain.1", R="fetal.brain.1", Gb="bg1", Rb="bg2"))
brain.fetalbrain.2color.loess = 
  normalizeWithinArrays(brain.fetalbrain.2color, method="loess")
par(mfrow=c(1,2)) # Set up a page with two figures next to each other
plotMA(brain.fetalbrain.2color) # Print the figures before
plotMA(brain.fetalbrain.2color.loess) # and after loess

# Part II. Probe annotation
biocLite("annotate")
library("annotate")
biocLite("hgu95av2.db") # Install, if not available
library("hgu95av2.db")
ls("package:hgu95av2.db") # Objects in this package
hgu95av2() # Description of an object
hgu95av2SYMBOL # Check specific object
get("39900_at",env=hgu95av2SYMBOL) # Accessing elements
hgu95av2SYMBOL[["39900_at"]]
hgu95av2SYMBOL$"39900_at"
symbols<-data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", ")) # Get all elements as data frame
symbols[c("1000_at","1001_at","1002_f_at","1003_s_at"),] # Get annotations for several elements
heasymbols<-as.list(hgu95av2SYMBOL) # Get all elements as list
length(symbols) # How many elements total
names(symbols)

# Part III. Accessing probe sequence data.
biocLite("hgu95av2probe")
library(hgu95av2probe)  # Opens library with probe sequence data.
print.data.frame(hgu95av2probe[1:10,]) # Prints probe sequences and their positions for first two Affy IDs
biocLite("Biostrings")
library(Biostrings)
pm <- DNAStringSet(hgu95av2probe$sequence)
names(pm) <- hgu95av2probe$Probe.Set.Name # Stores probe sequences as DNAStringSet object. See HT-Seq manual for more details.
head(pm) # Look what's there

# Part IV. Reading the NCBI's GEO microarray SOFT files in R/BioConductor
biocLite("GEOquery")
library(Biobase)
library(GEOquery)
#Download GDS file, put it in the current directory, and load it:
gds858 <- getGEO('GDS858', destdir=".")
#Or, open an existing GDS file (even if its compressed):
gds858 <- getGEO(filename='GDS858.soft.gz')
# Check metadata
Meta(gds858)$channel_count
Meta(gds858)$description
Meta(gds858)$feature_count
Meta(gds858)$platform
Meta(gds858)$sample_count
Meta(gds858)$sample_organism
Meta(gds858)$sample_type
Meta(gds858)$title
Meta(gds858)$type
# Check the expression data table
colnames(Table(gds858))
Table(gds858)[1:10,1:6]
# Convert to eset
eset <- GDS2eSet(gds858, do.log2=TRUE)
eset
featureNames(eset)[1:10]
sampleNames(eset)[1:10]
# Check pheno data
pData(eset)
pData(eset)$infection
pData(eset)$"genotype/variation"
# Check a subset of the data
eset["1320_at","GSM14504"]
# Get annotation data
Meta(gds858)$platform # Which platform?
#Download GPL file, put it in the current directory, and load it:
gpl96 <- getGEO('GPL96', destdir=".")
#Or, open an existing GPL file:
gpl96 <- getGEO(filename='GPL96.soft')
# What's in there
Meta(gpl96)$title
colnames(Table(gpl96))
Table(gpl96)[1:10,1:4]

# Cleanup, delete files
unlink(c("Su_raw_matrix.txt", "Su_mas5_matrix.txt", "Su_mas5calls.txt", "Su_mas5_trmean_norm.txt", "GDS858.soft.gz", "GPL96.soft"))
rm(list=ls()) # Remove all variables from the workspace