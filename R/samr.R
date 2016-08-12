library(affy) # Load affy package 
library(GEOquery)
#Download GDS file, put it in the current directory, and load it:
gds858 <- getGEO('GDS858', destdir=".") 
gds858<-getGEO(filename="GDS858.soft.gz") # If FTP doesn't work, read in from local file
eset <- GDS2eSet(gds858, do.log2=TRUE) # Convert the data to ESET object
help(ExpressionSet) # If needed, refresh your memory
pData(eset)$infection # Let's check the infection status
source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
biocLite("samr")
library(samr) # Load the library
library(limma) # To get access to normalizeQuantiles
table(pData(eset)$infection) # How many samples are in each infection status
selected<-grep("uninfected|FRD440",pData(eset)$infection) # Select indexes
pData(eset)$infection[selected] # Check if the right infection status was selected
y=c(rep(1,4),rep(2,4)) # Vector of outcome measurements
y # Visualize it
exprs.selected<-exprs(eset)[,selected] # Get a subset of expression values of the selected samples
exprs.selected.q<-normalizeQuantiles(exprs.selected) # Quantile normalize the data
genenames<-rownames(exprs(eset)) # Get row names = gene names or IDs
data<-list(x=exprs.selected.q,y=y, geneid=genenames,genenames=genenames,logged2=T)
samr.obj<-samr(data,resp.type="Two class unpaired",nperms=100)
names(samr.obj) # Get object names from samr.obj
delta.table<-samr.compute.delta.table(samr.obj, min.foldchange=1.5) # Compute thresholds for different deltas
delta.table # Look at the whole range
delta.table[delta.table[,"median FDR"] < 0.1,][1,] # Check delta corresponding to median FDR ~0.1
delta<-1.5 # Select the delta
samr.plot(samr.obj,delta) # Check SAM plot
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta,data,delta.table,min.foldchange=1.5) # Summarize significant genes
names(siggenes.table) # What data we have in the summary list
nrow(siggenes.table$genes.up) # How many upregulated genes
nrow(siggenes.table$genes.lo) # How many downregulated
# Or
siggenes.table$ngenes.up; siggenes.table$ngenes.lo
siggenes.table$genes.up # Check how table with the results look like
# Write the results in file
write.table(siggenes.table$genes.up,"genes.up.txt",sep='\t') 
write.table(siggenes.table$genes.lo,"genes.dn.txt",sep='\t')
# Extract up- and downregulated IDs
up.ids<-siggenes.table$genes.up[,"Gene ID"] # Upregulated
dn.ids<-siggenes.table$genes.lo[,"Gene ID"] # Downregulated
# Annotate Affy probe IDs
Meta(gds858)$platform # Check which platform do we have
gpl96<-getGEO('GPL96', destdir=".") # Get the data for this platform from GEO
Table(gpl96)[1:5,1:15] # Check the content of annotation data
Table(gpl96)[Table(gpl96)[,"ID"] %in% up.ids,c("ID","Gene Symbol","Gene Title")] # Extract annotation for up.ids
# Use biomaRt for annotation
biocLite("biomaRt") 
library(biomaRt)
listDatasets(useMart("ensembl")) # List available datasets
mart<-useMart("ensembl", dataset="hsapiens_gene_ensembl") # Load BIOMART dataset for homo sapiens
# Information - lists of filters and attributes
head(listFilters(mart), n=50) # Filters, these are our IDs we'll be subsetting the BIOMART annotations on
head(listAttributes(mart), n=50) # Attributes, these are annotations we would like to get for our IDs
Meta(gpl96)$title # Check which microarray we have, to selext the right attributes
attr<-listAttributes(mart) # Get all attributes as a table
attr[grep("affy",attr[,1]),] # Parse them for anything that looks like from affymetrix
# Get annotations from Biomart
genes.up<-getBM(attributes=c('affy_hg_u133a','external_gene_id','description'), filters='affy_hg_u133a', values=up.ids, mart=mart)#, uniqueRows=T)
# Clean workspace
unlink(c("GDS858.soft.gz","genes.up.txt","genes.dn.txt","GPL96.soft"))
