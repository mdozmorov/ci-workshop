###################################################
# Differential expression analysis
# Mikhail Dozmorov
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

#####################
# SAM (Significance Analysis of Microarrays)
##################### 

In this exercise we will be detecting differentially expressed genes between uninfected lung epithelial Calu-3 cells and cells infected with FRD440 strain of Pseudomonas aeruginosa.

Let's load necessary libraries, and the data themselves.

> library(affy) # Load affy package 
> library(GEOquery)
> #Download GDS file, put it in the current directory, and load it:
> gds858 <- getGEO('GDS858', destdir=".") 
> gds858<-getGEO(filename="GDS858.soft.gz") # If FTP doesn't work, read in from local file
> eset <- GDS2eSet(gds858, do.log2=TRUE) # Convert the data to ESET object
> help(ExpressionSet) # If needed, refresh your memory
> pData(eset)$infection # Let's check the infection status
> source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
> biocLite("samr")
> library(samr) # Load the library
> library(limma) # To get access to normalizeQuantiles

Let's check which infection statuses do we have

> table(pData(eset)$infection) # How many samples are in each infection status

Select indexes of the columns associated with "uninfected" or "FRD440" infection status. Check that the right infection status was selected.

> selected<-grep("uninfected|FRD440",pData(eset)$infection) # Select indexes
> pData(eset)$infection[selected] # Check if the right infection status was selected

We selected 4 "uninfected" samples and 4 "FRD440" samples. Let's make a vector of outcome measurements, defining groups of replicates. "Uninfected" samples are labeled as "1" and "FRD440" samples are labeled as "2".

> y=c(rep(1,4),rep(2,4)) # Vector of outcome measurements
> y # Visualize it

From the whole expression matrix we select a subset of columns corresponding to the selected indexes. We then quantile normalize the data.

> exprs.selected<-exprs(eset)[,selected] # Get a subset of expression values of the selected samples
> exprs.selected.q<-normalizeQuantiles(exprs.selected) # Quantile normalize the data

In order to perform the analysis SAM requires a list object containing a data object with expression values in a form of p by n matrix, one observation per column (missing values allowed),  a vector of length n of outcome measurements, vectors of gene names and gene IDs, both of length p, and a boolean indicating whether the data is log2-transformed. This object resembles ExpressionSet with slots having different data, but we assemble it as a list.
Now we have our expression matrix (exprs.selected.q), our vector of outcomes (y). Let's use row names of the expression matrix as both gene names and IDs, and assemble our data.

> genenames<-rownames(exprs(eset)) # Get row names = gene names or IDs
> data<-list(x=exprs.selected.q,y=y, geneid=genenames,genenames=genenames,logged2=T)

SAM can perform a variety types of analyses, specified in "resp.type" parameter. Each analysis type requires specific formatting of outcome vector and expression data. Refer to help(samr) for the details. Now, we're performing "Two class unpaired" analysis.

> samr.obj<-samr(data,resp.type="Two class unpaired",nperms=100)

Let's check what we have. Everything obvious, isn't it?

> names(samr.obj) # Get object names from samr.obj

Now we have to choose the delta value that is able to give the best compromise in terms of called genes, false genes and False Discovery Rate (FDR). In microarray analysis is very important to have statistically robust results, but we have keep in mind that too small sized results are not able to describe the biological meaning of the experiment. In any case, keeping the FDR < 10% (the number of false positives is < 10%) is pretty safe in most of the cases.
In general, defining the cut-off is a subjective choice and there is no absolute best way to do it.

> delta.table<-samr.compute.delta.table(samr.obj, min.foldchange=1.5) # Compute thresholds for different deltas
> delta.table # Look at the whole range

Let's select delta with median FDR <10% - subset the whole delta table and take the first row.

> delta.table[delta.table[,"median FDR"] < 0.1,][1,] # Check delta corresponding to median FDR ~0.1

How many genes are called significant? What's the median number of false positives? Is it ~10%, as expected?

Let's select the delta, and have a look at SAM plot.

> samr.plot(samr.obj,delta) # Check SAM plot

Do we have larger number of upregulated genes? Or downregulated? Let's have a look at them.

> siggenes.table<-samr.compute.siggenes.table(samr.obj,delta,data,delta.table,min.foldchange=1.5) # Summarize significant genes
> names(siggenes.table) # What data we have in the summary list
> nrow(siggenes.table$genes.up) # How many upregulated genes
> nrow(siggenes.table$genes.lo) # How many downregulated
# Or
> siggenes.table$ngenes.up; siggenes.table$ngenes.lo

Let's have a look at actual differentially expressed genes

> siggenes.table$genes.up # Check how table with the results look like

We can write them in a clipboard, then paste into Excel. Or export them to files.

# Write the results in clipboard
> write.table(siggenes.table$genes.up,"genes.up.txt",sep='\t') 
> write.table(siggenes.table$genes.lo,"genes.dn.txt",sep='\t')

These results show pretty much everything we need to know. But what about getting annotations to the Affy probe IDs? Let's extract those IDs and annotate them.

# Extract up- and downregulated IDs
> up.ids<-siggenes.table$genes.up[,"Gene ID"] # Upregulated
> dn.ids<-siggenes.table$genes.lo[,"Gene ID"] # Downregulated

We can simply get the platform data from GEO, and extract subsets of IDs.

# Annotate Affy probe IDs
> Meta(gds858)$platform # Check which platform do we have
> gpl96<-getGEO('GPL96', destdir=".") # Get the data for this platform from GEO
> Table(gpl96)[1:5,1:15] # Check the content of annotation data
> Table(gpl96)[Table(gpl96)[,"ID"] %in% up.ids,c("ID","Gene Symbol","Gene Title")] # Extract annotation for up.ids

Or we can use biomaRt.

# Use biomaRt for annotation
> biocLite("biomaRt") 
> library(biomaRt)
> listDatasets(useMart("ensembl")) # List available datasets
> mart<-useMart("ensembl", dataset="hsapiens_gene_ensembl") # Load BIOMART dataset for homo sapiens
# Information - lists of filters and attributes
> head(listFilters(mart), n=50) # Filters, these are our IDs we'll be subsetting the BIOMART annotations on
> head(listAttributes(mart), n=50) # Attributes, these are annotations we would like to get for our IDs
> Meta(gpl96)$title # Check which microarray we have, to selext the right attributes
> attr<-listAttributes(mart) # Get all attributes as a table
> attr[grep("affy",attr[,1]),] # Parse them for anything that looks like from affymetrix

Now we know that our IDs correspond to "affy_hg_u133a" attributes in Biomart. Let's extract annotations for them.

# Get annotations from Biomart
> genes.up<-getBM(attributes=c('affy_hg_u133a','external_gene_id','description'), filters='affy_hg_u133a', values=up.ids, mart=mart)#, uniqueRows=T)

Do that for the downregulated genes, and write the results to file.

Clean your workspace

> unlink(c("GDS858.soft.gz","genes.up.txt","genes.dn.txt","GPL96.soft"))

############
# Exercise 1
############

Using ALL dataset, find differentially expressed genes by limma and SAM. Start by importing and exploring the data.

> source("http://www.bioconductor.org/biocLite.R")
> biocLite("ALL") # Install ALL package
> library("ALL") # Load ALL package
> data("ALL") # Attach ALL data to the workspace
> ALL # Check what's in there

 got expression matrix with  128 samples and 12625 genes, and associated phenotypic data. We can looks at the parameters available to us in the phenotypic data, and the results of molecular biology testing for the 128 samples.

> names(pData(ALL)) # Check what phenotype data we have
> ALL$mol.biol # Results of molecular biology testing for the 128 samples

Ignoring the samples which came back negative on this test (NEG), most have been classified as having a translocation between chromosomes 9 and 22 (BCR/ABL), or a translocation between chromosomes 4 and 11 (ALL1/AF4).

For the purposes of this example, we are only interested in these two subgroups, so we will create a filtered version of the dataset using this as a selection criteria:

> eset <- ALL[, ALL$mol.biol %in% c("BCR/ABL", "ALL1/AF4")] # Select a subset of the data

The resulting variable, eset, contains just 47 samples - each with the full 12,625 gene expression levels.

Use eset@mol.bio as a factor for building matrix for lima, and for defining classes for SAM.

