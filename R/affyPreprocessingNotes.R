# http://jura.wi.mit.edu/bio/education/bioinfo2007/arrays/array_exercises_1R.html
# http://bioconductor.org/packages/release/bioc/vignettes/annotate/inst/doc/annotate.pdf
# http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#biocon_affypack
# http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/

Microarray data retrieval and preprocessing

Introduction

You'll be using a sample of expression data from a study using Affymetrix (one color) U95A arrays that were hybridized to tissues from fetal and human liver and brain tissue. Each hybridization was performed in duplicate. Many other tissues were also profiled but won't be used for these exercises.

What we'll be doing to analyze these data:

preprocess the raw data (summarize probe measurements into one measurement for each probe)
normalize data from these eight chips
calculate Absent/Present calls which attempts to label genes that are "expressed"
calculate expression ratios of genes between two different tissues
use a common statistical test to identify differentially expressed genes
flag low intensity data (most probably background noise)
cluster a differentially expressed subset of all genes to identify those with similar expression profiles
try to find what functions specific groups of genes (with similar expression profiles) have in common

Preliminary steps: Image analysis and calculation of expression value

As described in Su et al., 2002 [https://www.ncbi.nlm.nih.gov/pubmed/11904358?dopt=Abstract], human tissue samples were hybridized on Affymetrix (one-color) arrays and chips were scanned. For each tissue, at least two independent samples were hybridized to separate chips.
Scanned images were quantified (including measurement of background) using standard software.

Part 0. Preprocessing and normalization of Affymetrix expression data

Data from a probeset (a series of oligos designed to a specific gene target) needs to be summarized by calculating an expression values for that probeset. This can be done using a Bioconductor/R version of the methods in the Microarray Suite 5.0 (MAS5, the standard Affymetrix algorithm) used by the Affymetrix software suite.
See the manuals from Affymetrix[http://www.affymetrix.com/support/technical/manuals.affx] for more information about these processes, and the Statistical Algorithms Description Document[http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf] for the actual equations used.
You also have the option of performing other algorithms such as RMA, GCRMA, and dChip.
Download the starting data ("CEL" files) and unzip this file to get 8 files ending in ".CEL".
Start R on your laptop by clicking on the "R" icon. 
If you're using your own computer, install Bioconductor (if you haven't done so before) by starting R and then pasting or typing these commands. The installation will take a few minutes.

> source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
> biocLite() # Install or updated Bioconductor

Use the pull-down menu (File >> Change dir [Windows] or Misc >> Change working directory [Mac]) to go to the directory where you put the raw data (CEL files). Or use the following command to type in path to CEL files manually

> setwd("your/working/directory/with/CEL/files") # Set your working directory

Install, if necessary, and load the "library" that contains the Affymetrix microarray code we'll want to use with the command
                                                         
> biocLite("affy") # Install affy package, if absent
> library(affy) # Load affy package

Read the CEL files (first command below) and then summarize and normalize with MAS5 (second command below). This could take a few minutes.

> affy.data = ReadAffy(celfile.path="Su_CELs/")
> eset.mas5 = mas5(affy.data) # Summarize with MAS5 algorithm

The variable 'eset.mas5' contains normalized expression values for all probesets, along with other information. Let's check its parameters

> eset.mas5 # A summary of AffyBatch object
> head(exprs(eset.mas5)) # Peek on the expression matrix
> pData(eset.mas5) # Check experimental phenoData

The phenodata object contains all the information about the samples of our dataset. It is always important to keep track off all the information about our samples, since they can strongly affect the final results.

It is possible to visually explore the dataset using different methods.
The boxplot visualizes the distribution of the intensities in each array of the dataset.

> boxplot(log2(exprs(eset.mas5))) # Boxplot of log2-transformed data

It is also possible to visualize the image representing each CEL file (in this example, we visualize the first slide of the dataset only). 
In this way, it is possible to pinpoint technical problems occurring eventually only to one region of the array.

> image(affy.data[1]) # Check image of the first CEL file

Let's continue by getting the expression matrix (probesets/genes in rows, chips in columns).

> exprSet.nologs = exprs(eset.mas5) # Expression matrix
> colnames(exprSet.nologs) # List the column (chip) names
# Rename the column names if we want
> colnames(exprSet.nologs) = c("brain.1", "brain.2", 
                             "fetal.brain.1", "fetal.brain.2",
                             "fetal.liver.1", "fetal.liver.2", 
                             "liver.1", "liver.2")

At this time let's log-transform the expression values to get a more normal distribution. We have to remember we've done this when we calculate ratios. Logarithms can use any base, but base 2 is easiest when transforming ratios, since transformed 2-fold ratios up or down will be +1 or -1. As a result, we'll do all logs with base 2 to keep thing simplest.

> exprSet = log(exprSet.nologs, 2) # Log2-transformed expression

Optional. In place of mas5 normalization above, we could choose many other ways to preprocess the data, such as RMA (which outputs log2-transformed expression values)

> eset.rma <- justRMA(celfile.path="Su_CELs/") # RMA summarization of the CEL files

or GCRMA (which also outputs log2-transformed expression values)

> biocLite("gcrma") # Install GCRMA package, if absent
> library(gcrma) # Load GCRMA package
> eset.rma <- justGCRMA(celfile.path="Su_CELs/") # GCRMA summarization of the CEL files

or dChip (also known as MBEI; not log-transformed)
eset.dChip = expresso(affy.data, normalize.method="invariantset", 
                      bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong" # dCHIP summarization
					  
To print out our expression matrix (as with most data), we can use a command like

> write.table(exprSet, file="Su_mas5_matrix.txt", quote=F, sep="\t")

to get a tab-delimited file that we could view in Excel or a text editor.
While we're doing Affymetrix-specific preprocessing, let's calculate an Absent/Present call for each probeset.

# Run the Affy A/P call algorithm on the CEL files we processed above
> data.mas5calls = mas5calls(affy.data)
# Get the actual A/P calls
> data.mas5calls.calls = exprs(data.mas5calls)
# Print the calls as a matrix
> write.table(data.mas5calls.calls, file="Su_mas5calls.txt", quote=F, sep="\t")

Other methods of normalization available: ?justPlier, ?expresso

Part I. Normalization of expression data

Why normalize? Chips may have been hybridized to different amounts of RNA, for different amounts of time, with different batches of solutions, etc. Normalization should remove systematic biases and make any comparisons between chips more meaningful.
If you performed Part 0 above, your data is already normalized. But if you got an unnormalized expression from another source, how could you normalize it?
Download Su_raw_matrix.txt, an unnormalized matrix of the same Su data.
Read an expression matrix in the form of a tab-delimited text file that you created or downloaded above. The first line contains column headers, and the first row (without a header field) contains gene identifiers.

> exprSetRaw = read.delim("Su_raw_matrix.txt",sep="\ ") # Set separator as space

Calculate the trimmed mean of all expression values on each chip.
A trimmed mean calculates a summary value that is somewhere between the mean and median of the set of values.
We remove the top and bottom 2% of values (for example), and find the mean of the remaining values,
To get the trimmed mean of column 1 on our matrix, we can use a command like

> trmean.col.1 = mean(exprSetRaw[,1], trim=0.02) # Trimmed mean in the first column

We can apply this command to all columns at once with the "apply" command, where the '2' means to run the command on all columns (and '1' would be to do the same thing on all rows).

> trmean = apply(exprSetRaw, 2, mean, trim=0.02) # Trimmed mean in all columns
> trmean

If we're curious about the variation between genes on each chip, we can find the standard deviation for each chip

> sd = apply(exprSetRaw, 2, sd)
> sd

Or we could compare medians

> median = apply(exprSetRaw, 2, median)
> median

Divide all expression values by the mean for that chip, and multiply by a scaling factor, such as the mean of the trimmed means.

> mean.of.trmeans = mean(trmean)
> exprSet.trmean = exprSetRaw / trmean * mean.of.trmeans

Let's look at the data before and after our normalization

> boxplot(log2(exprSetRaw)) # Data before
> boxplot(log2(exprSet.trmean)) # And after trimmed mean normalization

Print the data that's been normalized by the global trimmed mean

> write.table(exprSet.trmean, file="Su_mas5_trmean_norm.txt", quote=F, sep="\t")

For further processing, rename your data

> exprSet = exprSet.trmean

Quantile normalization is often too extreme, but it's common for Affymetrix probe-level normalization (being part of MAS5). If you wish to use it, one way is to use a command from the "limma" package

> library(limma)
> exprSet.quantile = normalizeQuantiles(exprSet)

And visualize our work

> boxplot(log2(exprSet.quantile))

Another common normalization choice (generally for 2-color arrays) is loess. We combined the "brain.1" and "fetal.brain.1" expression data into a file that we can pretend is 2-color data (together with 0 values for "background"). Then we can normalize it with loess.
Start by downloading brain.fetalbrain.2color.data.txt into your working directory.
Read this fake 2-color file, indicating the columns for the Red and Green channels and the Red and Green backgrounds.

> brain.fetalbrain.2color = read.maimages("brain.fetalbrain.2color.data.txt", 
	columns=list(G="brain.1", R="fetal.brain.1", Gb="bg1", Rb="bg2"))

Now we have the two channels of data which we can loess normalize

> brain.fetalbrain.2color.loess = 
	normalizeWithinArrays(brain.fetalbrain.2color, method="loess")

How do you know if anything changed? Print an MA plot before and after normalization

> par(mfrow=c(1,2)) # Set up a page with two figures next to each other
> plotMA(brain.fetalbrain.2color) # Print the figures
> plotMA(brain.fetalbrain.2color.loess) # and after loess

How do the figures differ?

Part II. Probe annotation

Any given experiment typically involves a set of known identifers (probes
in the case of a microarray experiment). These identifers are typically unique
(for any manufacturer). But having Affymetrix IDs ("1001_at") as probe identifiers is not very helpful when we want to know which gene a probe is actually assesses. There are numerous annotation packages available on Bioconductor at http://www.bioconductor.org/packages/release/data/annotation/. We will consider the AFFymetrix human gene chip, hgu95av2, for our example.
We first load this chip's package and annotate.

> library("annotate")
> biocLite("hgu95av2.db") # Install, if not available
> library("hgu95av2.db")
> ls("package:hgu95av2.db") # Objects in this package

We see the listing of twenty or thirty different R objects in this package. Most
of them represent mappings from the identifiers on the Affymetrix chip to the
different biological resources and you can find out more about them by using the
R help system, since each has a manual page that describes the data together
with other information such as where, when and what files were used to construct
the mappings. Also, each meta-data package has one object that has the same
name as the package basename, in this case it is hgu95av2. This is function
and it can be invoked to find out some of the different statistics regarding the
mappings that were done. Its manual page lists all data resources that were
used to create the meta-data package.

> hgu95av2() # Description of an object

Now let's consider a specific object, say the hgu95av2SYMBOL object.

> hgu95av2SYMBOL # Check specific object

If we type its name we see that it is an R environment; all this means is that it
is a special data type that is efficient at storing and retrieving mappings between
symbols (the Affymetrix identifiers) and associated values (the gene symbols).
We can retrieve values from this environment in many different ways. Suppose that we are interested in finding the gene symbol for the Affymetrix probe, 39900_at. Then we can do it in any one of the following ways:

> get("39900_at",env=hgu95av2SYMBOL) # Accessing elements
> hgu95av2SYMBOL[["39900_at"]]
> hgu95av2SYMBOL$"39900_at"

If you want to get more than one object from an environment you also have
a number of choices. You can convert the object infto a data frame, and subset if by your set of IDs. You can extract them one at a time using either a for
loop or apply or eapply. It will be more efficient to use mget since it does the
retrieval using internal code and is optimized. You may also want to turn the
environment into a named list, so that you can perform different list operations
on it, this can be done using the function contents or as.list.

> symbols<-data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", ")) # Get all elements as data frame
> symbols[c("1000_at","1001_at","1002_f_at","1003_s_at"),] # Get annotations for several elements
> symbols<-as.list(hgu95av2SYMBOL) # Get all elements as list
> length(symbols) # How many elements total
> names(symbols)[1:15]

However, it is recommended to do as the last step, using manufacturer's IDs throughout the analysis.

Part III. Accessing probe sequence data.

> biocLite("hgu95av2probe")
> library(hgu95av2probe)  # Opens library with probe sequence data.
> print.data.frame(hgu95av2probe[1:10,]) # Prints probe sequences and their positions for first two Affy IDs
> biocLite("Biostrings")
> library(Biostrings)
> pm <- DNAStringSet(hgu95av2probe$sequence)
> names(pm) <- hgu95av2probe$Probe.Set.Name # Stores probe sequences as DNAStringSet object. See HT-Seq manual for more details.
> head(pm) # Look what's there

Part IV. Reading the NCBI's GEO microarray SOFT files in R/BioConductor

GEO expression omnibus, https://www.ncbi.nlm.nih.gov/geo/, is the largest repository of gene expression data. It may have the data you're looking for, saving you time and money to do the experiment yourself. 
We will discuss how to load GEO SOFT format microarray data from the Gene Expression Omnibus database (GEO) (hosted by the NCBI) into R/BioConductor. SOFT stands for Simple Omnibus Format in Text. There are actually four types of GEO SOFT file available:
GEO Platform (GPL)
These files describe a particular type of microarray. They are annotation files.
GEO Sample (GSM)
Files that contain all the data from the use of a single chip. For each gene there will be multiple scores including the main one, held in the VALUE column.
GEO Series (GSE)
Lists of GSM files that together form a single experiment.
GEO Dataset (GDS)
These are curated files that hold a summarised combination of a GSE file and its GSM files. They contain normalised expression levels for each gene from each sample (i.e. just the VALUE field from the GSM file).
As long as you just need the expression level then a GDS file will suffice. If you need to dig deeper into how the expression levels were calculated, you'll need to get all the GSM files instead (which are listed in the GDS or GSE file)

Installing GEOquery
Assuming you are running a recent version of BioConductor (1.8 or later) you should be able to install it from within R as follows:

> source("http://www.bioconductor.org/biocLite.R")
> biocLite("GEOquery")

Loading a GDS file with GEOquery
Here is a quick introduction to how to load a GDS file, and turn it into an expression set object:

> library(Biobase)
> library(GEOquery)
#Download GDS file, put it in the current directory, and load it:
> gds858 <- getGEO('GDS858', destdir=".")
#Or, open an existing GDS file (even if its compressed):
> gds858 <- getGEO(filename='GDS858.soft.gz')

The SOFT file is available in compressed form here GDS858.soft.gz, but GEOquery takes care of finding this file for you and unzipping it automatically.
There are two main things the GDS object gives us, meta data (from the file header) and a table of expression data. These are extracted using the Meta and Table functions. First lets have a look at the metadata:

> Meta(gds858)$channel_count
> Meta(gds858)$description
> Meta(gds858)$feature_count
> Meta(gds858)$platform
> Meta(gds858)$sample_count
> Meta(gds858)$sample_organism
> Meta(gds858)$sample_type
> Meta(gds858)$title
> Meta(gds858)$type

Useful stuff, and now the expression data table:

> colnames(Table(gds858))
> Table(gds858)[1:10,1:6]

Now, lets turn this GDS object into an expression set object (using base 2 logarithms) and have a look at it:

> eset <- GDS2eSet(gds858, do.log2=TRUE)
> eset
> featureNames(eset)[1:10]
> sampleNames(eset)[1:10]

GEOquery does an excellent job of extracting the phenotype data, as you can see:

> pData(eset)
> pData(eset)$infection
> pData(eset)$"genotype/variation"

As with any expression set object, its easy to pull out a subset of the data:

> eset["1320_at","GSM14504"]

In addition to loading a GDS file to get the expression levels, you can also load the associated platform annotation file. You can find this out from the GDS858 meta information:

> Meta(gds858)$platform # Which platform?

Now let's load up the GPL file and have a look at it (its a big file, about 12 MB, so this takes a while!):

#Download GPL file, put it in the current directory, and load it:
> gpl96 <- getGEO('GPL96', destdir=".")
#Or, open an existing GPL file:
> gpl96 <- getGEO(filename='GPL96.soft')

As with the GDS object, we can use the Meta and Table functions to extract information:

> Meta(gpl96)$title
> colnames(Table(gpl96))

Lets look at the first four columns, for the first ten genes:

> Table(gpl96)[1:10,1:4]

You can also use Bioconductor annotation file for GPL96/HG-U133A by using library(hgu133a).


