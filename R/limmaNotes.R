###################################################
# Differential expression analysis
# Mikhail Dozmorov
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

#####################
# Limma
##################### 

Limma is a software package for the analysis of gene expression microarray data, especially the use of linear models for analyzing designed experiments and the assessment of differential expression. The package includes pre-processing capabilities for two-color spotted arrays. The differential expression methods apply to all array platforms and treat Affymetrix, single channel and two channel experiments in a unified way. The methods are described in Smyth 2004 [https://www.ncbi.nlm.nih.gov/pubmed/15297296] and in the limma manual [http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf]. An illustrated introduction for the GUI packages can be found at WEHI [http://bioinf.wehi.edu.au/limma/index.html].

We will be analyzing tissue-specific differences. Prepare the data:

> source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
> biocLite("limma")
> library(limma)
> limmaUsersGuide() # Opens pdf manual for limma
> library(affy)
> eset.rma <- justRMA(celfile.path="../../Exercises/AffyPreprocessing/Su_CELs/") # RMA summarization of the CEL files
> pData(eset.rma) # Check what samples we have


There are two different ways to form the design matrix. We can either
1. create a design matrix which includes a coefficient for the treated vs wild type difference, or
2. create a design matrix which includes separate coefficients for wild type and mutant mice and then extract the difference as a contrast.

For the first approach, the treatment-contrasts parametrization, the design matrix should be as follows:

# Design matrix: Treatment-constrast parametrization
> a<-rep(0,length(pData(eset.rma)$sample)) # Create a vector of 0
> a[grep("liver", rownames(pData(eset.rma)), ignore.case=T)] <-1 # Mark "liver" conditions as "1"
> a # Check your work
> design <- cbind(Brain=1, LiverVsBrain=a) # Columnwise bind
> design # Check your work

     Brain LiverVsBrain
[1,]     1            0
[2,]     1            0
[3,]     1            0
[4,]     1            0
[5,]     1            1
[6,]     1            1
[7,]     1            1
[8,]     1            1

Here the first coefficient estimates the mean log-expression for brain tissue and plays the role of an intercept. The second coefficient estimates the difference between brain and liver cells. Differentially expressed genes can be found by

# Identifying differentially expressed genes
> fit <- lmFit(eset.rma, design)
> fit <- eBayes(fit)
> result <- topTable(fit, number=20, sort.by="B", adjust="BH", p.value=0.05) # Get top 20 differentially expressed genes
> result # Check your work

For the second approach, the group-means parametrization, the design matrix can be computed by

# Design matrix: separate group coefficients
> design <- cbind(Brain=c(rep(1,4),rep(0,4)), 
>                 Liver=c(rep(0,4),rep(1,4))) # Manually create design matrix
> design # Check

Or another way, less prone to typo errors

> design <- model.matrix(~0+factor(a)) # factor makes two levels, one for each group
> colnames(design) <- c("Brain", "Liver") # Label columns properly
> design # Check your work

Group-means parametrization should be converted into contrast matrix

> fit <- lmFit(eset.rma, design)
> cont.matrix <- makeContrasts(BrainvsLiver=Brain-Liver, levels=design) # Make matrix of contrasts
> cont.matrix # See what's inside
> fit2 <- contrasts.fit(fit, cont.matrix)
> fit2 <- eBayes(fit2)
> topTable(fit2, adjust="BH")

The above approaches for two groups extend easily to any number of groups. Suppose that we want to pairwise compare all four conditions. An appropriate design matrix can be created using

# Several groups
> a <- c(1, 1, 2, 2, 3, 3, 4, 4) # Four conditions, two replicates per condition
> design <- model.matrix(~0+factor(a)) # Now we have four levels for design matrix
> colnames(design) <- c("B", "fB", "fL", "L") # label columns
> design # Check your work

We create contrast matrix for three pairwise comparisons, for the sake of visualizing the results in a form of Venn diagram (can't take 4 comparisons). Finding differentially expressed genes are the same

> contrast.matrix <- makeContrasts(B-fB, L-fL, B-L, levels=design) # Make three contrasts
> fit <- lmFit(eset.rma, design)
> fit2 <- contrasts.fit(fit, contrast.matrix)
> fit2 <- eBayes(fit2)

Use decideTests function to have a summary of the results for Venn diagram, and visualize it.

We can save our results into a file:

> write.table(topTable(fit2,coef=1,number=1000,adjust.method="BH",p.value=0.05,lfc=log2(2)),"filename.txt",sep="\t") # vary coefficient to write corresponding results to a tab-separated file

