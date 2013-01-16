###################################################
# Graphics
# Mikhail Dozmorov
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

#####################
# Heatmaps
#####################

Using R to draw a Heatmap from Microarray Data

We will be using an example outlined in this paper, "Bioconductor: open software development for computational biology and bioinformatics", Gentleman et al. 2004. [PMID: 15461798]. We will reproduce their Figure 2 (http://genomebiology.com/2004/5/10/r80/figure/F2).

The original citation for the raw data we will be using is "Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival" by Chiaretti et al. Blood 2004. [PMID: 14684422]. The example dataset can be downloaded as the Bioconductor ALL package.

> source("http://www.bioconductor.org/biocLite.R")
> biocLite("ALL") # Install ALL package
> library("ALL") # Load ALL package
> data("ALL") # Attach ALL data to the workspace
> ALL # Check what's in there

We got expression matrix with  128 samples and 12625 genes, and associated phenotypic data. We can looks at the parameters available to us in the phenotypic data, and the results of molecular biology testing for the 128 samples.

> names(pData(ALL)) # Check what phenotype data we have
> ALL$mol.biol # Results of molecular biology testing for the 128 samples

Ignoring the samples which came back negative on this test (NEG), most have been classified as having a translocation between chromosomes 9 and 22 (BCR/ABL), or a translocation between chromosomes 4 and 11 (ALL1/AF4).

For the purposes of this example, we are only interested in these two subgroups, so we will create a filtered version of the dataset using this as a selection criteria:

> eset <- ALL[, ALL$mol.biol %in% c("BCR/ABL", "ALL1/AF4")] # Select a subset of the data

The resulting variable, eset, contains just 47 samples - each with the full 12,625 gene expression levels.

This is far too much data to draw a heatmap with, but we can do one for the first 100 genes using built-in R function:

> heatmap(exprs(eset[1:100,])) # An example of a heatmap

According to the BioConductor paper we are following, the next step in the analysis was to use the lmFit function (from the limma package) to look for genes differentially expressed between the two groups.

# Finding differentially expressed genes
> library("limma") # Load limma package

The fitted model object is further processed by the eBayes function to produce empirical Bayes test statistics for each gene, including moderated t-statistics, p-values and log-odds of differential expression.

> f <- factor(as.character(eset$mol.biol)) # Convert group data into factors
> design <- model.matrix(~f) # Create design matrix as a function of factors
> fit <- eBayes(lmFit(eset,design)) # Fit a model and obtain empirical Bayes test statistics

Let's look at the top 10 differentially expressed genes. The leftmost numbers are row indices, ID is the Affymetrix HGU95av2 accession number, logFC is the log ratio of expression, AverExpr is the log average expression, t the moderated t-statistic, and B is the log odds of differential expression (the more the better). We specify coef = 2 to see the genes from the second column of the design matrix.

> topTable(fit, coef = 2) # Look at the top 10 differentially expressed genes

Next, we select those genes that have adjusted p-values below 0.05, using a very stringent Holm method to select a small number (165) of genes.

> selected  <- p.adjust(fit$p.value[, 2], method = "holm") <0.05 # Indexes based on filtering adjusted p-values
> esetSel <- eset [selected, ] # Select a subset of genes

The variable esetSel has data on (only) 165 genes for all 47 samples (check). Now, again using generic R function, we can see the heatmap of truly differentially expressed genes, and play with colors.

> heatmap(exprs(esetSel)) # Heatmap of differentially expressed genes

We're almost there (check http://genomebiology.com/2004/5/10/r80/figure/F2), except they have added a red/blue banner across the top to really emphasize how the hierarchical clustering has correctly split the data into the two groups (10 and 37 patients).
To do that, we can use the heatmap function's optional argument of ColSideColors. I created a small function to map the eselSet$mol.biol values to red (#FF0000) and blue (#0000FF), which we can apply to each of the molecular biology results to get a matching list of colors for our columns:

# color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" } # Function to map group data to colord
patientcolors <- unlist(lapply(esetSel$mol.bio, color.map)) # Apply this function to group data
# heatmap(exprs(esetSel), col=topo.colors(100), ColSideColors=patientcolors) # Use these colors as column side colors

One subtle point in the previous examples is that the heatmap function has automatically scaled the colours for each row (i.e. each gene has been individually normalised across patients). This can be disabled using scale="none", which you might want to do if you have already done your own normalisation (or this may not be appropriate for your data):

> heatmap(exprs(esetSel), col=topo.colors(75), scale="none", ColSideColors=patientcolors, cexRow=0.5) # Disable color rescaling

Note cexRow parameter, used to decrease font size of row labels. cexCol does the same for column labels.

Generic R function does not provide us with an opportunity to see a key, or a legend to color coding. The good news is that heatmap.2 from the gplots library does. In fact, it offers a lot of other features, disabled here, check them with help(heatmap.2)

> heatmap.2(exprs(esetSel), col=topo.colors(75), scale="none", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5) # Heatmap with color key

By default, heatmap.2 will also show a trace on each data point (removed this with trace="none"). If you ask for a key (using key=TRUE) this function will actually give you a combined "color key and histogram", but that can be overridden (with density.info="none").

We can use different color schemes. Try using the functions bluered/redblue for a red-white-blue spread, or redgreen/greenred for the red-black-green colour scheme often used with two-colour microarrays:

> heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5) # Heatmap with red/green color scheme

If we're not satisfied with the default distance and hierarchical clustering metrics (euclidean/complete), we can try other options by trying all combinations, outputting the results into a pdf file and then select most visually pleasant/biologically relevang combination of metrics.

# Testing different metrics
> dist.methods<-c("euclidean",  "manhattan", "binary", "minkowski","canberra","maximum") # Define distance metrics
> hclust.methods<-c("ward", "single", "complete", "average", "mcquitty", "median", "centroid") # Clustering metrics
> dev.off() # Clear graphic window
> pdf("test.pdf") # Open PDF file as graphic device
> for (d in dist.methods) {
>   for (h in hclust.methods){
>     heatmap.2(exprs(esetSel),col=redgreen(75),distfun=function(x){dist(x,method=d)}, hclustfun=function(x){hclust(x,method=h)}, , scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, main=paste("Dist : ",d,"; Hclust : ",h)) # Heatmap with different distance/clustering metrics
>   }
> }
> dev.off()

After inspecting visual appearance of different heatmaps, delete the file.

> unlink("test.pdf") # Delete file

It is often important to see which samples belong to which subclusters. We can use cutree function to cut the dendrogram tree either by height or by the number of groups. Note one has to transpose the matrix, as dist function works on rows.

# Extracting subclustered groups
> colDendrogram<-(hclust(dist(t(exprs(esetSel)),method="euclidean"),method="complete")) # Make column dendrogram
> plot(tmp1) # visualized
> cutree(colDendrogram,h=25) # Sample-group assignment, h - height
> cutree(colDendrogram, k=2) # Same effect, k - number of groups
> sort(cutree(colDendrogram, k=2)) # Sort for easier group separation