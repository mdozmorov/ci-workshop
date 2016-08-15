# source("http://www.bioconductor.org/biocLite.R")
# biocLite("ALL") # Install ALL package
library("ALL") # Load ALL package
data("ALL") # Attach ALL data to the workspace
ALL # Check what's in there
names(pData(ALL)) # Check what phenotype data we have
ALL$mol.biol # Results of molecular biology testing for the 128 samples
eset <- ALL[, ALL$mol.biol %in% c("BCR/ABL", "ALL1/AF4")] # Select a subset of the data
heatmap(exprs(eset[1:100,])) # An example of a heatmap
# Finding differentially expressed genes
library("limma") # Load limma package
f <- factor(as.character(eset$mol.biol)) # Convert group data into factors
design <- model.matrix(~f) # Create design matrix as a function of factors
fit <- eBayes(lmFit(eset,design)) # Fit a model and obtain empirical Bayes test statistics
topTable(fit, coef = 2) # Look at the top 10 differentially expressed genes
selected  <- p.adjust(fit$p.value[, 2], method = "holm") <0.05 # Indexes based on filtering adjusted p-values
esetSel <- eset [selected, ] # Select a subset of genes
heatmap(exprs(esetSel)) # Heatmap of differentially expressed genes
heatmap(exprs(esetSel), col=topo.colors(100)) # Heatmap with topographical colors
color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" } # Function to map group data to colord
patientcolors <- unlist(lapply(esetSel$mol.bio, color.map)) # Apply this function to group data
heatmap(exprs(esetSel), col=topo.colors(100), ColSideColors=patientcolors) # Use these colors as column side colors
heatmap(exprs(esetSel), col=topo.colors(75), scale="none", ColSideColors=patientcolors, cexRow=0.5) # Disable color rescaling
library("gplots") # Load gplots package
heatmap.2(exprs(esetSel), col=topo.colors(75), scale="none", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5) # Heatmap with color key
heatmap.2(exprs(esetSel), col=redgreen(75), scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5) # Heatmap with red/green color scheme
# Testing different metrics
dist.methods<-c("euclidean",  "manhattan", "binary", "minkowski","canberra","maximum") # Define distance metrics
hclust.methods<-c("ward", "single", "complete", "average", "mcquitty", "median", "centroid") # Clustering metrics
dev.off() # Clear graphic window
pdf("test.pdf") # Open PDF file as graphic device
for (d in dist.methods) {
  for (h in hclust.methods){
    heatmap.2(exprs(esetSel),col=redgreen(75),distfun=function(x){dist(x,method=d)}, hclustfun=function(x){hclust(x,method=h)}, , scale="row", ColSideColors=patientcolors, key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5, main=paste("Dist : ",d,"; Hclust : ",h)) # Heatmap with different distance/clustering metrics
  }
}
dev.off()
unlink("test.pdf") # Delete file
# Extracting subclustered groups
colDendrogram<-(hclust(dist(t(exprs(esetSel)),method="euclidean"),method="complete")) # Make column dendrogram
plot(colDendrogram) # visualized
cutree(colDendrogram,h=25) # Sample-group assignment, h - height
cutree(colDendrogram, k=2) # Same effect, k - number of groups
sort(cutree(colDendrogram, k=2)) # Sort for easier group separation
