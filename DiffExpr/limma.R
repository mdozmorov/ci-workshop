source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
biocLite("limma")
library(limma)
limmaUsersGuide() # Opens pdf manual for limma
library(affy)
eset.rma <- justRMA(celfile.path="../../Exercises/AffyPreprocessing/Su_CELs/") # RMA summarization of the CEL files
pData(eset.rma) # Check what samples we have

# Design matrix: Treatment-constrast parametrization
a<-rep(0,length(pData(eset.rma)$sample)) # Create a vector of 0
a[grep("liver", rownames(pData(eset.rma)), ignore.case=T)] <-1 # Mark "liver" conditions as "1"
a # Check your work
design <- cbind(Brain=1, LiverVsBrain=a) # Columnwise bind
design # Check your work
# Identifying differentially expressed genes
fit <- lmFit(eset.rma, design)
fit <- eBayes(fit)
result <- topTable(fit, number=20, sort.by="B", adjust="BH", p.value=0.05) # Get top 20 differentially expressed genes
result # Check your work

# Design matrix: separate group coefficients
design <- cbind(Brain=c(rep(1,4),rep(0,4)), 
                Liver=c(rep(0,4),rep(1,4))) # Manually create design matrix
design # Check
design <- model.matrix(~0+factor(a)) # factor makes two levels, one for each group
colnames(design) <- c("Brain", "Liver") # Label columns properly
design # Check your work
fit <- lmFit(eset.rma, design)
cont.matrix <- makeContrasts(BrainvsLiver=Brain-Liver, levels=design) # Make matrix of contrasts
cont.matrix # See what's inside
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")

# Several groups
a <- c(1, 1, 2, 2, 3, 3, 4, 4) # Four conditions, two replicates per condition
design <- model.matrix(~0+factor(a)) # Now we have four levels for design matrix
colnames(design) <- c("B", "fB", "fL", "L") # label columns
design # Check your work
contrast.matrix <- makeContrasts(B-fB, L-fL, B-L, levels=design) # Make three contrasts
fit <- lmFit(eset.rma, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
result <- decideTests(fit2, adjust="BH", p.value=0.05, lfc=log2(2))
vennDiagram(result) # How genes differentially expressed in different conditions
vennDiagram(result,include="up") # Only upregulated
vennDiagram(result,include="down") # Or downregulated
write.table(topTable(fit2,coef=1,number=1000,adjust.method="BH",p.value=0.05,lfc=log2(2)),"filename.txt",sep="\t") # vary coefficient to write corresponding results to a tab-separated file
