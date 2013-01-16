###################################################
# Annotation analysis
# Cory Giles
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

#####################
# Annotation packages
##################### 

source("http://bioconductor.org/biocLite.R")

biocLite("hgu133a.db")

library(hgu133a.db)

# The hgu133a package is one of many microarray annotation packages that provides
# mappings from probes to many other data types.
ls("package:hgu133a.db")

# Annotation packages are transparent to the storage mechanism used, meaning
# they can pull data from local or even remote databases. Usually, they store
# their information in a SQLite database and pull in information only when 
# queried. 

# You can convert a particular annotation mapping to a list with as.list:
probe2entrez <- as.list(hgu133aENTREZID)
probe2entrez[1:5]

# Or the reverse mapping:
entrez2probe <- as.list(revmap(hgu133aENTREZID))
entrez2probe[1:5]

# If you have a set of probes, you can use these mappings:
probes <- c("215134_at", "221987_s_at", "214740_at", "205903_s_at")
entrez <- unlist(probe2entrez[probes])

# The "unlist" method above only works when you have a "many-to-one" or "one-to-one"
# mapping of keys to values. If you have "one-to-many" or "many-to-many", you will 
# get back several mapped values for each key:
entrez2probe[entrez]

########################################################################
# genefilter & collapseRows - 
#   removing genes with certain properties before DE analysis
######################################################################## 

# Filtering any unimportant genes before DE analysis is quite desirable because
# it increases the statistical power of your DE analysis. In other words, a smaller
# gene set suffers less of a hit from multiple testing correction.

# Another potential use for filtering is reducing the feature space for prediction
# analysis. This will cause models to be trained faster and probably to have
# less variance.

# Of course, the tradeoff in both cases is that the genes you filter 
# may be false negatives: that they may be truly differentially expressed or 
# may carry unique predictive power.

library(ALL)
library(genefilter)
data(ALL)

# Looks like log2-transformed, not quantile normalized
summary(exprs(ALL)[,1:3])

# kOverA filters out probes which do not have at least "k" experiments
# in which they are over expression level "A". The idea here is that 
# the variability of low-expressed genes is harder to distinguish from
# the baseline noise level.
ka.filter <- filterfun(kOverA(5, 4.0))
ix.ka <- genefilter(exprs(ALL), ka.filter)
nrow(exprs(ALL))
sum(ix.ka)
# you could now alter the original ALL structure by
# exprs(ALL) <- exprs(ALL)[ix.ka,]

# nsFilter stands for "non-specific filter". Its main purpose is
# to remove probes which do not map to genes or probes which do
# not have GO categories.
filtered <- nsFilter(ALL, require.entrez=T)

filtered$filter.log

filtered$eset # The new ExpressionSet

# Lots of other filters exist. For example, you can 
# filter by Anova or T-tests with "Anova" or "colttests"
# filters. See some more options:
ls("package:genefilter")

# collapseRows from the WGCNA package serves a similar purpose.
# It takes an expression matrix, a vector of probe IDs, and a vector
# of corresponding category IDs (usually Entrez Gene or Ensembl IDs)
# and collapses all the probes mapping to a given gene into one row.

# It has many possible methods for collapsing the rows (in fact, an entire
# paper was written just about this one function:
# * http://www.biomedcentral.com/1471-2105/12/322
# By default, rows are collapsed by max mean.

library(WGCNA)
library(hgu95av2.db)

M <- exprs(ALL)
probes <- rownames(M)
xx <- as.list(hgu95av2SYMBOL)
symbols <- unlist(xx[probes])

clr <- collapseRows(M, symbols, probes)
names(clr)
M.genes <- clr$datETcollapsed

###########################################
# Differential expression review with limma
###########################################

# Before proceeding into GOEA and other enrichment analyses,
# we need a list of probes and differential expression p-values.
# Here, we review limma usage with the ALL dataset. 

# The ALL dataset has many covariates, which you can view with
# head(pData(ALL))
# The ALL dataset contains either B-cell or T-cell leukemias, and 
# various cytogenetic abnormalities. Here, we compare expression in
# B-cell leukemia patients which have undergone BCR/ABL transition 
# with B-cell leukemia patients which have no assigned molecular biology.

# The major novel aspect of this review compared with the previous
# exercises is that here we are looking at a multifactorial design.

library(limma)

pd <- pData(ALL)
# limma disallows "operator" characters in factor names
mol.biol <- sub("[/-]", ".", pd$mol.biol)
groups <- paste(substr(pd$BT,1,1), mol.biol, sep=".")
groups <- factor(groups, levels=unique(groups))

design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)

M <- normalizeQuantiles(exprs(ALL))
fit <- lmFit(M, design)

contrast <- makeContrasts(
    BCR.ABLvsNEG=B.BCR.ABL-B.NEG,
    levels=design)

head(design)
contrast

fit.c <- contrasts.fit(fit, contrast)
fit.c <- eBayes(fit.c)

tt <- topTable(fit.c, number=nrow(exprs(ALL)))

p.values <- tt$adj.P.Val
names(p.values) <- tt$ID

result <- decideTests(fit.c)
sig.probes <- rownames(result)[result!=0]

# (Feel free to also check out vennDiagram on fit.c)

# What were the symbols of these significant probes?
xx <- as.list(hgu95av2SYMBOL)
sig.symbols <- unlist(xx[sig.probes])
# A few genes show up multiple times, and conversely
# a few probes don't map to genes.
sort(table(sig.symbols)) 

# The same sort of process can be used to map significant probes
# to entrez IDs, gene names, pubmed IDs, etc, etc.

######################################################
# geneplotter - plotting genomic coordinates of probes
######################################################

# geneplotter plots the location of probes on chromosomes. Internally,
# it uses the hgu95av2CHR, hgu95av2CHRLOC, and hgu95av2CHRLOCEND (or the
# analogous maps for a different annotation) to plot probes from an annotation.

biocLite("geneplotter")

library(geneplotter)

chrloc <- buildChromLocation("hgu95av2")
cPlot(chrloc, useChroms=c(1:22, "X", "Y"))
cColor(sig.probes, "red", chrloc)

# Neat looking, but not exceedingly informative. X chromosome looks
# possibly depleted. We would need more rigorous methods to know for sure...

###################################
# Gene Ontology Enrichment Analysis
###################################

biocLite("topGO")

library(ALL)
library(topGO)
data(ALL)

sel.fn <- function(p.vals) { p.vals < 0.01 }
affyLib <- paste(annotation(ALL), "db", sep=".")
go.data <- new("topGOdata",
               ontology="BP", allGenes=p.values, geneSel=sel.fn,
               nodeSize=10, # search GO terms with >= 10 genes
               annot=annFUN.db, affyLib=affyLib)

# A variety of algorithms and ranking statistics are available for running the 
# actual enrichment step.
# A list of algorithms and tests can be found by calling
# whichAlgorithms() and whichTests()
result <- runTest(go.data, algorithm="classic", statistic="fisher")
result.01 <- runTest(go.data, algorithm="weight01", statistic="fisher")
result

# lots of possible combinations of algorithms and statistics:
whichAlgorithms()
whichTests()

# View the top results from multiple methods in tabular format
GenTable(go.data, result, result.01)

# Hmm, platelets. A little literature review shows that ABL is activated by PDGF,
# and that BCR-ABL inhibitor ponatinib causes platelet dysfunction. Lots of interesting
# questions arise here. Are we looking at cancer biology or medication effects or both?

# Show results in the context of the GO DAG
# (Very small text. This can be easier to read if you save into an external file
# with png() and dev.off())
# You may need to have graphviz headers and binaries
# installed on your system for this to work.
biocLite("Rgraphviz")

library("Rgraphviz")

showSigOfNodes(go.data, score(result), firstSigNodes=3, useInfo="all")

#######################################
# Creating a custom enrichment analysis
#######################################

# At its core, enrichment analysis is very simple. You have a list of genes,
# a list of categories, and a list of gene-category associations.

# The question is, given a list of genes (for example, the top differentially
# expressed genes from an experiment), what categories are these genes
# statistically enriched for? The standard test for this is Fisher Exact.

# This simplicity means enrichment analysis is very broadly applicable. Any time
# you have a data source that annotates genes with categories, you can use
# enrichment analysis.

# In the following, we use REACTOME DB, which is a pathway database like KEGG,
# as our source of gene annotations.

# Fisher Exact test with contingency table
ct <- matrix(c(35,165,150,12000),byrow=T,nrow=2)
ct
fisher.test(ct, alternative="greater")

library(topGO)
library(hgu95av2.db)
library(reactome.db)

probe2entrez <- as.list(hgu95av2ENTREZID)
sel.probes <- names(p.values[p.values < 0.05])
all.probes <- names(p.values)
sel.entrez <- unlist(probe2entrez[sel.probes])
all.entrez <- unlist(probe2entrez[all.probes])

entrez2path <- as.list(reactomeEXTID2PATHID)
sel.counts <- table(unlist(entrez2path[sel.entrez]))
all.counts <- table(unlist(entrez2path[all.entrez]))
n.sel <- length(sel.entrez)
n.all <- length(all.entrez)
pathways <- names(sel.counts)
pathway.names <- as.list(reactomePATHID2NAME)

result <- lapply(pathways, function(p) {
                 sel <- sel.counts[[p]]
                 all <- all.counts[[p]]
                 pathway.name <- pathway.names[[p]]
                 m <- matrix(c(sel, all, n.sel-sel, n.all-all),
                             byrow=T,nrow=2)
                 p.value <- fisher.test(m, alternative="greater")$p.value
                 data.frame(PathwayID=p, 
                            Pathway=pathway.name,
                            N.Selected=sel,
                            N.In.Category=all,
                            p.value=p.value)
})

result <- as.data.frame(do.call(rbind, result))
result <- result[order(result$p.value),]

result$p.value <- p.adjust(result$p.value,
                           method="holm")

sum(result$p.value < 0.05)
# Note a relatively small number of significantly enriched terms.
# In enrichment analysis, this is usually less of a problem, because you
# are usually interested in the relative ranks of the terms rather than
# which are significant per se. Although obviously, if an enrichment is
# not significant, you should take it with an extra grain of salt.

head(result)
# Encouragingly, the results give a similar message to our GOEA.
# Independent corroboration is important!

############
# Exercise 1
############

# Using the ALL dataset, do the following:

# 1. Choose a contrast that interests you. 
# For more information about all the pData fields, you can see:
#    * http://www.bioconductor.org/packages/release/data/experiment/manuals/ALL/man/ALL.pdf

# 2. Use genefilter to select only probes which map to genes.

# 3. Encode your contrast of interest into a limma design and contrast matrix.
#    Use limma alongside the filtered expression matrix 
#    to find the differentially expressed probes. 

# 4. Perform Gene Ontology enrichment analysis on the differentially expressed probes.
#    Try using the molecular function ("MF") or cellular component ("CC") ontologies.
#    Try performing GOEA using only upregulated or only downregulated genes. (Recall
#    that this information can be found using the "decideTests" function.)

# 5. Go back to the gene filtering stage, and use genefilter to filter out
#    probes with kOverA(10, 5.0). Re-run the rest of the pipeline. How do the results
#    change?

############
# Exercise 2
############

# Write a function that takes a vector of significantly 
# differentially expressed probes and an annotation name,
# and uses Fisher's Exact test to test whether any chromosomes
# are significantly enriched or depleted for differentially expressed
# genes. Be sure to correct for multiple testing.

# To do this, you might find the "get" function useful. For example:

get(paste("hgu95av2", "CHRLOC", sep=""))

# returns the CHRLOC map for the hgu95av2 annotation package, assuming
# it has been loaded already by library().
