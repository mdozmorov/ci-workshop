###################################################
# Network analysis
# Cory Giles
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

#############################
# Manipulating generic graphs
#############################

library(igraph)

g <- graph(1:10)
plot(g)

g <- graph(1:10, directed=F)
plot(g)

# In igraph, the "barabasi.game" function models a "scale-free" network.
# A scale-free network is a network whose connectivity follows a power-law distribution.
# In other words, there are a few "hubs" with high connectivity and many more
# non-hub nodes. Biological networks and the internet are among the many real-world
# examples of networks that follow this pattern.

# The Barabasi model builds scale-free networks by "preferential attachment". In other words,
# when a new node is added, it is more likely to be added to a node which has more
# pre-existing connections. In PA networks, "the rich get richer".

# (For more information about scale-free networks, see:  
#   Barabási, Albert-László; Albert, Réka. (October 15, 1999). 
#   "Emergence of scaling in random networks". Science 286 (5439): 509-512.)

# Other generators of random networks with interesting properties in igraph:
# erdos.renyi.game, degree.sequence.game, callaway.traits.game, establishment.game

g <- barabasi.game(50, directed=F)
plot(g)

# Nodes have IDs, names (which can be used to refer to them in 
# function calls), and labels (which are used for plotting).
# By default the node names and labels are equal to their IDs.
# They can also be assigned arbitrary attributes using this assignment syntax.

V(g)

V(g)[1]

V(g)$name <- c(letters,LETTERS)[1:50]

V(g)

V(g)$label <- c(letters,LETTERS)[1:50]

plot(g)

E(g)

# Different metrics of network connectivity and topology:

# Returns a vector representing the percentage of nodes having each degree
# ("degree" = number of connections)
degree.distribution(g)

# Returns the degree of each node
degree(g)

average.path.length(g)

# Use Google PageRank to estimate the importance of each vertex:
pr <- page.rank(g, V(g))

# Check if the PageRank matches up with your intuitions
sort(pr$vector, decreasing=T)[1:5]
plot(g)

# Adjacency matrix manipulation:
# Recall that self-multiplying an adjacency matrix N times
# returns a matrix describing the number of N+1 length paths
# between each pair of nodes.
A <- get.adjacency(g)

A[1:5,1:5]

(A %*% A)[1:5,1:5]

(A %*% A %*% A)[1:5,1:5]

(A %*% A %*% A %*% A)[1:5,1:5]

# (Q: why does the diagonal disappear on odd-numbered exponentiations?)

# Shortest paths, using Dijkstra's algorithm:
# Get lengths of all shortest paths or specific shortest paths
shortest.paths(g)

shortest.paths(g,1:5,1:5)

# Get the actual node sequence for a shortest path
get.shortest.paths(g, "a", "Q")

ls("package:igraph")

##################################################
# Construct and view a simple coexpression network
##################################################

library(ALL)
library(hgu95av2.db)
library(graph)
library(Rgraphviz)
library(preprocessCore)
library(WGCNA) # Provides both collapseRows and other WGCNA functions
library(multicore)

data(ALL)

# Quantile normalize experiments 
M <- normalize.quantiles(exprs(ALL))
dimnames(M) <- dimnames(exprs(ALL))

probes <- rownames(M)
xx <- as.list(hgu95av2SYMBOL)
symbols <- unlist(xx[probes])
M.genes <- collapseRows(M, symbols, probes)$datETcollapsed

# Construct correlation graph
corGraph <- function(X,k=50) { 
    edges <- mclapply(1:ncol(X), function(i) {
             # coerce from matrix to vector
             cors <- as.numeric(cor(X[,i],X)) 
             # Start at 2 b/c element 1 will be the node itself
             order(-cors)[2:k+1]
        }, mc.preschedule=T)
    nodes <- colnames(X)
    # graph.adjlist (adjacency list) is an igraph constructor 
    # that takes a list of source nodes,
    # each containing a vector of target nodes.
    g <- graph.adjlist(edges)
    V(g)$name <- nodes
    V(g)$label <- nodes
    class(g) <- c("igraph", "corGraph")
    g
}

# Just plot coexpression network for first few genes
g <- corGraph(t(M.genes)[,1:75],k=5)

# Not very scale-free! This is due to our approach in generating edges...
plot(g)

write.graph(g, "cox.dot", format="dot")

# Note that other simple approaches for constructing a coexpression network
# are by generating all-vs-all Pearson correlations and thresholding at
# a) Some statistical significance, or
# b) Some simple Pearson threshold,

###############################################
# Simple coexpression-based function prediction
###############################################

# Now, let's constuct a coexpression network for all genes,
# returning the top 20 most highly correlated genes for each gene.
# This time using Entrez ID instead of symbol.

probes <- rownames(M)
xx <- as.list(hgu95av2ENTREZID)
entrezIds <- unlist(xx[probes])
M.genes <- collapseRows(M, entrezIds, probes)$datETcollapsed
g <- corGraph(t(M.genes), k=20)

# This may take a few minutes, depending on your computer's speed.

# In the following section, we are going to synthesize
# what we have learned about annotation, prediction, 
# and network analysis to construct a simple application 
# to predict the function of unknown genes. The purpose
# of this exercise is to learn about AFP, and also to get a feel
# for the development process of a small real-world R app.

# The idea is simple: two genes which have high pairwise
# correlation coefficients are likely to share function.
# In other words, for any gene, we can look at its coexpression
# neighbors to infer its function.

# Guilt-by-association based functional prediction is being
# actively used to prioritize unknown transcripts. The most
# popular methods are more sophisticated and use many different
# kinds of data, such as PPI networks and homology. 
# However, they follow the same general principles.

# Our function will have 3 steps:
# 1. Find the coexpressed genes for the query gene 
#    (look them up in the corGraph).
# 2. Run GOEA on the coexpressed genes.
# 3. Return the significant functions in a readable data frame.

# To achieve this, we need to write a function that does GOEA,
# and a function that wraps up the graph search and GOEA.
# (Note that you could use topGO instead of writing your
# own function.)

library(org.Hs.eg.db)

go2eg <- as.list(org.Hs.egGO2EG)

# Unfortunately, "revmap" does not work correctly for 
# org.Hs.egGO2EG. So we have to travel deep into the structure
# of org.Hs.egGO to get the information of interest.

as.list(org.Hs.egGO[1:5])

# Notice the use of the "block" structure in defining eg2go,
# which can be useful to avoid polluting the global namespace.

eg2go <- {
    xx <- as.list(org.Hs.egGO)
    result <- lapply(xx, function(geneMappings) {
        if (!is.na(geneMappings)[1])
            sapply(geneMappings, function(m) {m$GOID})
    })
    names(result) <- names(xx)
    result
}

goea.background <- function(genes) {
    bg <- table(unlist(eg2go[genes]))
    attr(bg, "ngenes") <- length(genes)
    bg
}

goea <- function(genes, background) {
    # Given a list of human entrez IDs, return a vector of 
    # GO IDs and p-values.
    # Note that even though eg2go countains NULL values,
    # "table" removes them.
    sel.counts <- table(unlist(eg2go[genes]))
    n.sel <- length(genes)
    n.all <- attr(background, "ngenes")
    pvals <- sapply(names(sel.counts), function(term) {
        print(term)
        sel <- sel.counts[[term]]
        all <- background[[term]]
        if (sel > n.sel) {
            #print(genes)
            print(sel.counts)
            stop()
        }
        ct <- matrix(c(sel, all, n.sel-sel, n.all-all),
                     byrow=T, nrow=2)
        fisher.test(ct, alternative="greater")$p.value
    })
    pvals <- p.adjust(pvals, method="holm")
}

# This will extend the "predict" function to operate 
# on "corGraph" types. This is called "extending an S3 method", 
# and we will talk more about it next session.
predict.corGraph <- function(gr, genes, pval.cutoff=0.05) {
    background <- goea.background(V(gr)$name)
    # "do.call(rbind, ...) is a common way to merge data frames
    result <- do.call(rbind, 
        lapply(genes, function(gene) {
               coexNodeIds <- neighbors(gr, gene, mode="out")
               coexGenes <- V(gr)$name[coexNodeIds]
               p.values <- goea(coexGenes, background)
               p.values <- p.values[p.values<pval.cutoff]
               data.frame(GeneID=rep(gene, length(p.values)),
                          GOID=names(p.values),
                          PValue=p.values)
        })
    )
    # Add textual GO and gene names
    go2term <- lapply(as.list(GOTERM), Term)
    result$GOTerm <- factor(unlist(go2term[result$GOID]))
    genename <- as.list(org.Hs.egGENENAME)
    result$GeneName <- factor(unlist(genename[result$GeneID]))

    # Reset rownames to autoincremented integer
    rownames(result) <- NULL 
    result[order(result$PValue),
           c("GeneID", "GeneName", "GOID", "GOTerm", "PValue")]
}

allGenes <- V(g)$name
result <- predict(g, allGenes)

head(result)

# These functions could be improved by allowing the user
# to specify the organism, probeset, etc.

# Most importantly, to know how much we should trust these results,
# we would have to run performance metrics like precision, recall,
# F-measure, AUC, lift, etc on our results. Probably lift is the 
# best criterion considering the purpose of function prediction. 
# AUC is also very popular in the AFP field. Recall is not very 
# important in this instance. Unfortunately, we lack time to do 
# this, but the necessary gold standards are readily available (the 
# true annotations of genes with GO terms).


# For more information about automated function prediction using coexpression data, see:

# * The "MouseFunc" competition: http://hugheslab.med.utoronto.ca/supplementary-data/mouseFunc_I
# * Guan Y, Myers CL, Hess DC, Barutcuoglu Z, Caudy AA, Troyanskaya OG. (2008) Predicting gene function in a hierarchical context with an ensemble of classifiers. Genome Biol. Suppl 1:S3.
# * Dozmorov MG, Giles CB, Wren JD, 2011. "Predicting Gene Ontology from a global meta-analysis of 1-color microarray experiments." BMC Bioinformatics.

# Recently, a series of papers has been published which shows
# limitations to the method. See the following references:

# * Gillis J, Pavlidis P, 2012. "'Guilt by association' is the exception rather than the rule in gene networks." PLoS Comput Biol.
# * Gillis J, Pavlidis P, 2011. "The role of indirect connections in gene networks in predicting function." Bioinformatics.

##################################################################
# WGCNA
##################################################################
# This section adapted from Horvath S, Weighted Network Analysis. 
# Applications in Genomics and Systems Biology.
# http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Book/RcodeChapter12NEW.pdf

# The purpose of WGCNA is to find network modules from 
# coexpression data. The expression of each module can then be
# summarized as an "eigengene", which can simplify understanding
# of the network and things like annotation analysis.

# In the following, we choose the beta coefficient which maximizes 
# the scale-free criterion of the resulting coexpression network.
# The resulting adjacency matrix will be 
# = A + (A * A) + (A * A * A) + ... + A^beta
# where A is the original unweighted or weighted adjacency matrix 
# In the jargon, this is called soft thresholding using the 
# power adjacency function. A "simple" threshold as described
# above is called a step adjacency function.

# For more explanation, see:
# * http://labs.genetics.ucla.edu/horvath/GeneralFramework/NetworkConstruction.pdf
# * http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/OverviewWGCNA.pdf
# * Langfelder and Horvath 2008, "WGCNA: an R package for weighted
#   correlation network analysis." BMC Bioinformatics.

# To partially explain the motivation for a power adjacency function,
# recall that for an unweighted (binary) adjacency matrix A, 
# each slot in the matrix A^n shows the number of paths of length n
# between two nodes. However, in WGCNA, the original adjacency matrix A
# is not binary, it is filled with the absolute value of pairwise Pearson
# correlation coefficients.

###################################################################

allowWGCNAThreads()

# All WGCNA functions require genes in columns
M.genes <- t(M.genes)

powers <- c(1:10, seq(12,30,2))
sft <- pickSoftThreshold(M.genes, powerVector=powers)
sft

# The optimal beta is a tradeoff: 
# Beta too high -> false negative edges.
# Beta too low -> false positive edges.
# The decision is made primarily by looking at "SFT.R.sq", which quantifies
# how "scale-free" a network is. It is thought that the more scale-free the coexpression
# network, the more biologically plausible it is. 

# One common heuristic is to choose the smallest value of beta
# which results in a SFT R^2 > 0.85.
# Before making this decision we also multiply by the negative signum of the "slope".
# Read the docs for more details, but essentially if a slope is positive, the network
# has more hub genes than non-hub genes, which is considered implausible. 

# Note the similarity to our old friend, the bias-variance tradeoff.
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale-free topology, signed R^2", 
     main=paste("Scale Independence"))

# Note how increased beta -> lower connectivity. This occurs because
# beta is essentially the number of self-multiplications of the adjacency
# matrix. Powers of numbers < 1 rapidly approach 0.
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", main="Mean Connectivity")

sft.criterion <- sft$fitIndices$SFT.R.sq * -sign(sft$fitIndices$slope)
beta <- sft$fitIndices$Power[sft.criterion>=0.85][1]

# Instead of picking a beta to multiply A by, you can instead search 
# for a "hard" Pearson cutoff to apply to all pairs, 
# also using this scale-free criterion.
hft <- pickHardThreshold(M.genes, T, RsquaredCut=0.85)
hft.criterion <- hft$fitIndices$SFT.R.sq * -sign(hft$fitIndices$slope)
pearson.cut <- hft$fitIndices$Cut[hft.criterion>=0.85][1]
pearson.cut

# The adjacency matrix (filled with double precision floats) 
# is of size N^2 and uses 8 * N^2 bytes of RAM, 
# where N is the number of probes.
# In this case, with N=8805, this is only ~600 MB, which is easily
# attainable by most PCs. But larger matrices take considerably more space
# -- the Affymetrix U133a Plus 2.0, which has 22000 probes, takes 4GB.
# A RNAseq-derived adjacency matrix for all UCSC knownGenes (~70000) 
# would take ~40GB. These sizes are not feasible on most PCs. 

# So WGCNA has the ability to perform coexpression clustering "blockwise", 
# without loading the entire adjacency matrix into RAM. 
# It works in two phases:
# 1) A computationally inexpensive crude clustering into "blocks" of an appropriate
#    size for your PC.
# 2) Sequential fine-grained hierarchical clustering of each block.

# Another advantage of the function "blockwiseModules", is that the 
# user does not have to directly interact with the adjacency matrix 
# and can let WGCNA find the modules itself.

# However, there are several keyword arguments to this function,
# such as saveTOMs,which allow you to save, load, or set parameters for the 
# Topological Overlap Matrix (which is the adjacency matrix after 
# several self-multiplications as described above).
# See ?blockwiseModules for a list of options.

mods <- blockwiseModules(M.genes, corType="pearson",
                         maxBlockSize=5000, networkType="unsigned",
                         power=beta, minModuleSize=30)
names(mods) 

# WGCNA gives each cluster a "color".
# A module eigengene (ME) is a way of summarizing the expression level of each
# "color" in each sample. An eigengene is the first principal component of the submatrix
# containing the genes in each cluster. In short, it is a way of summarizing the expression
# level of a cluster within a sample.
head(mods$MEs)

# mods$colors shows which cluster each gene is assigned to
# Note that clusters can diverge (widely) in size.
# The "grey" color is used for genes not belonging to any module.
head(mods$colors)
table(mods$colors)

# If there is time, we will also discuss another important feature of WGCNA,
# correlating modules with traits.

######################
# Additional resources
######################

# * The KEGGgraph package has KEGG in graph format. Unfortunately,
#   KEGG is no longer free, so the data is old.

# * RBGL (R Boost Graph Library), an alternative to igraph, is a R wrapper to the popular C++ Boost Graph Library.

# * Although igraph is good for manipulating graphs, it is not great for visualization. Cytoscape is a very full-featured Java program that can visualize graphs you have created with igraph and exported to an appropriate format.

###########
# Exercises
###########

# 1. Gene Ontology is a directed acyclic graph whose edges can be
#    found in GOBPCHILDREN, GOCCCHILDREN, and GOMFCHILDREN.
#    (You must first run "library(GO.db)" to access 
#    these annotation objects.)

#    Load GO into an igraph graph. Visualize it. Examine
#    its network topology properties. 
#    What nodes are considered most important by PageRank?

# 2. Modify the AFP program in one of the following ways:
#    a. Modify "cor.graph" to rank coexpression partners 
#       by the absolute value of their correlation coefficient 
#       instead of considering those with positive correlation 
#       coefficients only.
#    b. Modify "cor.graph" to use a hard correlation threshold 
#       (i.e., to consider all coexpression partners with a 
#       correlation coefficient above a user-specified threshold, 
#       perhaps defaulting to 0.7). 

# 3. Develop a performance evaluation (precision, recall, AUC, and/or lift) 
#    for the AFP program. If you use precision, recall, or lift, remember
#    that you need some kind of negative control to determine the baseline
#    performance. The easiest way to do this is to randomize the coexpression
#    graph. If you use AFP, you would have to modify the code so that p-values are returned for ALL GO-Gene pairs.
