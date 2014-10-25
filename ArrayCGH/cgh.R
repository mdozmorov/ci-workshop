###################################################
# Array CGH
# Cory Giles
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

library(CGHcall)

# Cervical cancer data
data(Wilting)

# Crucial information for an aCGH experiment is:
# - probe ID
# - genomic location of probe (chrom, start, end)
# - probe logratio between case and control
head(Wilting)

# Transform data.frame into Biobase AnnotatedDataFrame.
# This function can take a data.frame or an appropriately
# formatted TSV file.
cgh.raw <- make_cghRaw(Wilting)

# Note the data structure
cgh.raw
attributes(cgh.raw)

# The vaguely-named "preprocess" function:
# - Filters data w/ missing position information
# - Can impute missing values with kNN algorithm
# - Removes probes with too many missing values
cgh.pre <- preprocess(cgh.raw)

# Note that "preprocess" is a very generic name.
# For names like these, it can be make clearer code to use the 
# fully qualified package namespace. 
# So, the following line is equivalent:
cgh.pre <- CGHcall::preprocess(cgh.raw)

# This function (obviously) normalizes the CGH data, using
# either median or mode.
# In this case, the normalization has a straightforward
# interpretation: it reflects the assumptions that most probes
# will have a copy number of 2 in both case and control
# (and thus the mean copy number logratio should be very 
# close to zero)
cgh.norm <- CGHcall::normalize(cgh.pre, 
                               method="median", smoothOutliers=T)

# This wraps the DNAcopy algorithm and performs segmentation.
# Note the "clen" and "relSDlong" parameters, which can be used
# to penalize short segments from being called. This is useful
# because excessively short segmentation can cause plots to 
# look noisy.
cgh.seg <- segmentData(cgh.norm, method="DNAcopy")

# Check the ?postsegnormalize for details.
# Essentially, the average segment (rather than the average
# probe as above) has its logratio normalized to zero.
cgh.postseg <- postsegnormalize(cgh.seg)

plot(cgh.postseg[,1])

# This function runs the expectation-maximization (EM) algorithm to estimate
# the probability of a gain or loss of each segment.

# This function takes a parameter called "cellularity"
# which is the fraction of each sample that is comprised of tumor cells.
# CGHcall can use this information, if you have it, to correct for
# different fractions of tumor cells in each sample. 

# This function also takes "organism" and "build" parameters,
# so that CGHcall knows what the genome "background" is.
# It defaults to "human" and "GRCh37", which is ok for this 
# dataset.
result <- CGHcall(cgh.postseg, nclass=5)

# This turns the raw list data structure output by "CGHcall"
# into a "cghCall" object, which can be used for plotting, etc.
# Part of the cellularity normalization is also performed
# in this function.
# (This is probably an example of bad API design.)
result.exp <- ExpandCGHcall(result, cgh.postseg)

# Plot individual samples.
plot(result.exp[,1])

plot(result.exp[,2])

# Plot mean probabilities of gain or loss across samples.
summaryPlot(result.exp)

# Summarize the frequency of called gains or losses 
# across all samples.
frequencyPlotCalls(result.exp)

# Dig into the actual calls.
# The calls are coded with numbers:
# -2 = double deletion
# -1 = single deletion
# 0 = no change
# 1 = copy number gain
# 2 = amplification
head(calls(result.exp))
table(calls(result.exp))

head(copynumber(result.exp))

# These posterior probabilities can be useful for 
# downstream analyses:

head(probloss(result.exp))

head(probgain(result.exp))

head(probnorm(result.exp))

# Surprisingly, AdCA does not cluster by itself:
heatmap(cor(calls(result.exp)))

# The next step in many analyses is to take called loci
# and find whether there are significant differences between
# groups of samples. This can be done using significance
# testing and multiple hypothesis correction techniques
# we've discussed in earlier sessions.

# In this case, of course, the sample size is insufficient to
# search for differences between AdCA and SCC (it is only
# a subset of the original data published in Wilting et al, 2006).

######################
# Additional resources
######################

# This has been a very basic introduction to aCGH analysis.
# For more information, check out the DNAcopy and CGHcall
# papers, which were published in Bioinformatics:

# * Venkatraman & Olshen, 2007. "A faster circular binary segmentation algorithm for the analysis of array CGH data." Bioinformatics.
# * van de Wiel et al, 2007. "CGHcall: calling aberrations for array CGH tumor profiles."

# Also, there are several broad overviews 
# of the aCGH technology in general:

# * Kallionemi, 2008. "CGH microarrays and cancer." Curr Opin Biotechnol. 
# * http://www.cs.cmu.edu/~epxing/Class/10810-05/Lecture11.pdf
# * http://www.ncbi.nlm.nih.gov/pubmed/19696946
