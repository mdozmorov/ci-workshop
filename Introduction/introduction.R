###################################################
# Introductory tutorial to R/Bioconductor
# Cory Giles
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

##############################
# Introduction to R data types
##############################

# Vectors
# "<-" is the assignment operator. "=" can also be used, although it has
# some subtle differences which usually don't matter.
numbers <- 3:7
numbers[1:3]
# c(...) is alternative syntax for vector creation
numbers[c(1,3,5)]

length(numbers)
sum(numbers)
prod(numbers)
min(numbers)
max(numbers)
cumsum(numbers)
cumprod(numbers)
log(numbers)

# Vectors must have homogeneous type
# character, integer, numeric (floating point), logical
myVector <- c(1, "foo", 2.2, 0)
as.character(myVector)
as.integer(myVector)
as.numeric(myVector)

v <- c(1,0,2)
as.logical(v)

# Vectorized arithmetic
# Note that R does not have a scalar type. A "scalar" is simply a vector
# of length 1.

v <- 1:10
v * 3
v + 2
v - 5
v / 2
v ^ 2
v %% 2

# Logical vectors can be used for indexing, ...
ix <- (v %% 2) != 0
v[ix]

# But integer vectors retrieve elements by index:
ix <- (v %% 2)
v[ix]

# Matrices
m <- matrix(1:9, nrow=3, byrow=3)
m
m[2:3,1:2]
rownames(m) <- c("a","b","c")
rownames(m)
colnames(m) <- c("x","y","z")
m
nrow(m)
ncol(m)
dim(m)
t(m)

?matrix

# The class of any object can be inspected with "class"
class(m)
class(v)

# Lists are containers for heterogeneous objects

lst <- list(1,c("foo","bar"),c(3.5,5))
lst
lst[3]
lst[[3]]
names(lst) <- c("a","b","c")
lst$b
names(lst)

# Data frames are lists with the additional restriction that 
# each list element must have the same length.
# In other words, data frames, like matrices, must be "rectangular".
# Data frames are often used to represent expression or phenotype data.

df <- data.frame(Gene1=c(1,5,3), 
                 Gene2=c(2,4,7), Gene3=c(8,3,1))
df[1:2,]
df$Gene2
# apply a function across rows or columns
apply(df, 1, mean)
apply(df, 2, mean)
# If a data.frame has only numeric elements, it can be easily coerced
# to a matrix:
m <- as.matrix(df)
heatmap(m)

# Factors represent categorical variables in efficient form

groups <- c("red","green","red","orange","green")
groups.f <- factor(groups)
groups.f
levels(groups.f)
groups.f=="red"
groups.f[groups.f=="red"]

# Tables are *frequency* tables, representing occurrences or co-occurrences
# of variables in a single or multivariate distribution.
x <- sample(1:5, 50, replace=T)
table(x)
y <- sample(1:3, 50, replace=T)
table(x,y)

# prop.table (proportion table) finds probabilities or conditional probabilities
# for a probability distribution.
txy <- table(x,y)
prop.table(txy)
prop.table(txy,margin=1)
prop.table(txy,margin=2)

# Exercise!
# Before going further, think: how would you implement this marginalization?

# --------------------------------------------------------------------------

# Implementing prop.table marginalization:
pt <- table(x,y)
pt / apply(pt,1,sum) # row-wise
pt / apply(pt,2,sum) # column-wise

# Note that a similar technique is also used to standardize a gene expression matrix. 
# For each column:
# - subtract the mean expression from each gene
# - then divide by the standard deviation or root mean square
# - see ``?scale`` and R's built-in implementation by typing ``prop.table``

##################
# R package system
##################

# There are two main package systems of interest to biologists:
# CRAN, which has generic programming-related packages and statistics packages
# Bioconductor, which has biology-specific packages

# CRAN uses install.packages
install.packages("sos")
# "library" loads an installed package into memory
library(sos)
# "findFn" from the "sos" package lets you search for functions of interest
findFn("ReadAffy")

# Bioconductor must first be loaded from source
# This line must be run every R session that you intend to install BioC packages,
# unless you put it in ~/.Rprofile
source("http://bioconductor.org/biocLite.R")
#Windows users can find their R User (home) directory by:
Sys.getenv("R_USER")

# sourcing biocLite.R enables the use of the "biocLite" function, which installs
# packages from Bioconductor.
biocLite("WGCNA")

###################
# The R help system
###################

biocLite("affy")
library(affy)

# view help for a particular function
?read.table
# search help for a substring
??read
# to list all the defined variables in a package:
library(affy)
ls("package:affy")
# to view the source code for a function, simply type its name
ReadAffy

# To list or view vignettes (expository PDFs about packages)
vignette()
vignette("affy")


##############################
# Creating and using functions
##############################

# Creating your own function
PI <- 3.14
area <- function(radius) {
    # Functions can "capture" data external to the function. This is called a **closure**.
    # The last item in a function body is implicitly returned.
    2 * PI * radius
}
area(5)

# In the ordinary case, R functions are pass-by-value:
x <- 5
redefineX <- function(newX) {
    x <- newX
    print(x)
}
redefineX(7)
x

# However, a double arrow redefines a variable in the global scope. Use it sparingly.
redefineX <- function(newX) {
    x <<- newX
}
redefineX(7)
x

# The %apply functions (apply, lapply, sapply, etc.) apply functions to a collection
# apply we have already met above. For those familiar with
# functional programming languages, these functions are 
# analogous to the "map" function in Python, or "mapcar" in Lisp.

# (Note, however, that R has "Map", "Reduce", and "Filter"
# functions also).

doubleMe <- function(x) {
    x * 2
}
# sapply returns a vector:
sapply(1:5, doubleMe)
# lapply returns a list:
lapply(1:5, doubleMe)

# R has normal control-flow structures
if (5 > 7) { 
    print("A")
} else {
    print("B")
}

# R also has foreach and while loops.
# However, R has many optimizations for vectorization, 
# so use of these loops is discouraged unless it is necessary.
for (i in 1:5) {
    print(i)
}

i <- 1
while (i < 10) {
    print(i)
    i <<- i + 1
}

##########################
# Interacting with your PC
##########################

# Save and load arbitrary data structures to file
df <- data.frame(A=c(1,2,3), B=c("a","b","c"))
myData=list(myVector=c(1,2,3), 
            myData=df)

save(myData, file="myData.Rda")
load("myData.Rda")

write.csv(df, file="myData.csv")
df <- read.csv("myData.csv")

# Interact with the file system
getwd()
# Set working directory
# setwd("~") 
file.exists("myData.Rda")
unlink("myData.Rda")

# Inspect workspace
ls()

#########################################
# More operations on vectors and matrices
#########################################

# Sample 5 elements from uniform random distribution
v1 <- runif(5)
v2 <- runif(5)
v1
v2
pmin(v1,v2)
cor(v1,v2)

m1 <- matrix(runif(25), nr=5)
m2 <- matrix(runif(25), nr=5)

# Dot product
v1 %*% v2
m1 %*% v1
m1 %*% m2

# Outer product
outer(v1,v2)

arr <- outer(m1,v1)

class(arr)
# Note that this returns a 3-dimensional matrix-like structure,
# called an array. An array is implemented as a large vector,
# which stores additional information indicating the dimensions.

# In many situations, an array behaves like a vector or matrix.
# Basic arithmetic is the same:
arr * 3

arr ^ 2

# However, many functions which take both vectors and matrices,
# such as "cor", do not take arrays:

cor(arr, v1) # will produce an error

# Arrays can be constructed by the "array" function, and more
# information about them can be obtained by "?array".
# Arrays are seldom used in bioinformatics applications, but 
# can be advantageous for, e.g., time series data.

##################
# Basic statistics
##################

# Simulate 10000 trials of 100 "coin flips"
# (sampling random numbers from the binomial distribution)
nheads <- rbinom(10^4, 100, 0.5)

hist(nheads)
# Create a kernel density estimate (KDE) for this distribution
# (A KDE is a nonparametric estimate of a PDF from a sample. See
# ?density for available kernels). In simple terms, a KDE is 
# similar to a "smoothed" histogram.
dn <- density(nheads)
plot(dn)

# There are many distributions available in R: 
# binom, normal, beta, gamma, hyper(geometric), logis(tic), pois(son),
# and many more.

# Each distribution has 4 associated functions, identifiable by a prefix:
# d* - density function (pdf)
# p* - cumulative density function (cdf)
# q* - quantile function
# r* - random sampling function

# Example: standard normal

# PDF of standard normal evaluated at zero
dnorm(0)
# PDF of standard normal integrated from negative infinity to zero
# (i.e., CDF at zero)
pnorm(0) 

# Plot CDF and PDF together in the same window
x <- seq(-6,6,by=0.25)
par(mfrow=c(2,1))
plot(x, dnorm(x))
plot(x, pnorm(x))

dev.off()

# The quantile function q* is the inverse of the CDF.
# It returns the value of "x" for a given quantile.
# Thus, the domain of q* is [0,1]. Anything else will return NA.
qnorm(0.25)
qnorm(0.5)
qnorm(1.25)

# To view documentation on a particular distribution, call for help on any
# of the 4 distribution-associated functions.
?rpois

# And once again, r* :
hist(rnorm(1000))

########################################
# Basic Bioconductor for Expression Data
########################################

# Install the ALL package
biocLite("ALL")
# Load the ALL package
library(ALL)
# Load a dataset that comes with ALL package
data(ALL)

# ALL is a dataset containing microarray expression data and 
# phenotype data for some Acute Lymphoblastic Leukemia patients.
# It has already been partially preprocessed.
ALL

attributes(ALL)
# "pData" (or "phenoData") returns a data frame containing phenotype data 
# for each array 
pd <- pData(ALL)
head(pd)
tail(pd)

# "experimentData" shows metadata about the whole experiment
experimentData(ALL)

# The "annotation" function shows the probeset used for this microarray.
annotation(ALL)

# Perhaps most importantly, "exprs" returns a matrix of observed gene 
# expression data, with probes as rows and samples as columns.
m <- exprs(ALL)

# View a subset of the data. 
# Depending on the particular dataset, different normalization procedures
# will have been applied. Note that this one has been log-transformed.
# (Viewing "m" without subsetting first is unwise because it 
# is very large and printing it will be quite slow.)
m[1:5,1:5]

# View the dimensions of the expression matrix
dim(m)

# Plot a small heatmap. Note that R performs hierarchical clustering 
# on both axes for you.
heatmap(m[1:10,1:10])

# Quantile-quantile plot
qqplot(m[1,],m[2,])

# matplot is another useful way to quickly plot matrix data
first5probes <- t(m[1:5,])
matplot(first5probes, type="l")

# Tumor suppressor p53 is a gene that plays a key role in many cancers.
# One of its probes is "1939_at". Get its expression and key statistics
# about its distribution:
p53 <- m["1939_at",]
p53
quantile(p53)
summary(p53)

# Find the Pearson correlation between p53 and all probes. 
# "cor" can be used to find correlation between two vectors, two matrices,
# or a vector and a matrix. If matrices are used, variables are columns,
# not rows. Hence we must transpose the matrix with the "t" function before
# finding the correlation.
# We must also convert it from a matrix to vector for the next step.
cors <- cor(p53,t(m))
cors <- as.vector(cors)
names(cors) <- rownames(m)

# Show the most highly correlated probes. Note that p53 has two other probes:
# "1974_s_at" and "31618_at". What do these results tell you about the nature 
# of probes?
sorted.cors <- sort(cors, decreasing=T)
sorted.cors[1:10]
# "match" is like .index or .indexOf in Python/Java/C++
match(c("1974_s_at", "31618_at"), names(sorted.cors))

######################
# Additional Exercises
######################

# Try to do the following exercises yourself or with neighbors. 
# The facilitators will also be available for help. 
# Do not feel compelled to do them all -- choose the ones that interest you.

############
# Exercise 1
############

# a. Write a recursive function that, given an integer, computes the Fibonacci number for that integer.
# The definition of the Fibonacci sequence is:
# F(0) = 0
# F(1) = 1
# F(n) = F(n-2) + F(n-1)

# b. Install the "memoise" package from CRAN. 
# View the package documentation at: http://cran.r-project.org/web/packages/memoise/memoise.pdf
# Wrap the Fibonacci function above in a memoiser so that it is more efficient.

############
# Exercise 2
############

# Another common cancer-related gene is MYC. MYC has the following probes in the ALL dataset:
myc.probes <- c("1827_s_at", "1936_s_at", "1973_s_at", "37724_at")

# Plot all four probes together on a scatter plot.

# Which of these probes has the highest and lowest mean expression? 
# Which has the highest and lowest variance?

# View the expression of these four genes in a heatmap. How are they hierarchically clustered?

# Are these probes more highly correlated with one another than the p53 probes?

# The "cor" function has three kinds of correlation: pearson, spearman, and kendall.
# Which of these three methods of correlation causes these probes to have the 
# highest relative correlation with one another? 
