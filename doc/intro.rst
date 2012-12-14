==============================
Introduction to R/Bioconductor
==============================

.. topic:: Mikhail Dozmorov & Cory Giles

    Oklahoma Medical Research Foundation, OK, USA

.. image:: img/logo-R.png
    :scale: 175%
    :align: center

.. image:: img/bioconductor.jpg
    :scale: 175%
    :align: right

Welcome & Schedule
==================

TODO

What is R?
==========

- A programming language
- A data analysis and visualization environment
- A platform for new statistical techniques
- Widely used for microarray and sequencing analysis

A Brief History of R
====================

http://digitheadslabnotebook.blogspot.com/2010/01/r-type-system.html

- Started as S/SPlus
- Lisp influences

R versus other programming languages
====================================

- Algol family syntax (C/C++)
- Vectorized computation (Matlab)
- Dynamically typed (Python, Ruby, Perl)
- Strong libraries (Python)
- Excellent visualization tools (Mathematica)
- Functional (Lisp, OCaml, Haskell)

A whirlwhind tour of R
======================

Data types:

- vectors
- matrices
- lists
- data frames
- factors

Getting around in the R environment:

- Help and documentation system
- functions and classes
- I/O and string manipulation

Installing and using R packages

R data types: vector
====================

Creating:

.. code-block:: r
    
    > numbers <- 3:7
    > letters <- c("a","b","c")

Subscripting:

.. code-block:: r
    
    > numbers[1:3]
    [1] 3 4 5
    > numbers[c(1,3,5)]
    [1] 3 5 7

baz
===

Basic descriptive statistics:

.. code-block:: r

    > length(numbers)
    [1] 5
    > sum(numbers)
    [1] 25

See also: *prod*, *min*, *max*, *cumsum*, *cumprod*

R data types: vector
====================

Optional named vector indices:

.. code-block:: r
    
    > names(numbers) <- c("red","green","blue","orange","purple")
    > numbers
    red  green   blue orange purple 
    3      4      5      6      7
    > numbers[c("red","purple")]
    red purple
     3      7


Vectors must have a homogeneous type:

.. code-block:: r
    
    > myVector <- c(1, "foo", 2)
    > myVector
    [1] "1" "foo" "2"

Vectors are 1-indexed and scalars are just vectors of length 1

R data types: basic vectorized math
===================================

.. code-block:: r

    > v <- 1:5
    > v * 3
    [1] 3 6 9 12 15
    > v + 2
    [1] 3 4 5 6 7
    > v - 5
    [1] -4 -3 -2 -1 0
    > v / 2
    [1] 0.5 1.0 1.5 2.0 2.5

Exponentiation and modulus:

.. code-block:: r

    > v ^ 2
    [1] 1 4 9 16 25
    > v %% 2
    [1] 1 0 1 0 1

R data types: operations on two (or more) vectors
=================================================

Dot product:

.. code-block:: r

    > v %.% v
        [,1]
    [1,] 55

pairwise min max
cor

R data types: matrices
======================

- Essentially 2-dimensional vectors
- Also require homogeneous type

Creation:

.. code-block:: r
    
    > m <- matrix(1:9, nrow=3)
         [,1] [,2] [,3]
    [1,]    1    4    7
    [2,]    2    5    8
    [3,]    3    6    9

Subscripting:

.. code-block:: r
    
    > m[2:3,1:2]
         [,1] [,2]
    [1,]    2    5
    [2,]    3    6 

R data types: matrices
======================

Like vectors, named indices can be assigned and subscripted (with *rownames* and *colnames* instead of *names*)

.. code-block:: r
    
    > colnames(m) <- c("red", "green", "blue")
    > m.sub <- m[,c("red","blue")]
    > m.sub
        red blue
    [1,]   1    7
    [2,]   2    8
    [3,]   3    9

Transposition:

.. code-block:: r

    > t(m.sub)
         [,1] [,2] [,3]
    red     1    2    3
    blue    7    8    9

R data types: matrices
======================

Dimensions:

.. code-block:: r

    > nrow(m.sub)
    [1] 3
    > ncol(m.sub)
    [1] 2
    > dim(m.sub)
    [1] 2 3

The type of any object can be inspected using *class*:

.. code-block:: r
    
    > class(m)
    [1] "matrix"
    > class(1:5)
    [1] "numeric"
    > class("apple")
    [1] "character"

R data types: lists
===================

Lists are containers for *heterogeneous* objects.

.. code-block:: r

    > lst <- list(1,c("foo","bar"),c(3.5,5))
    > lst
    [[1]]
    [1] 1

    [[2]]
    [1] "foo" "bar"

    [[3]]
    [1] 3.5 5

    > lst[3]
    [[3]]
    [1] 3.5 5
    > lst[[3]]
    [1] 3.5 5

R data types: lists
===================

Named indexing:

.. code-block:: r

    > names(lst) <- c("a","b","c")
    > lst$b
    [1] "foo" "bar"

R data types: data frame
========================

- A *data frame* is a "square" list. 
- This is the data type usually used to manipulate experimental results.
- Many methods that apply to matrices (*nrow*, *ncol*, *dim*) also apply to data frames

.. code-block:: r

    > df <- data.frame(Gene1=c(1,5,3), 
        Gene2=c(2,4,7), Gene3=c(8,3,1))
    > df[1:2,]
      Gene1 Gene2 Gene3
    1     1     2     8
    2     5     4     3
    > df$Gene2
    [1] 2 4 7
    # apply a function across rows or columns
    > apply(df, 1, mean)
    [1] 3.666667 4.000000 3.666667
    > apply(df, 2, mean)
       Gene1    Gene2    Gene3 
    3.000000 4.333333 4.000000

R data types: factors
=====================

- Represents a categorical variable in computationally efficient form 
- A sort of hybrid between character and integer
- Similar to a C++ or Java **enumeration**

R data types: tables
====================

- AKA **frequency** tables, across one or more dimension

.. code-block:: r

    > x <- sample(1:5, 50, replace=T)
    > table(x)
    # Your results may be slightly different
     1  2  3  4  5 
     12 14  8 10  6
    > y <- sample(1:3, 50, replace=T)
    > table(x,y)
           y
    x   1 2 3
      1 5 4 3
      2 4 4 6
      3 4 1 3
      4 5 1 4
      5 0 3 3

R tables and Bayesian statistics
================================

- R tables can be used to find probabilities and conditional probabilities for variables
- *prop.table* with no arguments gives the joint probability distribution of two or more variables

.. code-block:: r

    > prop.table(table(x,y))
           y
    x      1    2    3
      1 0.10 0.08 0.06
      2 0.08 0.08 0.12
      3 0.08 0.02 0.06
      4 0.10 0.02 0.08
      5 0.00 0.06 0.06 

R tables and Bayesian statistics
================================

- *prop.table* takes a keyword argument, *margin*, which marginalizes the joint probability distribution across a dimension, returning a conditional probability table (CPT)

.. code-block:: r

    > prop.table(table(x,y), margin=1)
       y
    x           1         2         3
      1 0.4166667 0.3333333 0.2500000
      2 0.2857143 0.2857143 0.4285714
      3 0.5000000 0.1250000 0.3750000
      4 0.5000000 0.1000000 0.4000000
      5 0.0000000 0.5000000 0.5000000
    > prop.table(table(x,y), margin=2)
       y
    x            1          2          3
      1 0.27777778 0.30769231 0.15789474
      2 0.22222222 0.30769231 0.31578947
      3 0.22222222 0.07692308 0.15789474
      4 0.27777778 0.07692308 0.21052632
      5 0.00000000 0.23076923 0.15789474

- **Challenge: can you implement this marginalization yourself?** (Solution on next slide)

Manual marginalization of prop.table: solution
==============================================

.. code-block:: r

    > pt <- table(x,y)
    > pt / apply(pt,1,sum) # row-wise
        y
    x           1         2         3
      1 0.4166667 0.3333333 0.2500000
      2 0.2857143 0.2857143 0.4285714
      ... 
    > pt / apply(pt,2,sum) # column-wise
        y
    x            1          2          3
      1 0.27777778 0.30769231 0.15789474
      2 0.22222222 0.30769231 0.31578947
      ...


- This technique is also used to standardize a gene expression matrix. For each column:
    - subtract the mean expression from each gene
    - then divide by the standard deviation or root mean square
    - see ``?scale`` and R's built-in implementation by typing ``prop.table``

R documentation system
======================

.. code-block:: r

    # view help for a particular function
    > ?read.table
    # search help for a substring
    > ??read
    # to list all the defined variables in a package:
    > library(affy)
    > ls("package:affy")
    # to view the source code for a function simply type its name
    > ReadAffy

List or view vignettes:

.. code-block:: r
    
    > vignette()
    > vignette("affy")

R package system: CRAN
======================

CRAN is the repository for general (non biology-related) R packages. CRAN can be found at http://cran.r-project.org/

Searching for and installing new packages:

.. code-block:: r

    # Install a package from CRAN
    > install.packages("sos")
    > library(sos)
    # Search CRAN (will open a browser window)
    > findFn("ReadAffy")

R package system: Bioconductor
==============================

Bioconductor (BioC) is a repository for biology-related R packages. Browse it at http://bioconductor.org/

.. code-block:: r

    # load Bioconductor system (must be done each session)
    > source("http://bioconductor.org/biocLite.R")
    # install a specific package
    > biocLite("affy")
    # install the core set of packages
    > biocLite()

Summary: Part One
=================

- R is a dynamic programming language and data analysis environment
- R's type system does not distinguish between collections and scalars
- R's type system has several types specialized for analyzing datasets:
    - Vectors
    - Matrices
    - Lists
    - Data Frames
    - Factors
    - Tables
- General R packages can be installed from CRAN, and biology-related packages can be installed from Bioconductor

Part Two: Basic Data Analysis and Statistics in R
=================================================

- Statistical methods in R
- Power analysis and microarray experimental design

Statistics with R
=================

.. code-block:: r

    # First, load a dataset
    > library(ALL)
    > data(ALL)
    > expression <- exprs(ALL) 
    > myc <- expression["1827_s_at",]

    > mean(myc) # see also min, max, median
    8.372506

    > quantile(myc)
        0%      25%      50%      75%     100% 
    7.688711 8.210462 8.339376 8.557417 9.582462 

    > summary(myc)
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    7.689   8.210   8.339   8.373   8.557   9.582 


Experimental design
===================

Why do a microarray study?

- Class Comparison – given classes with membership known a priori, ﬁnd genes showing differences between classes.
- Class Prediction - build a model characterizing known classes, and use the model to predict the class status of future samples
- Class Discovery - identify subsets of samples based on their clustering behavior

Power analysis: Basics
======================

- Question: How many samples will be needed to achieve a desired false negative rate (FNR)?
- Depends on the statistical test that will be used to calculate differential expression (t-test, ANOVA, etc.)
- R has built-in `power.t.test`, 
- There is *always* a tradeoff between Type I and Type II errors

Power analysis: Parameters
==========================

- n <- Sample size to collect or available
- delta <- log2(FC) effect size
- sigma <- "noise" or variability within groups
- alpha <- type I error rate (FPR)
- power <- (1 - FNR)

For ANOVA:

- number of treatment groups

Uncategorized
=============

- The issue of correctly modeling ordered categorical data is complicated
- For example, Grade I-IV glioma: the simplest approach is to treat all intervals as equal
- **Think carefully** before you use this assumption!
