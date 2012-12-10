==============================
Introduction to R/Bioconductor
==============================

Welcome & Schedule
==================

.. include:: schedule.rst

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

Installing and using R packages:

- CRAN
- Bioconductor

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

Also note:

- Vectors are 1-indexed
- Scalars are just vectors of length 1

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

Apply a function across rows or columns:

.. code-block:: r

    > apply(df, 1, mean)
    [1] 3.666667 4.000000 3.666667
    > apply(df, 2, mean)
       Gene1    Gene2    Gene3 
    3.000000 4.333333 4.000000

R documentation system
======================

.. code-block:: r

    # view help for a particular function
    > ?read.table
    # search help for a substring
    > ??read

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

Bioconductor is a repository for biology-related R packages. Browse it at http://bioconductor.org/

.. code-block:: r

    # load Bioconductor system (must be done each session)
    > source("http://bioconductor.org/biocLite.R")
    # install the package
    > biocLite("affy")

Summary: Part One
=================

Part Two: Basic Data Analysis and Statistics in R
=================================================
