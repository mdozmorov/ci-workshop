###################################################
# Intermediate R
# Cory Giles
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

###################
# Grouping datasets
###################

data(iris)
iris

tapply(iris$Sepal.Length, iris$Species, FUN=mean)
aggregate(Sepal.Length ~ Species, data=iris, FUN=mean)
 and even:
aggregate(Sepal.Length * Sepal.Width ~ Species, 
          data=iris, FUN=mean)

# See also:
# * http://www.slideshare.net/jeffreybreen/grouping-summarizing-data-in-r
# * the **doBy** package. 

######################################
# Flexible data manipulation with plyr
######################################

# plyr naming conventions: 
# (input type) + (output type) + "ply"
# "a" = array
# "d" = data.frame
# "l" = list
# "m" = matrix
# "_" = no output returned

ddply(iris, 'Species', function(df) {
      c(Count=nrow(df), Sepal.Width=mean(df$Sepal.Width))
})

# Can group by more than one field
data(mtcars)
head(mtcars)
ddply(mtcars, .(cyl, gear), function(df) {
      c(Count=nrow(df), MPG=mean(df$mpg))
})

#########################################
# sqldf - manipulate data frames with SQL
#########################################

library(sqldf); data(iris)
head(iris,3)
sqldf("SELECT Species, AVG(Sepal_Length) 
      FROM iris GROUP BY Species")

##########################################
# SQL API - accessing SQL databases from R
##########################################

biocLite("GEOmetadb")
library(GEOmetadb)
# Download GEOmetadb file into current directory
getSQLiteFile()

install.packages("RSQLite")
library(RSQLite)
driver <- dbDriver("SQLite")
# For a network (MySQL, MSSQL, etc.) DB, you could
# specify host, port, etc, here
db <- dbConnect(driver, dbname="GEOmetadb.sqlite")
dbListTables(db)
# SELECT statement
dbGetQuery(db, "SELECT ID,gds,update_date FROM gds LIMIT 5")

# Create and modify an in-memory database
db <- dbConnect(SQLite(), dbname=":memory:")
dbGetQuery(db, "create table person (name text, age integer)")
dbBeginTransaction(conn)
people <- data.frame(name=c("ann","bill","charlie"),
                     age=c(33,54,27))
dbGetPreparedQuery(db, "insert into person values (?,?)",
                   bind.data=people)
dbCommit(db)
dbGetQuery(db, "select count(*) from person")
 
#######################
# Multicore programming
#######################

# Multicore is a simple way to thread R jobs
# (we have already seen this in the network session)

biocLite("multicore")
library(multicore)
mclapply(1:10, rnorm)

# Benefits: 
# * Multithreading, not multiprocessing, so spawning is fast
# * Parallelized function automatically has access to all workspace data (contra MP)

# Caveats:
# * To avoid performance loss from context switching, each job should be fairly large (good rule of thumb: > 500 ms/job)

# An alternative model for parallelism is called:
# MPI (Message-passing interface)
# - Multiprocessing, not multithreading
# - Disadvantage: Memory must be explicitly shared
# - Advantage: MPI can be run on multiple cores of same workstation or on multiple nodes

# Installing it requires platform-specific steps (you must
# install MPI headers). We will not show an example here,
# but remember if you have very large jobs, too large for one
# machine, the "Rmpi" package can be handy.

#############################################
# Warnings, assertions and exception handling
#############################################

mySqrt <- function(x) {
    if (is.na(x)) {
        # "throw" or "raise" an exception
        stop("mySqrt does not accept NAs")
    } else if (x < 0) {
        # Display a warning
        warning("Taking square root of a negative number...")
        complex(imaginary=sqrt(-x))
    } else {
        sqrt(x)
    }
}

mySqrt(5)

mySqrt(-5)

mySqrt(NA)

# Similar to "assert" in other languages
stopifnot(F)

stopifnot(T)

# A try-catch statement attempts to execute the block
# in the first argument. When it encounters a warning or error,
# if there is a handler, the handler is called instead of 
# normally executing the warning/error.

# Especially note the differences between examples 1 and 2.

tryCatch({
    result <- mySqrt(-5)
    print(result)
}, warning=function(w) { print("caught a warning") },
   error=function(e) { print("caught an error") })

tryCatch({
    result <- mySqrt(-5)
    print(result)
}, error=function(e) { print("caught an error") })

tryCatch({
    result <- mySqrt(NA)
    print(result)
}, warning=function(w) { print("caught a warning") },
   error=function(e) { print("caught an error") })

# The "finally" argument executes regardless of whether an 
# error was raised, and whether it was caught or uncaught.

tryCatch({
    result <- mySqrt(NA)
}, warning=function(w) { print("caught a warning") },
   error=function(e) { print("caught an error"),
   finally={ print("hello") })

# Beyond tryCatch, there is also a "naked" try statement, without
# a corresponding "catch", but it is rarely useful.

#############################
# Object-oriented programming
#############################

# An "object" is a collection of data and associated functions (sometimes called methods). 
# In other languages, objects are sometimes called classes or structs.
# A key idea of objects is to enforce "encapsulation" -- 
# in other words, you don't have to know the actual data structures in an object 
# in order to use it.

# Taken to an extreme, you don't access any data except via "getter" and "setter" methods (common in e.g,. Java).

# R's has two (rather ad hoc) OO systems -- **S3** and **S4**:

# S3 is simpler and more dynamic, but has less "enforcement". Single dispatch.
# S4 is more complex and is closer to OO in other languages. Multiple dispatch.

# (The nomenclature comes from the proprietary S language version numbers, which R clones)

############
# S3 Classes
############

# Many familiar R methods are S3 methods:

methods("plot")
methods(class="lm")

# S3 Classes:
# "Inheritance" is achieved in S3 by having a vector of classes:
fit <- glm(rpois(100,1) ~ 1, family="poisson")
class(fit)

# When a method is called, methods are searched in left-to-right order. 
methods("residuals")
methods("model.matrix")

# Creating a custom class and associated class methods is easy:
bill <- list(name="bill", age=35, gender="M")
class(bill) <- "person"
print.person <- function(p) {
    cat("name:", p$name, "\n")
    cat("age:", p$age, "\n")
    cat("gender:", p$gender, "\n")
}
print(bill)

# Limitations of S3 methods:

# In the previous example, note that a "person" is not required to actually 
# have the "name", "age", or "gender" attributes. Also, "age" is not required
# to be a integer/numeric, name is not required to be character, etc...
# This is known as a dynamic type system, and can easily lead to runtime errors.

# Also, A S3 class method is dispatched (only) on its first argument. 
# Thus, you cannot easily write a method that enforces types on all arguments as in C:

"
char nthchar(char* str, int idx) {
    return str[idx];
}  
"

############
# S4 Classes
############ 

rep <- representation(name="character", 
                      age="numeric", gender="character")
# note that representation is simply 
# implemented as a list
rep
setClass("Person", rep)

# Instantiating a S4 class
bill <- new("Person", name="Bill", age=35, gender="M")
bill
bill <- new("Person", name="Bill", age=35, gender="M")
bill@gender
bill@age <- 30

# Creating methods:
setGeneric("birth.year", function(p) standardGeneric("birth.year"))
setMethod("birth.year", signature("Person"), function(p) {
          current.year <- as.numeric(substr(Sys.time(), 1, 4))
          current.year - p@age})
birth.year(bill)

###########################
# Creating a custom package
###########################

# Creating custom R packages can be useful for distributing
# your work to others, or even to just keep in-house code
# modular and understandable.

# The basic package creation function is:

package.skeleton()

# We will discuss this in detail in class, time permitting,
# including such topics as:
# - basic package layout
# - including native code
# - literate programming with knitr or Sweave

#################
# Version control
#################

# Although not strictly related to R, version control is an extremely important tool for software development and reproducible research.

# If there is time, we will discuss the basics of git, a popular distributed version control system (DVCS).

###########
# Exercises
###########

# 1. Making use of R's object-oriented programming and exception handling abilities, choose the workflow from all of our sessions that most interested you, and turn it into a small R package for your own use. 
