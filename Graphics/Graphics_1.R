###################################################
# Graphics
# Mikhail Dozmorov
# Workshop at Rajiv Gandhi Center for Biotechnology
# January 16-18, 2013
# Material is public domain
###################################################

#####################
# Graphics overview
##################### 


source("http://bioconductor.org/biocLite.R") # Import biocLite() function into R environment
##Graphical Procedures
# R provides comprehensive graphics utilities for visualizing and exploring scientific data. In addition, several powerful graphics environments extend these utilities. These include the grid, lattice and ggplot2 packages. The grid package is part of R's base distribution. It provides the low-level infrastructure for many graphics packages, including lattice and ggplot2. Extensive and systematic information on graphics utilities in R can be found on the Graphics Task Page [http://cran.r-project.org/web/views/Graphics.html] and Paul Murrell's book R Graphics [http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html].

## Base Graphics
# The following list provides an overview of some very useful plotting functions in R's base graphics. To get familiar with their usage, it is recommended to carefully read their help documentation with ?myfct as well as the help on the functions ?plot and ?par. 
# plot: generic x-y plotting
# barplot: bar plots
# boxplot: box-and-whisker plot
# hist: histograms
# pie: pie charts
# dotchart: cleveland dot plots
# image, heatmap, contour, persp: functions to generate image-like plots
# qqnorm, qqline, qqplot: distribution comparison plots
# pairs, coplot: display of multivariant data

## Lattice  [ Manuals: lattice [http://lmdvr.r-forge.r-project.org/figures/figures.html], Intro [http://www.his.sunderland.ac.uk/~cs0her/Statistics/UsingLatticeGraphicsInR.htm] ] 
# The lattice package developed by Deepayan Sarkar implements in R the Trellis graphics system from S-Plus. The environment greatly simplifies many complicated high-level plotting tasks, such as automatically arranging complex graphical features in one or several plots. The syntax of the package is similar to R's base graphics; however, high-level lattice functions return an object of class "trellis", that can be either plotted directly or stored in an object. The command library(help=lattice) will open a list of all functions available in the lattice package, while ?myfct and example(myfct) can be used to access and/or demo their documentation. Important functions for accessing and changing global parameters are: ?lattice.options and ?trellis.device.

## ggplot2 [ Manuals: ggplot2 [http://ggplot2.org/], Intro [http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html]
# ggplot2 is another more recently developed graphics system for R, based on the grammar of graphics theory. The environment streamlines many graphics routines for the user to generate with minimum effort complex multi-layered plots. Its syntax  is centered around the main ggplot function, while the convenience function qplot provides many shortcuts. The ggplot function accepts two arguments: the data set to be plotted and the corresponding aesthetic mappings provided by the aes function. Additional plotting parameters such as geometric objects (e.g. points, lines, bars) are passed on by appending them with '+' as separator. For instance,  the following command will generate a scatter plot for the first two columns of the iris data frame: ggplot(iris, aes(iris[,1], iris[,2])) + geom_point(). A list of the available geom_* functions can be found on [http://ggplot2.org/]. The settings of the plotting theme can be accessed with the command theme_get(). Their settings can be changed with the opts()function. To benefit from the many convenience features built into ggplot2, the expected input data class is usually a data frame where all labels for the plot are provided by the column titles and/or grouping factors in additional column(s).

## The following graphics sections demonstrate how to generate different types of plots first with R's base graphics device and then with the lattice and ggplot2 packages.   

## Scatter Plots
# Base Graphics
y <- as.data.frame(matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))) # Generates a sample data set.
plot(y[,1], y[,2]) # Plots two vectors (columns) in form of a scatter plot against each other.
plot(y[,1], y[,2], type="n", main="Plot of Labels"); text(y[,1], y[,2], rownames(y)) # Prints instead of symbols the row names.
plot(y[,1], y[,2], pch=20, col="red", main="Plot of Symbols and Labels"); text(y[,1]+0.03, y[,2], rownames(y)) # Plots both, symbols plus their labels.
grid(5, 5, lwd = 2) # Adds a grid to existing plot.
op <- par(mar=c(8,8,8,8), bg="lightblue"); plot(y[,1], y[,2], type="p", col="red", cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1, lwd=4, pch=20, xlab="x label", ylab="y label", main="My Main", sub="My Sub"); par(op) # This command demonstrates the usage of the most important plotting parameters. The 'mar' argument specifies the margin sizes around the plotting area in this order: c(bottom, left, top, right). The color of the plotted symbols can be controled with the 'col' argument. The plotting symbols can be selected with the 'pch' argument, while their size is controlled by the 'lwd' argument. A selection palette for 'pch' plotting symbols can be opened with the command 'example(points)'. As alternative, one can plot any character string by passing it on to 'pch', e.g.: pch=".". The font sizes of the different text components can be changed with the 'cex.*' arguments. Please consult the '?par' documentation for more details on fine tuning R's plotting behavior.
plot(y[,1], y[,2]); myline <- lm(y[,2]~y[,1], data=y[,1:2]); abline(myline, lwd=2) # Adds a regression line to the plot.
summary(myline) # Prints summary about regression line.
plot(y[,1], y[,2], log="xy") # Generates the same plot as above, but on log scale.
plot(y[,1], y[,2]); text(y[1,1], y[1,2], expression(sum(frac(1,sqrt(x^2*pi)))), cex=1.3) # Adds a mathematical formula to the plot.
plot(y) # Produces all possible scatter plots for all-against-all columns in a matrix or a data frame. The column headers of the matrix or data frame are used as axis titles.
pairs(y) # Alternative way to produce all possible scatter plots for all-against-all columns in a matrix or a data frame.
# biocLite("scatterplot3d")
library(scatterplot3d)
scatterplot3d(y[,1:3], pch=20, color="red") # Plots a 3D scatter plot for first three columns in y.
# biocLite("geneplotter")
library(geneplotter)
smoothScatter(y[,1], y[,2]) # Same as above, but generates a smooth scatter plot that shows the density of the data points.

# lattice
# biocLite("lattice")
library(lattice)
xyplot(1:10 ~ 1:10) # Simple scatter plot. 
xyplot(1:10 ~ 1:10 | rep(LETTERS[1:5], each=2), as.table=TRUE) # Plots subcomponents specified by grouping vector after '|' in separate panels. The argument as.table controls the order of the panels. 
myplot <- xyplot(Petal.Width ~ Sepal.Width | Species , data = iris); print(myplot) # Assigns plotting function to an object and executes it.
xyplot(Petal.Width ~ Sepal.Width | Species , data = iris, layout = c(1, 3, 1)) # Changes layout of individual plots.
# Change plotting parameters
show.settings() # Shows global plotting parameters in a set of sample plots.
default <- trellis.par.get(); mytheme <- default; names(mytheme) # Stores the global plotting parameters in list mytheme and prints its component titles. 
mytheme["background"][[1]][[2]] <- "grey" # Sets background to grey
mytheme["strip.background"][[1]][[2]] <- "transparent" # Sets background of title bars to transparent.
trellis.par.set(mytheme) # Sets global parameters to 'mytheme'.
show.settings() # Shows custom settings.
xyplot(1:10 ~ 1:10 | rep(LETTERS[1:5], each=2), as.table=TRUE, layout=c(1,5,1), col=c("red", "blue"))
trellis.par.set(default)

# ggplot2
# biocLite("ggplot2")
library(ggplot2)
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point() # Plots two vectors (columns) in form of a scatter plot against each other.
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point(aes(color = Species), size=4) # Plots larger dots and colors them with default color scheme.
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point(aes(color = Species), size=4) + ylim(2,4) + xlim(4,8) + scale_color_manual(values=rainbow(10)) # Colors dots with custom scheme.
ggplot(iris, aes(Sepal.Length, Sepal.Width, label=1:150)) + geom_text() + opts(title = "Plot of Labels") # Print instead of symbols a set of custom labels and adds a title to the plot.
ggplot(iris, aes(Sepal.Length, Sepal.Width, label=1:150)) + geom_point() + geom_text(hjust=-0.5, vjust=0.5) # Prints both symbols and labels.
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point() + opts(panel.background=theme_rect(fill = "white", colour = "black")) # Changes background color to white.
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point() + stat_smooth(method="lm", se=FALSE) # Adds a regression line to the plot.
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_point() + coord_trans(x = "log2", y = "log2") # Plots axes in log scale. 

## Line Plots
# Base Graphics
dev.off() # Clear display
y <- as.data.frame(matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))) # Generates a sample data set.
plot(y[,1], type="l", lwd=2, col="blue") # Plots a line graph instead.
split.screen(c(1,1)); plot(y[,1], ylim=c(0,1), xlab="Measurement", ylab="Intensity", type="l", lwd=2, col=1); for(i in 2:length(y[1,])) { screen(1, new=FALSE); plot(y[,i], ylim=c(0,1), type="l", lwd=2, col=i, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n") }; close.screen(all=TRUE) # Plots a line graph for all columns in data frame 'y'. The 'split.screen()' function is used in this example in a for loop to overlay several line graphs in the same plot. More details on this topic are provided in the 'Arranging Plots' section.
# A very nice line plot function for time series data is available in the Mfuzz library.

# lattice
xyplot(Sepal.Length ~ Sepal.Width | Species, data=iris, type="a", layout=c(1,3,1))
parallel(~iris[1:4] | Species, iris) # Plots data for each species in iris data set in separate line plot.
parallel(~iris[1:4] | Species, iris, horizontal.axis = FALSE, layout = c(1, 3, 1)) # Changes layout of plot.

# ggplot2
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_line(aes(color=Species), size=1) # Plots lines for three samples in the 'Species' column in a single plot.
ggplot(iris, aes(Sepal.Length, Sepal.Width)) + geom_line(aes(color=Species), size=1) + facet_wrap(~Species, ncol=1)  # Plots three line plots, one for each sample in Species column. 

## Bar Plots
# Base Graphics
y <- as.data.frame(matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))) # Generates a sample data set.
barplot(as.matrix(y[1:4,]), ylim=c(0,max(y[1:4,])+0.1), beside=T) # Generates a bar diagram for the first four rows in y. The barplot() function expects as input format a matrix, and it uses by default the column titles for labeling the bars. 
text(labels=round(as.vector(as.matrix(y[1:4,])),2), x=seq(1.5, 13, by=1)+sort(rep(c(0,1,2), 4)), y=as.vector(as.matrix(y[1:4,]))+0.02) # Adds corresponding values on top of each bar. 
ysub <- as.matrix(y[1:4,]); myN <- length(ysub[,1]); mycol1 <- gray(1:(myN+1)/(myN+1))[-(myN+1)]; mycol2 <- sample(colors(),myN); barplot(ysub, beside=T, ylim=c(0,max(ysub)*1.2), col=mycol2, main="Bar Plot", sub="data: ysub"); legend("topright", legend=row.names(ysub), cex=1.3, bty="n", pch=15, pt.cex=1.8, col=mycol2, ncol=myN) # Generates a bar plot with a legend for the first four rows of the input matrix 'y'. To plot the diagram in black and white, set the 'col' argument to 'col=col1". The argument 'ncol' controls the number of columns that are used for printing the legend.
par(mar=c(10.1, 4.1, 4.1, 2.1)); par(xpd=TRUE); barplot(ysub, beside=T, ylim=c(0,max(ysub)*1.2), col=mycol2, main="Bar Plot"); legend(x=4.5, y=-0.3, legend=row.names(ysub), cex=1.3, bty="n", pch=15, pt.cex=1.8, col=mycol2, ncol=myN) # Same as above, but places the legend below the bar plot. The arguments 'x' and 'y' control the placement of the legend and the 'mar' argument specifies the margin sizes around the plotting area in this order: c(bottom, left, top, right).
bar <- barplot(x <- abs(rnorm(10,2,1)), names.arg = letters[1:10], col="red", ylim=c(0,5)); stdev <- x/5; arrows(bar, x, bar, x + stdev, length=0.15, angle = 90); arrows(bar, x, bar, x + -(stdev), length=0.15, angle = 90) # Creates bar plot with standard deviation bars. The example in the 'barplot' documentation provides a very useful outline of this function.

# lattice
y <- matrix(sample(1:10, 40, replace=TRUE), ncol=4, dimnames=list(letters[1:10], LETTERS[1:4])) # Creates a sample data set.
barchart(y, auto.key=list(adj = 1), freq=T, xlab="Counts", horizontal=TRUE, stack=FALSE, groups=TRUE) # Plots all data in a single bar plot. The TRUE/FALSE arguments control the layout of the plot.
barchart(y, col="grey", layout = c(2, 2, 1), xlab="Counts", as.table=TRUE, horizontal=TRUE, stack=FALSE, groups=FALSE) # Plots the different data components in separate bar plots. The TRUE/FALSE and layout arguments allow to rearrange the bar plots in different ways. 

# ggplot2
# (A) Sample Set: the following transforms the iris data set into a ggplot2-friendly format   
iris_mean <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=mean) # Calculates the mean values for the aggregates given by the Species column in the iris data set. 
iris_sd <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=sd) # Calculates the standard deviations for the aggregates given by the Species column in the iris data set. 
convertDF <- function(df=df, mycolnames=c("Species", "Values", "Samples")) { myfactor <- rep(colnames(df)[-1], each=length(df[,1])); mydata <- as.vector(as.matrix(df[,-1])); df <- data.frame(df[,1], mydata, myfactor); colnames(df) <- mycolnames; return(df) } # Defines function to convert data frames into ggplot2-friendly format.
df_mean <- convertDF(iris_mean, mycolnames=c("Species", "Values", "Samples")) # Converts iris_mean. 
df_sd <- convertDF(iris_sd, mycolnames=c("Species", "Values", "Samples")) # Converts iris_sd.
limits <- aes(ymax = df_mean[,2] + df_sd[,2], ymin=df_mean[,2] - df_sd[,2]) # Define standard deviation limits.

# (B) Bar plots of data stored in df_mean 
ggplot(df_mean, aes(Samples, Values, fill = Species)) + geom_bar(position="dodge") # Plots bar sets defined by 'Species' column next to each other.
ggplot(df_mean, aes(Samples, Values, fill = Species)) + geom_bar(position="dodge") + coord_flip() + opts(axis.text.y=theme_text(angle=0, hjust=1)) # Plots bars and labels sideways. 
ggplot(df_mean, aes(Samples, Values, fill = Species)) + geom_bar(position="stack") # Plots same data set as stacked bars.
ggplot(df_mean, aes(Samples, Values)) + geom_bar(aes(fill = Species)) + facet_wrap(~Species, ncol=1) # Plots data sets below each other.
ggplot(df_mean, aes(Samples, Values, fill = Species)) + geom_bar(position="dodge") + geom_errorbar(limits, position="dodge") # Generates the same plot as before, but with error bars.

# (C) Customizing colors 
library(RColorBrewer); display.brewer.all() # Select a color scheme and pass it on to 'scale_*' arguments.
ggplot(df_mean, aes(Samples, Values, fill=Species, color=Species)) + geom_bar(position="dodge") + geom_errorbar(limits, position="dodge") + scale_fill_brewer(palette="Greys") + scale_color_brewer(palette = "Greys") # Generates the same plot as before, but with grey color scheme.  
ggplot(df_mean, aes(Samples, Values, fill=Species, color=Species)) + geom_bar(position="dodge") + geom_errorbar(limits, position="dodge") + scale_fill_manual(values=c("red", "green3", "blue")) + scale_color_manual(values=c("red", "green3", "blue")) # Uses custom colors passed on as vectors.


## Pie Charts
# Base Graphics
y <- table(rep(c("cat", "mouse", "dog", "bird", "fly"), c(1,3,3,4,2))) # Creates a sample data set.
pie(y, col=rainbow(length(y), start=0.1, end=0.8), main="Pie Chart", clockwise=T) # Plots a simple pie chart.
pie(y, col=rainbow(length(y), start=0.1, end=0.8), labels=NA, main="Pie Chart", clockwise=T); legend("topright", legend=row.names(y), cex=1.3, bty="n", pch=15, pt.cex=1.8, col=rainbow(length(y), start=0.1, end=0.8), ncol=1) # Same pie chart as above, but with legend.

# ggplot2
df <- data.frame(variable=rep(c("cat", "mouse", "dog", "bird", "fly")), value=c(1,3,3,4,2)) # Creates a sample data set.
ggplot(df, aes(x = "", y = value, fill = variable)) + geom_bar(width = 1) + coord_polar("y", start=pi / 3) + opts(title = "Pie Chart") # Plots pie chart.
ggplot(df, aes(x = variable, y = value, fill = variable)) + geom_bar(width = 1) + coord_polar("y", start=pi / 3) + opts(title = "Pie Chart") # Plots wind rose pie chart.

## Heatmaps
# Base Graphics
y <- matrix(rnorm(50), 10, 5, dimnames=list(paste("g", 1:10, sep=""), paste("t", 1:5, sep="")))
image(y) 
# lattice
library(lattice); library(gplots)
y <- lapply(1:4, function(x) matrix(rnorm(50), 10, 5, dimnames=list(paste("g", 1:10, sep=""), paste("t", 1:5, sep=""))))

## Plot single heatmap:
levelplot(y[[1]])

## Arrange several heatmaps in one plot
x1 <- levelplot(y[[1]], col.regions=colorpanel(40, "darkblue", "yellow", "white"), main="colorpanel")
x2 <- levelplot(y[[2]], col.regions=heat.colors(75), main="heat.colors")
x3 <- levelplot(y[[3]], col.regions=rainbow(75), main="rainbow")
x4 <- levelplot(y[[4]], col.regions=redgreen(75), main="redgreen")
print(x1, split=c(1,1,2,2))
print(x2, split=c(2,1,2,2), newpage=FALSE)
print(x3, split=c(1,2,2,2), newpage=FALSE)
print(x4, split=c(2,2,2,2), newpage=FALSE)

## Venn Diagrams: Computation of Venn intersects for 2-20 samples and plotting of 2-5 way Venn diagrams.
# Computation of Venn intersects
source("overLapper.R") # Imports several functions from the overLapper.R script for computing Venn intersects and plotting Venn diagrams (old version: vennDia.R). These functions are relatively generic and scalable by supporting the computation of Venn intersects of 2-20 or more samples. The upper limit around 20 samples is unavoidable because the complexity of Venn intersects increases exponentially with the sample number n according to this relationship: (2^n) - 1. A useful feature of the actual plotting step is the possiblity to combine the counts from several Venn comparisons with the same number of test sets in a single Venn diagram. 
#The overall workflow of the method is to first compute for a list of samples sets their Venn intersects using the overLapper function, which organizes the result sets in a list object. Subsequently, the Venn counts are computed and plotted as bar or Venn diagrams. The current implementation of the plotting function, vennPlot, supports Venn diagrams for 2-5 sample sets. To analyze larger numbers of sample sets, the Intersect Plot methods often provide reasonable alternatives. These methods are much more scalable than Venn diagrams, but lack their restrictive intersect logic. Additional Venn diagram resources are provided by limma, gplots, vennerable, eVenn, VennDiagram, shapes, C Seidel [http://www.pangloss.com/seidel/Protocols/venn.cgi] and Venny [http://bioinfogp.cnb.csic.es/tools/venny/index.html. 
setlist <- list(A=sample(letters, 18), B=sample(letters, 16), C=sample(letters, 20), D=sample(letters, 22), E=sample(letters, 18), F=sample(letters, 22, replace=T)) # To work with the overLapper function, the sample sets (here six) need to be stored in a list object where the different compontents are named by unique identifiers, here 'A to F'. These names are used as sample labels in all subsequent data sets and plots.
sets <- read.delim("sets.txt"); setlistImp <- lapply(colnames(sets), function(x) as.character(sets[sets[,x]!="", x])); names(setlistImp) <- colnames(sets) # Example how a list of test sets can be imported from an external table file stored in tab delimited format. Such a file can be easily created from a spreadsheet program, such as Excel. As a reminder, copy & paste from external programs into R is also possible (see read.delim function).
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets"); OLlist; names(OLlist) # With the setting type="vennsets", the overLapper function computes all Venn Intersects for the six test samples in setlist and stores the results in the Venn_List component of the returned OLlist object. By default, duplicates are removed from the test sets. The setting keepdups=TRUE will retain duplicates by appending a counter to each entry. When assigning the value "intersects" to the type argument then the function will compute Regular Intersects instead of Venn Intersects. Both analyses return a present-absent matrix in the Intersect_Matrix component of OLlist. Each overlap set in the Venn_List data set is labeled according to the sample names provided in setlist. For instance, the composite name 'ABC' indicates that the entries are restricted to A, B and C. The separator used for naming the intersect samples can be specified under the sep argument. By adding the argument cleanup=TRUE, one can minimize formatting issues in the sample sets. This setting will convert all characters in the sample sets to upper case and remove leading/trailing spaces.

# Bar plot of Venn counts
olBarplot(OLlist=OLlist, horiz=T, las=1, cex.names=0.6, main="Venn Bar Plot") # Generates a bar plot for the Venn counts of the six test sample sets. In contrast to Venn diagrams, bar plots scale to larger numbers of sample sets. The layout of the plot can be adjusted by changing the default values of the argument: margins=c(4,10,3,1). The minimum number of counts to consider in the plot can be set with the mincount argument (default is 0). The bars themselves are colored by complexity levels using the default setting: mycol=OLlist$Complexity_Levels.
# 2-way Venn diagrams
setlist2 <- setlist[1:2]; OLlist2 <- overLapper(setlist=setlist2, sep="_", type="vennsets"); OLlist2$Venn_List; counts <- sapply(OLlist2$Venn_List, length); vennPlot(counts=counts) # Plots a non-proportional 2-way Venn diagram. The main graphics features of the vennPlot function can be controlled by the following arguments (here with 2-way defaults): mymain="Venn Diagram": main title; mysub="default": subtitle; ccol=c("black","black","red"): color of counts; lcol=c("red","green"): label color; lines=c("red","green"): line color; mylwd=3: line width; ccex=1.0: font size of counts; lcex=1.0: font size of labels. Note: the vector lengths provided for the arguments ccol, lcol and lines should match the number of their corresponding features in the plot, e.g. 3 ccol values for a 2-way Venn diagram and 7 for a 3-way Venn diagram. The argument setlabels allows to provide a vector of custom sample labels. However, assigning the proper names in the original test set list is much more effective for tracking purposes. 
# 3-way Venn diagrams
setlist3 <- setlist[1:3]; OLlist3 <- overLapper(setlist=setlist3, sep="_", type="vennsets"); counts <- list(sapply(OLlist3$Venn_List, length), sapply(OLlist3$Venn_List, length)); vennPlot(counts=counts, mysub="Top: var1; Bottom: var2", yoffset=c(0.3, -0.2)) # Plots a non-proportional 3-way Venn diagram. The results from several Venn comparisons can be combined in a single Venn diagram by assigning to the count argument a list with several count vectors. The positonal offset of the count sets in the plot can be controlled with the yoffset argument. The argument setting colmode=2 allows to assign different colors to each count set. For instance, with colmode=2 one can assign to ccol a color vector or a list, such as ccol=c("blue", "red") or ccol=list(1:8, 8:1).
# 4-way Venn diagrams
setlist4 <- setlist[1:4]; OLlist4 <- overLapper(setlist=setlist4, sep="_", type="vennsets"); counts <- list(sapply(OLlist4$Venn_List, length), sapply(OLlist4$Venn_List, length)); vennPlot(counts=counts, mysub="Top: var1; Bottom: var2", yoffset=c(0.3, -0.2)) # Plots a non-proportional 4-way Venn diagram. The setting type="circle" returns a incomplete 4-way Venn diagram as circles. This representation misses two overlap sectors, but is sometimes easier to navigate than the default ellipse version.
# 5-way Venn diagrams
setlist5 <- setlist[1:5]; OLlist5 <- overLapper(setlist=setlist5, sep="_", type="vennsets"); counts <- sapply(OLlist5$Venn_List, length); vennPlot(counts=counts, ccol=c(rep(1,30),2), lcex=1.5, ccex=c(rep(1.5,5), rep(0.6,25),1.5)) # Plots a non-proportional 5-way Venn diagram.

# Export and other utilities
OLexport <- as.matrix(unlist(sapply(OLlist5[[4]], paste, collapse=" "))); write.table(OLexport, file="test.xls", col.names=F, quote=F, sep="\t") # Exports intersect data in tabular format to a file.
OLexport <- data.frame(Venn_Comp=rep(names(OLlist5[[4]]), sapply(OLlist5[[4]], length)), IDs=unlist(OLlist5[[4]])); write.table(OLexport, file="test.xls", row.names=F, quote=F, sep="\t") # Same as above, but exports to an alternative tabular format.
tapply(counts, OLlist5[[3]], function(x) rev(sort(x))) # Sorts the overlap results within each complexity level by their size. This allows to identify the sample set combinations with the largest intersect within each complexity level.
sapply(names(setlist), function(x) table(setlist[[x]])[table(setlist[[x]])!=1]) # Command to identify and count duplicated objects in the original sample set object setlist. In the given example, only set 'F' contains duplications. Their frequency is provided in the result.
vennPlot(counts, mymain="", mysub="", ccol="white", lcol="white") # Returns an empty Venn diagram without counts or labels.

# Typical analysis routine for sets of differentially expressed genes (DEGs)
ratio <- matrix(sample(seq(-5, 5, by=0.1), 100, replace=T), 100, 4, dimnames=list(paste("g", 1:100, sep=""), paste("DEG", 1:4, sep="")), byrow=T) # Creates a sample matrix of gene expression log2 ratios. This could be any data type!
setlistup <- sapply(colnames(ratio), function(x) rownames(ratio[ratio[,x]>=1,])); setlistdown <- sapply(colnames(ratio), function(x) rownames(ratio[ratio[,x]<=-1,])) # Identifies all genes with at least a two fold up or down regulation and stores the corresponding gene identifiers in setlistup and setlistdown, respectively..
OLlistup <- overLapper(setlist=setlistup, sep="_", type="vennsets"); OLlistdown <- overLapper(setlist=setlistdown, sep="_", type="vennsets"); counts <- list(sapply(OLlistup$Venn_List, length), sapply(OLlistdown$Venn_List, length)); vennPlot(counts=counts, ccol=c("red", "blue"), colmode=2, mysub="Top: DEG UP; Bottom: DEG Down", yoffset=c(0.3, -0.2)) # Performs Venn analysis for the four sets stored in setlistup and setlistdown. The argument setting colmode=2 allows to assign different colors to each count set. For instance, with colmode=2 one can assign to ccol a color vector or a list, such as ccol=c("blue", "red") or ccol=list(1:8, 8:1).

## Intersect Plots: Collection of methods for analyzing and visualizing intersect relationships among large numbers of sample sets.
# All Possible Intersects: for comprehensive overlap analyses of 2-20 or more sample sets
source("overLapper.R") # To compute all possible intersects among more than two samples at the same time - such as common in AB, ABC, ABCD and so on - the overLapper function can be used. Note: processing more than 20 sets with this function can take a long time, because the complexity of all possible intersects increases with the number of sample sets n according to this relationship: (2^n) - 1 (includes self comparisons). This means there are 63 possible intersects for n=6, 1,023 for n=10, 32,767 for n=15 and 33,554,431 for n=25. With the current implementation, the computation time is about 0.1 second for n=10 and 2.5 seconds for n=15.
setlist <- lapply(1:6, function(x) sample(letters, 20)); names(setlist) <- LETTERS[seq(along=setlist)] # To work with the overLapper function, the sample sets need to be stored in a list object where the different compontents are named by unique identifiers: A, B, C, ..., Z.
OLlist <- overLapper(setlist=setlist, complexity=1:length(setlist), sep="-", type="intersects"); OLlist; names(OLlist) # Computes all possible intersects for the samples stored in 'setlist'. The results are returned as a list where each overlap component is labeled by the corresponding sample names. For instance, the name 'A-B-C' indicates that the entries are common in samples A, B and C. The seperator used for naming the intersect samples can be specified under the 'sep' argument. The complexity level range to consider can be controlled with the 'complexity' argument, which can have values from 2 to the total number of samples. 
OLlist[[2]]; OLlist[[3]] # Returns the corresponding intersect matrix and complexity levels.
OLexport <- as.matrix(unlist(sapply(OLlist[[4]], paste, collapse=" "))); write.table(OLexport, file="test.xls", col.names=F, quote=F, sep="\t") # Exports intersect data in tabular format to a file.
counts <- sapply(OLlist[[4]], length) # Counts the number of elements in each intersect component.
tapply(counts, OLlist[[3]], function(x) rev(sort(x))) # Sorts the overlap results within each complexity level by their size. This allows to identify the sample set combinations with the largest intersect within each complexity level.
x11(height=12, width=8); olBarplot(OLlist=OLlist, horiz=T, las=1, cex.names=0.6, main="Intersect Bar Plot") # Plots the counts as bar diagram. More details on this function are provided in the Venn diagram section. 

## Histograms
#Base Graphics
x <- rnorm(100); hist(x, freq=FALSE); curve(dnorm(x), add=TRUE) # Plots a random distribution as histogram and overlays curve ('freq=F': ensures density histogram; 'add=T': enables 'overplotting'). 
plot(x<-1:50, dbinom(x,size=50,prob=.33), type="h") # Creates pin diagram for binomial distribution with n=50 and p=30.

# ggplot2 (See http://docs.ggplot2.org/current/geom_histogram.html for more details)
ggplot(iris, aes(x=Sepal.Width)) + geom_histogram(aes(fill = ..count..), binwidth=0.2) # Plots histogram for second column in 'iris' data set
ggplot(iris, aes(x=Sepal.Width)) + geom_histogram(aes(y = ..density.., fill = ..count..), binwidth=0.2) + geom_density() # Density representation.

# Density Plots
plot(density(rnorm(10)), xlim=c(-2,2), ylim=c(0,1), col="red") # Plots a random distribution in form of a density plot. 
split.screen(c(1,1)); plot(density(rnorm(10)), xlim=c(-2,2), ylim=c(0,1), col="red"); screen(1, new=FALSE); plot(density(rnorm(10)), xlim=c(-2,2), ylim=c(0,1), col="green", xaxt="n", yaxt="n", ylab="", xlab="", main="",bty="n"); close.screen(all = TRUE) # Several density plots can be overlayed with the 'split.screen' function. In ths example two density plots are plotted at the same scale ('xlim/ylim') and overlayed on the same graphics device using 'screen(1,new=FALSE))'.

# Box Plots
y <- as.data.frame(matrix(runif(300), ncol=10, dimnames=list(1:30, LETTERS[1:10]))) # Generates a sample data set. 
boxplot(y, col=1:length(y)) # Creates a box plot. A box plot (also known as a box-and-whisker diagram [https://en.wikipedia.org/wiki/File:Boxplot_vs_PDF.png]) is a graphical representation of a five-number summary, which consists of the smallest observation, lower quartile, median, upper quartile and largest observation.

# ROC Curves
# A variety of libraries are available for these tasks: ROCR [http://rocr.bioinf.mpi-sb.mpg.de/], ROC [http://bioconductor.org/packages/release/bioc/html/ROC.html].

## Miscellaeneous Plotting Utilities
plot(1:10); lines(1:10, lty=1) # Superimposes lines onto existing plots. The usage of plotted values will connect the data points. Different line types can be specified with argument "lty = 2". More on this can be found in the documentation for 'par'.
plot(x <- 1:10, y <- 1:10); abline(-1,1, col="green"); abline(1,1, col="red"); abline(v=5, col="blue"); abline(h=5, col="brown") # Function 'abline' adds lines in different colors to x-y-plot.

# Arranging Plots
x11(height=6, width=12, pointsize=12) # Defines the size of the plotting device. 
par(mfrow=c(2,3)); for(i in 1:6) { plot(1:10) } # With the argument 'par(mfrow=c(nrow,ncol))' one can define how several plots are arranged next to each other on the same plotting device.
nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), c(3,7), c(5,5), respect=TRUE); layout.show(nf); for(i in 1:3) { barplot(1:10) } # The 'layout()' function allows to devide the plotting device into variable numbers of rows and columns with the column-widths and the row-heights specified in the respective arguments.
split.screen(c(1,1)); plot(density(rnorm(10)), xlim=c(-2,2), ylim=c(0,1), col="red"); screen(1, new=FALSE); plot(density(rnorm(10)), xlim=c(-2,2), ylim=c(0,1), col="green", xaxt="n", yaxt="n", ylab="", xlab="", main="",bty="n"); close.screen(all = TRUE) # Possibility to overlay independent plots with 'split.screen' function. In ths example two density plots are plotted at the same scale ('xlim/ylim') and overlayed on the same graphics device using 'screen(1,new=FALSE))'. Split.screen is also very useful for generating multiple plots next to each other on a single device by setting the layout setting to several screens, e.g.: split.screen(c(2,2))

# Customize X-Y Axes
op <- par(mar=c(6,6,6,8)); barplot(1:10, names.arg=letters[1:10], cex.names=1.5, col=3, yaxt = "n") # Creates bar plot without y-axis and wider margins around plot. 
axis(2, las = 1, at=seq(0,10,0.5), lwd=2, cex.axis=1.5); axis(3, las = 1, lwd=2, cex.axis=1.5); axis(4, las = 1, at=seq(0,10,2), labels=month.name[1:6], lwd=2, cex.axis=1.5); par(op) # Adds three custom axes.
curve(sin(x), -pi, pi, axes=F, ylab="", xlab=""); axis(1, pos=0); axis(2, pos=0) # Plots a function in Cartesian-style coordinate system.

# Color Selection Utilities (see also: Earl Glynn's Color Chart, http://research.stowers-institute.org/efg/R/Color/Chart/)
y <- as.data.frame(matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))) # Generates a sample data set.
cat("palette():\n"); palette(); palette(rainbow(20, start=0.1, end=0.2)); cat("palette(rainbow(20, start=0.1, end=0.2)):\n"); palette(); palette("default") # Prints and modifies the color palette which is used when the argument 'col=' has a numeric index. The last step sets the palette back to its default setting.
par(mfrow=c(1,2)); barplot(as.matrix(y[1:8,]), beside=T, col=1:8); palette(rainbow(8, start=0.1, end=0.8)); barplot(as.matrix(y[1:8,]), beside=T, col=1:8); palette("default") # Example to call the default color palette with a numeric vector in the 'col=' argument'. In the second plot a modified palette is called the same way.
example(rainbow) # The function rainbow() allows to select various predefined color schemes.
barplot(as.matrix(y[1:10,]), beside=T, col=rainbow(10, start=0.1, end=0.2)) # Example to select 10 colors with the rainbow() function. The start and end values need to be between 0 and 1. The wider their distance the more diverse are the resulting colors.
barplot(as.matrix(y[1:10,]), beside=T, col=gray(c(1:10)/10)) # The gray() function allows to select any type of gray shades by providing values from 0 to 1.
colors()[1:10]; col2rgb(colors()[1:10]) # The colors() function allows to select colors by their name. The col2rgb() can translates them into the RGB color code.
library(RColorBrewer); example(brewer.pal) # Shows how to select color schemes with the RColorBrewer library.

# Adding Text
y <- as.data.frame(matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))) # Generates a sample data set.
plot(y[,1], y[,2]); text(0.2,0.2,"My Text") # Adds text at defined coordinates of plot.
plot(y[,1], y[,2]); mtext("My Text", 4) # Adds text to one of the four margin areas (side: 1-4) of the plot.
plot(y[,1], y[,2], pch=20, col="red", main="Plot of Symbols and Labels"); text(y[,1]+0.03, y[,2], rownames(y)) # Plots a scatter plot with the corresponding labels (here row names) next to the data points.

# Interactive Functions
locator() # Allows user to click on existing graph with left mouse button. System returns the corresponding x-y-coordinates after clicking on right mouse button.
text(locator(1), "My_Outlier", adj=0) # Adds text "My_Outlier" to graph at position where left mouse button is clicked.

# Saving Graphics to Files
jpeg("test.jpeg"); plot(1:10, 1:10); dev.off() # After the 'jpeg("test.jpeg")' command all graphs are redirected to the file "test.jpeg" in JPEG format. To export images with the highest quality, the default setting "quality = 75" needs to be changed to 100%. The actual image data are not written to the file until the 'dev.off()' command is executed!
pdf("test.pdf"); plot(1:10, 1:10); dev.off() # Same as above, but for pdf format. The pdf and svg formats provide often the best image quality, since they scale to any size without pixelation.
png("test.png"); plot(1:10, 1:10); dev.off() # Same as above, but for png format.
postscript("test.ps"); plot(1:10, 1:10); dev.off() # Same as above, but for PostScript format.
unlink("test.*")
