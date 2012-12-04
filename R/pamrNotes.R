Part I. Prediction analysis of microarrays for R.

Original description of the method
Tibshirani R, Hastie T, Narasimhan B, Chu G. Diagnosis of multiple cancer
types by shrunken centroids of gene expression. Proc Natl Acad Sci U S A. 2002
May 14;99(10):6567-72. PubMed PMID: 12011421; PubMed Central PMCID: PMC124443.

We will be using the same GSM858 dataset, but slightly modified. First two rows contain sample IDs and infection status:

# GSM14498	GSM14499	GSM14500	GSM14501	GSM14513	GSM14514	GSM14515	GSM14516	GSM14506	GSM14507	GSM14508	GSM14502	GSM14503	GSM14504	GSM14505	GSM14509	GSM14510	GSM14511	GSM14512
# uninfected	uninfected	uninfected	uninfected	FRD875	FRD875	FRD875	FRD875	FRD1234	FRD1234	FRD1234	FRD1	FRD1	FRD1	FRD1	FRD440	FRD440	FRD440	FRD440

First two columns contain Affy IDs and gene names:

# 1007_s_at	DDR1
# 1053_at	RFC2
# 117_at	HSPA6
# 121_at	PAX8
# 1255_g_at	GUCA1A

The rest is a matrix of corresponding expression values. These should be log2 transformed and normalized, the numbers in the example are quantile normalized.
As a result we have 19 samples, classified as either "uninfected", or infected with FDR8789, FDR1234, FDR1, or FDR440. Let's read it in:

> response.data<-pamr.from.excel("gds858.q.pamr.txt",21, sample.labels=T)
> summary(response.data) # Get the idea what's inside

Note we have to specify the TOTAL number of columns, 21, that is 19 samples plus 2 columns for Affy IDs/gene names.

response.data$x is a matrix containing the expression profiles. It consists of 19 columns corresponding to samples and 22283 rows corresponding to genes. response.data$y contains the class labels of the 19 samples: uninfected, FDR8789, FDR1234, FDR1, or FDR440.

> response.train<-pamr.train(response.data) # Train the set based on provided class labels
> response.train # See what's inside

With this command you train PAM on the dataset. Typing response.train you get an output showind dependence of the number of nonsero (predictive) genes from the threshold, and corresponding misclassification error.
A more reliable error estimate than the number of misclassifications on the training set is  the 10-fold cross validation error:

> response.cv<-pamr.cv(response.train,response.data) # Cross-validate the classifier
> response.cv # Check how cross-validation performs

The output of this function looks very similar to the table above. The numbers in the last column are now the summed errors of all 10 cross validation steps. The CV error usually is bigger than the training error. Explain this gap.
The results of cross validation can be visualized by

> pamr.plotcv(response.cv) # Plot the cross-validated error curves

You will get two figures. In both, the x-axis represents different values of threshold (corresponding to different numbers of nonzero genes as shown on top of each figure) and the y-axis shows the number of misclassifications. The upper figure describes the whole dataset, the lower one describes each class individually. Explain the behaviour at the right
tail of the lower figure.
Using the results of cross validation, choose a threshold value t as a tradeoff between a small number of genes and a good generalization accuracy.

t<-2.5 # Select the threshold

In the next steps, vary t through a range of values and observe how the plots and figures change.

The function pamr.plotcen() plots the shrunken class centroids for each class, for genes surviving the threshold for at least one class.

> pamr.plotcen(response.train,response.data,t) # Plot the class centroids

This function prints a 5 x 5 confusion table.

> pamr.confusion(response.cv,t) # Compute the confusion matrix for a particular model

4 samples belong to FRD1 class. Only 1 is correctly classified, others are classified as belonging to tother classes. 4 samples in FRD440 class, on the other hand, are always correctly classified, resulting in 0 error rate. How many samples in "uninfected" class are correctly identified?
To get a visual impression of how clearly the two classes are separated by PAM, we plot the cross-validated sample probabilities:

> pamr.plotcvprob(response.train,response.data,t) # Plot the cross-validated class probabilities by class

The 19 samples (x-axis) are plotted against the probabilities to belong to different classes, color-coded. For each sample you see two small circles: the color of the class shows the probability that this sample belongs to the corresponding class and the other color indicating it belongs to the other class. A sample is put into that class for which probability exceeds 0.5.
For each gene surviving the threshold we get a figure showing the expression level of this gene over the whole set of samples by using the following command:

> pamr.geneplot(response.train,response.data,5) # Make a gene plot of the most significant genes

We used threshold = 5 to plot the smaller number of genes, otherwise there will be "Error in plot.new() : figure margins too large". We will use this threshold to get more information about top predictive genes.

> pamr.listgenes(response.train,response.data,5,genenames=T) # List the significant genes

Part II. Prediction using support vector machines (SVM)

Load e1071 package and prepare the data.

> install.packages("e1071") # install e1071 package, ?install packages for help
> library(e1071)
> data<-t(response.data$x) # Columns = genes, rows = samples
> labels<-response.data$y # Class labels

Don't get confused with the dimensions of the gene expression matrix. Usually columns are samples and rows are genes. But the SVM software
needs it just the other way round: the 19 rows stand for the samples and the 22283 columns for the genes.

Traninig and cross validation. We will begin with the simple linear kernel.

> svm.model<-svm(data,labels,type="C-classification",kernel="linear") # ?svm for help

Let's compute the training error.

> predicted <- predict(svm.model, data) # predict labels of training data
> sum(predicted != labels) # count differences
> table(true=labels, pred=predicted) # confusion matrix

The linear kernel separates the training set without errors! But how is its ability to predict unseen samples? We investigate by 10-fold cross validation

> svm.cross <- svm(data,labels, type="C-classification",kernel="linear", cross=10)

Let's check what results do we have, total accuracy, and accuracies during cross validation

> names(svm.cross)
> svm.cross$tot.accuracy
> svm.cross$accuracies

Try different kernel functions (polynomial, radial) with several parameters (degree, gamma) and values of the error weight (cost). What is the best result you can get?

Randomized experiment is always necessary. Now we will assign the samples randomly into classes and investigate how SVM will react.

# Randomized test
> labels.rand<-sample(response.data$y,19,replace=F)
> svm.rand <- svm(data, labels.rand, type="C-classification",kernel="linear)" #"polynomial", degree="2")
> predicted.rand<-predict(svm.rand,data) # predict labels of training data
> sum(predicted.rand != labels) # count differences
> table(true=labels, pred=predicted.rand) # confusion matrix

Now train the SVM again with cross validation! Compare the CV error on the biological and on the randomized data.

> svm.cross.rand <- svm(data,labels.rand, type="C-classification",kernel="linear", cross=10)
> svm.cross.rand$tot.accuracy
> svm.cross.rand$accuracies

How to select the most informative genes. For simplicity, we compare uninfected vs. FRD440 classes. 

> selected<-grep("uninfected|FRD440",response.data$y)
> data<-t(response.data$x[,selected])
> labels<-response.data$y[selected]
> table(labels) # How many samples in each class

We can select the most variable genes, which would be expected to reflect the differences between classes, and test SVM on them. But the right way of doing that is to use a t-statistic to select the genes with the most impact on classification. See the formula (1) image. Legend:
n : number of samples (= 8)
n+ : number of positive samples (= 4)
n- : number of negative samples (= 4)
mu+ : mean of expression levels in positive samples
mu- : mean of expression levels in negative samples
sigma2+ : variance in positive class
sigma2- : variance in negative class

Since R is good at matrix computations, you don't have to compute the t-value for each gene individually. It's more efficient to do it like this:

> FRD440<-data[labels == "FRD440",] # Extract the data for "positive" class
> mu.pos<-apply(FRD440,2,mean) # compute mean of each column (gene)
> var.pos<-apply(FRD440,2,var) # compute variance of each column (gene)

mu.pos and var.pos are vectors with 22283 entries. The entry of mu.pos at position i is the mean expression level of gene i over all samples with a positive label. Analogous for var.pos. Repeat the commands above on the negative samples and get vectors mu.neg and var.neg.
Use these building blocks to compute formula (1). This will give you a vector of length 22283: each entry is the t-score of the corresponding gene. Call this vector tscores. Now we order the vector tscores and select the 100 genes with highest t-score. Check the top 30 genes - do they look familiar?

> tscores<-(abs(mu.pos-mu.neg)/sqrt((4-1)*var.pos+(4-1)*var.neg))*sqrt(4*4*(8-2)/8)
> index <- order(tscores, decreasing=TRUE)
> data.sel <- data[,index[1:100]]
> response.data$genenames[index[1:30]]  # Let's look at the top genes

The matrix data.sel contains 19 rows (samples) and 100 columns (the selected genes). Train a SVM on this reduced dataset with different kernels and parameters and compare the results to those obtained with all genes. How do you explain the differences?
Vary the number of selected genes. How many genes do you need to still get a reasonable CV error?

