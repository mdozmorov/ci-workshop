library(pamr)
response.data<-pamr.from.excel("gds858.q.pamr.txt",21, sample.labels=T)
summary(response.data) # Get the idea what's inside
pamr.menu(response.data) # All automatic
response.train<-pamr.train(response.data) # Train the set based on provided class labels
response.train # See what's inside
response.cv<-pamr.cv(response.train,response.data) # Cross-validate the classifier
response.cv # Check how cross-validation performs
pamr.plotcv(response.cv) # Plot the cross-validated error curves
t<-2.5 # Select the threshold
pamr.plotcen(response.train,response.data,t) # Plot the class centroids
pamr.confusion(response.cv,t) # Compute the confusion matrix for a particular model
pamr.plotcvprob(response.train,response.data,t) # Plot the cross-validated class probabilities by class
pamr.geneplot(response.train,response.data,5) # Make a gene plot of the most significant genes
pamr.listgenes(response.train,response.data,5,genenames=T) # List the significant genes


install.packages("e1071") # install e1071 package, ?install packages for help
library(e1071)
data<-t(response.data$x) # Columns = genes, rows = samples
labels<-response.data$y # Class labels
# Training and cross validation
svm.model<-svm(data,labels,type="C-classification",kernel="linear") # ?svm for help
predicted <- predict(svm.model, data) # predict labels of training data
sum(predicted != labels) # count differences
table(true=labels, pred=predicted) # confusion matrix
svm.cross <- svm(data,labels, type="C-classification",kernel="linear", cross=10)
names(svm.cross)
svm.cross$tot.accuracy
svm.cross$accuracies
# Randomized test
labels.rand<-sample(response.data$y,19,replace=F)
svm.rand <- svm(data, labels.rand, type="C-classification",kernel="linear)" #"polynomial", degree="2")
predicted.rand<-predict(svm.rand,data) # predict labels of training data
sum(predicted.rand != labels) # count differences
table(true=labels, pred=predicted.rand) # confusion matrix
svm.cross.rand <- svm(data,labels.rand, type="C-classification",kernel="linear", cross=10)
svm.cross.rand$tot.accuracy
svm.cross.rand$accuracies
# How to select the most informative genes
selected<-grep("uninfected|FRD440",response.data$y)
data<-t(response.data$x[,selected])
labels<-response.data$y[selected]
table(labels) # How many samples in each class
FRD440<-data[labels == "FRD440",] # Extract the data for "positive" class
mu.pos<-apply(FRD440,2,mean) # compute mean of each column (gene)
var.pos<-apply(FRD440,2,var) # compute variance of each column (gene)
uninfected<-data[labels == "uninfected",] # Extract the data for "negative" class
mu.neg<-apply(uninfected,2,mean) # compute mean of each column (gene)               
var.neg<-apply(uninfected,2,var) # compute variance of each column (gene)
tscores<-(abs(mu.pos-mu.neg)/sqrt((4-1)*var.pos+(4-1)*var.neg))*sqrt(4*4*(8-2)/8)
index <- order(tscores, decreasing=TRUE)
data.sel <- data[,index[1:100]]
response.data$genenames[index[1:30]]  # Let's look at the top genes               
