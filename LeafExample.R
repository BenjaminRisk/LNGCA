#------------------
# Benjamin Risk
# brisk@emory.edu
# Apply LCA to analysis of leaf feature data
# This code reproduces the analysis and figures in Section 7 of the LNGCA manuscript by Risk, Matteson, and Ruppert
#---------------------------------------------
# Leaf Data available at https://archive.ics.uci.edu/ml/machine-learning-databases/00288/
# 'Evaluation of Features for Leaf Discrimination', Pedro F.B. Silva, Andre R.S. Marcal, Rubim M. Almeida da Silva (2013). Springer Lecture Notes in Computer Science, Vol. 7950, 197-204.


library(ggplot2)
library(ProDenICA)
source('lcaFunctions.R')

set.seed(54321)
nRestart=10 # set equal to 10 to reproduce analysis. The complete analysis takes approximately 6.5 minutes with nRestart=10. 
# If restarts=0, the analysis takes 30 seconds, although one can see some differences owing to local optima.


leaf = read.csv('leaf.csv',header = FALSE) #downloaded from https://archive.ics.uci.edu/ml/machine-learning-databases/00288/
names(leaf)[1] = 'Species';
names(leaf)[2] = 'Sample No';
leaf$Species = as.character(leaf$Species)
names(leaf)[3:16] = c('Eccentricity','Aspect Ratio','Elongation','Solidity','Stochastic Convexity','Isoperimetric Factor','Max Indent Depth','Lobedness','Intensity','Contrast','Smoothness','Skewness','Uniformity','Entropy')
leafCont = leaf[,3:16]
scLeafCont = scale(leafCont)

# label species into three categories that will be useful for plotting later on:
indices = rep('All Others',nrow(scLeafCont));
indices[leaf$Species=='8'] = 'Linear 2';
indices[leaf$Species=='31' | leaf$Species=='34'] = 'Linear 1'

svd.leaf = svd(scLeafCont)
sc.ev.leaf = scale(svd.leaf$u)


## PCA:
dataPCA = data.frame('PC1' = sc.ev.leaf[,1], 'PC2' = sc.ev.leaf[,2],indices)
ggplot(dataPCA,aes(x=PC1, y=PC2, color=indices))+ geom_point(size=rel(5))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(3)),legend.position='none',plot.title=element_text(size=rel(3.5)))+xlab('PC1 Score') + ylab('PC2 Score')+ggtitle("PCA")
 
# PCA+Infomax:
est.PI2 = infomaxICA(X = sc.ev.leaf[,1:2],n.comp = 2,restarts=nRestart,whiten=FALSE)
estInfomax.S = order.likelihood(S = est.PI2$S,positive.skew = TRUE, distribution = 'logistic')
dataPI2 = data.frame('IC1'=estInfomax.S[,1],'IC2'=estInfomax.S[,2],'indices'=indices)

ggplot(dataPI2,aes(x=IC1, y=IC2, color=indices))+ geom_point(size=rel(3))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(2)),legend.position='none',plot.title=element_text(size=rel(2.5)))+xlab('IC1 Score')+ylab('IC2 Score')+ggtitle(expression(paste("PCA+Infomax: Q*=2",sep='')))

est.PI5 = infomaxICA(X = sc.ev.leaf[,1:5],n.comp=5,restarts=nRestart,whiten=FALSE)
estInfomax.S5 = order.likelihood(S = est.PI5$S,positive.skew=TRUE,distribution='logistic')
dataPI5 = data.frame('IC1'=estInfomax.S5[,1],'IC2'=estInfomax.S5[,2],indices)

ggplot(dataPI5,aes(x=IC1, y=IC2, color=indices))+ geom_point(size=rel(3))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(2)),legend.position='none',plot.title=element_text(size=rel(2.5)))+xlab('IC1 Score')+ylab('IC2 Score')+ggtitle(expression(paste("PCA+Infomax: Q*=5",sep='')))

## PCA+ProDenICA:
est.PP2 = mProDenICA(X = sc.ev.leaf[,1:2],restarts = nRestart,whiten = FALSE,df=8)


estProDenICA.S = order.likelihood(S = est.PP2$s,distribution='tiltedgaussian',df=8,positive.skew=TRUE)
dataPP2 = data.frame('IC1'=estProDenICA.S[,1],'IC2'=estProDenICA.S[,2],'indices'=indices)

ggplot(dataPP2,aes(x=IC1, y=IC2, color=indices))+ geom_point(size=rel(5))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(3)),legend.position='none',plot.title=element_text(size=rel(3.5)))+xlab('IC1 Score')+ylab('IC2 Score')+ggtitle("PCA+ProDenICA")

# PCA ProDenICA with 5 components, plotting first two
# this takes a few minutes:
est.PP5 = mProDenICA(X = sc.ev.leaf[,1:5],restarts = nRestart,whiten = FALSE,df=8)

estProDenICA.S5 = order.likelihood(S = est.PP5$s,distribution='tiltedgaussian',df=8,positive.skew=TRUE)
dataPP5 = data.frame('IC1'=estProDenICA.S5[,1],'IC2'=estProDenICA.S5[,2],'indices'=indices)
cor(estProDenICA.S5,estProDenICA.S)

ggplot(dataPP5,aes(x=IC1, y=IC2, color=indices))+ geom_point(size=rel(3))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(2)),legend.position='none',plot.title=element_text(size=rel(2.5)))+xlab('IC1 Score')+ylab('IC2 Score')+ggtitle(expression(paste("PCA+ProDenICA: Q*=5",sep='')))

#--------------------------------------------
## LogisLCA
est.logis.LC2 = mlcaFP(xData = scLeafCont, n.comp=2, restarts.pbyd=nRestart, restarts.dbyd=nRestart, distribution='logistic')

estLogisLC2.S = order.likelihood(S = est.logis.LC2$S,distribution='logistic',positive.skew=TRUE)
dataLL2 = data.frame(LC1 = estLogisLC2.S[,1], LC2=estLogisLC2.S[,2], indices)

ggplot(dataLL2,aes(x=LC1, y=LC2, color=indices))+ geom_point(size=rel(3))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(2)),legend.position='none',plot.title=element_text(size=rel(2.5)))+xlab('LC1 Score')+ylab('LC2 Score')+ggtitle(expression(paste("Logis-LCA: Q*=2",sep='')))

## Examine similarity of LCs when estimating 5 components:
est.logis.LC5 = mlcaFP(xData = scLeafCont, n.comp=5, restarts.pbyd=nRestart, restarts.dbyd=nRestart, distribution='logistic')
est.logis.LC5$loglik
estLogisLC5.S = order.likelihood(S = est.logis.LC5$S,distribution='logistic',positive.skew=TRUE)
cor(estLogisLC5.S,estLogisLC2.S)
dataLL5 = data.frame(LC1 = estLogisLC5.S[,1],LC2 = estLogisLC5.S[,2],indices)

ggplot(dataLL5,aes(x=LC1, y=LC2, color=indices))+ geom_point(size=rel(3))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(2)),legend.position='none',plot.title=element_text(size=rel(2.5)))+xlab('LC1 Score')+ylab('LC2 Score')+ggtitle(expression(paste("Logis-LCA: Q*=5",sep='')))


#-------------------------------------------
## SplineLCA

#takes approximately one minute to run
a= Sys.time()
est.spline.LC2 = mlcaFP(xData = scLeafCont, n.comp=2, restarts.pbyd=nRestart, restarts.dbyd=nRestart, distribution='tiltedgaussian',df=8)
Sys.time() - a
# Some of the initializations may reach Max iter, but that's okay because we use multiple initializations and the tolerances are relatively low


estSplineLC2.S = order.likelihood(S = est.spline.LC2$S, distribution='tiltedgaussian',df=8,positive.skew=TRUE)
dataSL = data.frame(LC1 = estSplineLC2.S[,1],LC2 = estSplineLC2.S[,2], 'indices' = indices)

ggplot(dataSL,aes(x=LC1, y=LC2, color=indices))+ geom_point(size=rel(5))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(3)),legend.position='none',plot.title=element_text(size=rel(3.5)))+xlab('LC1 Score')+ylab('LC2 Score')+ggtitle("Spline-LCA")

## Examine robustness to the number of components:
a = Sys.time()
est.spline.LC5 = mlcaFP(xData = scLeafCont, n.comp=5, restarts.pbyd=nRestart, restarts.dbyd=nRestart, distribution='tiltedgaussian', df=8)
Sys.time() - a

estSplineLC5.S = order.likelihood(est.spline.LC5$S,positive.skew = TRUE,distribution='tiltedgaussian',df=8)
cor(estSplineLC5.S,estSplineLC2.S)
dataSL5 = data.frame(LC1 = estSplineLC5.S[,1],LC2 = estSplineLC5.S[,2],indices)
ggplot(dataSL5,aes(x=LC1, y=LC2, color=indices))+ geom_point(size=rel(3))+theme(axis.text = element_text(size=rel(2)),axis.title=element_text(size=rel(2)),legend.position='none',plot.title=element_text(size=rel(2.5)))+xlab('LC1 Score')+ylab('LC2 Score')+ggtitle(expression(paste("Spline-LCA: Q* = 5",sep='')))


