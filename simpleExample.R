# Benjamin Risk
# brisk@emory.edu
# Simple example demonstrating use of mlcaFP, which is the main function for
# estimating likelihood components. 
# For a detailed example see LeafExample.R
library(ProDenICA)
source('~/Dropbox/DimensionReduction/Programs/RcodeWebSupplement/lcaFunctions.R')

#------------------------
# Example 1:  t-distribution with df = 3
# In this example, LogisLCA and SplineLCA perform similarly
simData = SimLCA(n.samples = 1000, distribution='t',snr = 0.2, nu.vector=c(5,5)) #simulates LCA model with 2 LCs


# eigenvec for pre-whitening: uses left eigenvectors from the svd of X.center
# sqrtprec than using the square root of the precision matrix.
est.LogisLCA = mlcaFP(xData=simData$sim.X, n.comp=2, distribution='logistic',restarts.pbyd=5,whiten='eigenvec',maxit = 2000)
frobICA(S1 = est.LogisLCA$S, S2 = simData$sim.S, standardize=TRUE)

estM = est.M.ols(est.LogisLCA$S,simData$sim.X)
frobICA(M1=estM,M2=simData$sim.Ms,standardize=TRUE)


set.seed(123)
est.SplineLCA = mlcaFP(xData = simData$sim.X, n.comp = 2, distribution='tiltedgaussian', df=8,restarts.pbyd=5,whiten='eigenvec')
frobICA(S1 = est.SplineLCA$S, S2 = simData$sim.S,standardize=TRUE)
frobICA(M1=est.M.ols(est.SplineLCA$S,simData$sim.X),M2=simData$sim.Ms,standardize=TRUE)




#---------------------------------------------
#---------------------------------------------------
# Example 2: sub-gaussian mixture of normals
# you can play with the snr and note that PCA+ProDenICA does pretty well for large snr,
# but does poorly for low snr; Spline-LCA does well for both.
# LogisLCA does poorly for sub-gaussians

simData = SimTwoNorms(n.samples = 1000, distribution='mix-sub',snr = 0.2)
hist(simData$sim.S[,1],breaks=20)

est.Infomax = infomaxICA(X = simData$sim.X, n.comp=2, whiten = TRUE,maxit=2000,restarts = 5)
frobICA(S1 = est.Infomax$S, S2 = simData$sim.S, standardize=TRUE) #note the PMSE in the paper is actually the square of frobICA()

est.LogisLCA = mlcaFP(xData=simData$sim.X, n.comp=2, distribution='logistic',restarts.pbyd=5,whiten='eigenvec',maxit = 2000)
frobICA(S1 = est.LogisLCA$S, S2 = simData$sim.S, standardize=TRUE)

est.ProDenICA = mProDenICA(X = simData$sim.X, n.comp=2, whiten=TRUE, G = 'GPois', df=8,restarts = 10)
frobICA(S1 = est.ProDenICA$s, S2 = simData$sim.S,standardize=TRUE)

est.SplineLCA = mlcaFP(xData = simData$sim.X, n.comp = 2, distribution='tiltedgaussian',whiten='eigenvec',df=8,restarts.pbyd=10)
frobICA(S1 = est.SplineLCA$S, S2 = simData$sim.S,standardize=TRUE)

