###################################################################
# R codes for simulation setting 1 in single-subject simulations
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(irlba)
library(R.matlab)

source("00_utils.R")

###################################################################
# General simulations settings
###################################################################

# The three given SNR levels
my_SNR = c(0.4,1.5,3)
# Number of replications
rep = 100
# Variance level of inactive points in components
var_level = 0
# Overall seed
seed = 2023
# Number of time points in each simulation
my_nTR = 50
# AR structure in each simulation
my_phi = 0.47
# The convergence criterion
my_eps = 1e-6
# The maximum number of iterations
my_maxit = 1000
# Whether it's a noisy ICA model
my_noisyICA = TRUE

###################################################################
# Simulations - Low SNR
###################################################################

S_PMSE1 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))
M_PMSE1 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))
TIME1 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))

load("../../Data/sim123_sparse_lowSNR.RData")
for (i in 1:rep) {
  xmat = xmat_list[[i]]
  smat = smat_list[[i]]
  mmat = mmat_list[[i]]

  # standardize xmat
  xmat_center = scale(xmat,center = T,scale = F)

  ######################################
  # Sparse ICA - single restart
  ######################################
  select_sparseICA= BIC_sparseICA_Rcpp(xData = xmat,n.comp = 3,U.list = NULL,whiten = "eigenvec",
                                       eps = my_eps,maxit = my_maxit,method = "C",irlba = TRUE,
                                       verbose = FALSE,BIC_plot = FALSE,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu

  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat,n.comp = 3,nu = my_nu,whiten = "eigenvec",method = "C",
                           restarts = 1,irlba = TRUE,eps = my_eps,maxit = my_maxit,verbose = F)
  sparse.end.time = Sys.time()

  #my_S_sparse = matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE1$sparse[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE1$sparse[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  TIME1$sparse[i] = as.numeric(sparse.end.time - sparse.start.time)

  ######################################
  # Sparse ICA - 40 restarts
  ######################################
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat,n.comp = 3,nu = my_nu,whiten = "eigenvec",method = "C",
                           restarts = 40,irlba = TRUE,eps = my_eps,maxit = my_maxit,verbose = F)
  sparse.end.time = Sys.time()

  #my_S_sparse = matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE1$sparse40[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE1$sparse40[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  TIME1$sparse40[i] = as.numeric(sparse.end.time - sparse.start.time)

  ######################################
  # Infomax ICA - single restart
  ######################################
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit)
  infomax.end.time = Sys.time()

  #my_S_infomax = matchICA(my_infomaxICA$S,smat)
  #my_M_infomax = est.M.ols(my_S_infomax,xmat_center)
  S_PMSE1$infomax[i] = sqrt(frobICA(S1=my_infomaxICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$infomax[i] = sqrt(frobICA(M1=my_infomaxICA$M,M2=mmat,standardize = TRUE))
  TIME1$infomax[i] = as.numeric(infomax.end.time - infomax.start.time)

  ######################################
  # Infomax ICA - 40 restarts
  ######################################
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 40)
  infomax.end.time = Sys.time()

  #my_S_infomax = matchICA(my_infomaxICA$S,smat)
  #my_M_infomax = est.M.ols(my_S_infomax,xmat_center)
  S_PMSE1$infomax40[i] = sqrt(frobICA(S1=my_infomaxICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$infomax40[i] = sqrt(frobICA(M1=my_infomaxICA$M,M2=mmat,standardize = TRUE))
  TIME1$infomax40[i] = as.numeric(infomax.end.time - infomax.start.time)

  ######################################
  # Fast ICA - single restart
  ######################################
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method = "C")
  fast.end.time = Sys.time()

  #my_S_fast = matchICA(my_fastICA$S,smat)
  #my_M_fast = est.M.ols(my_S_fast,xmat_center)
  S_PMSE1$fast[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$fast[i] = sqrt(frobICA(M1=my_fastICA$A,M2=mmat,standardize = TRUE))
  TIME1$fast[i] = as.numeric(fast.end.time - fast.start.time)

  ######################################
  # Fast ICA - 40 restarts
  ######################################
  fast.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,restarts = 40,method = "C")
  fast.end.time = Sys.time()

  #my_S_fast = matchICA(my_fastICA$S,smat)
  #my_M_fast = est.M.ols(my_S_fast,xmat_center)
  S_PMSE1$fast40[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE1$fast40[i] = sqrt(frobICA(M1=my_fastICA$A,M2=mmat,standardize = TRUE))
  TIME1$fast40[i] = as.numeric(fast.end.time - fast.start.time)

  ######################################
  # Sparse ICA - EBM
  ######################################
  filename = paste0("../../Results/EBM/low/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(rnorm(47*1089),nrow = 1089))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  S_PMSE1$ebm[i] = sqrt(frobICA(S1=my_S_EBM,S2=smat,standardize = TRUE))
  M_PMSE1$ebm[i] = sqrt(frobICA(M1=my_M_EBM,M2=mmat,standardize = TRUE))
  TIME1$ebm[i] = as.numeric(dat$tEnd)

  ######################################
  # SparseFastICA
  ######################################
  filename = paste0("../../Results/sparsefast/low/estS_",i,".mat")
  dat = readMat(filename)
  #my_S_sfi = matchICA(t(dat$myS),smat)
  #my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_PMSE1$sfi[i] = sqrt(frobICA(S1=t(dat$myS),S2=smat,standardize = TRUE))
  M_PMSE1$sfi[i] = sqrt(frobICA(M1=t(dat$A),M2=mmat,standardize = TRUE))
  TIME1$sfi[i] = as.numeric(dat$tEnd)
}

###################################################################
# Simulations - Medium SNR
###################################################################
S_PMSE2 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))
M_PMSE2 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))
TIME2 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))

load("../../Data/sim123_sparse_mediumSNR.RData")
for (i in 1:rep) {
  xmat = xmat_list[[i]]
  smat = smat_list[[i]]
  mmat = mmat_list[[i]]
  
  # standardize xmat
  xmat_center = scale(xmat,center = T,scale = F)
  
  ######################################
  # Sparse ICA - single restart
  ######################################
  select_sparseICA= BIC_sparseICA_Rcpp(xData = xmat,n.comp = 3,U.list = NULL,whiten = "eigenvec",
                                       eps = my_eps,maxit = my_maxit,method = "C",irlba = TRUE,
                                       verbose = FALSE,BIC_plot = FALSE,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat,n.comp = 3,nu = my_nu,whiten = "eigenvec",method = "C",
                           restarts = 1,irlba = TRUE,eps = my_eps,maxit = my_maxit,verbose = F)
  sparse.end.time = Sys.time()
  
  #my_S_sparse = matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE2$sparse[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE2$sparse[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  TIME2$sparse[i] = as.numeric(sparse.end.time - sparse.start.time)
  
  ######################################
  # Sparse ICA - 40 restarts
  ######################################
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat,n.comp = 3,nu = my_nu,whiten = "eigenvec",method = "C",
                           restarts = 40,irlba = TRUE,eps = my_eps,maxit = my_maxit,verbose = F)
  sparse.end.time = Sys.time()
  
  #my_S_sparse = matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE2$sparse40[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE2$sparse40[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  TIME2$sparse40[i] = as.numeric(sparse.end.time - sparse.start.time)
  
  ######################################
  # Infomax ICA - single restart
  ######################################
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit)
  infomax.end.time = Sys.time()
  
  #my_S_infomax = matchICA(my_infomaxICA$S,smat)
  #my_M_infomax = est.M.ols(my_S_infomax,xmat_center)
  S_PMSE2$infomax[i] = sqrt(frobICA(S1=my_infomaxICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$infomax[i] = sqrt(frobICA(M1=my_infomaxICA$M,M2=mmat,standardize = TRUE))
  TIME2$infomax[i] = as.numeric(infomax.end.time - infomax.start.time)
  
  ######################################
  # Infomax ICA - 40 restarts
  ######################################
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 40)
  infomax.end.time = Sys.time()
  
  #my_S_infomax = matchICA(my_infomaxICA$S,smat)
  #my_M_infomax = est.M.ols(my_S_infomax,xmat_center)
  S_PMSE2$infomax40[i] = sqrt(frobICA(S1=my_infomaxICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$infomax40[i] = sqrt(frobICA(M1=my_infomaxICA$M,M2=mmat,standardize = TRUE))
  TIME2$infomax40[i] = as.numeric(infomax.end.time - infomax.start.time)
  
  ######################################
  # Fast ICA - single restart
  ######################################
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method = "C")
  fast.end.time = Sys.time()
  
  #my_S_fast = matchICA(my_fastICA$S,smat)
  #my_M_fast = est.M.ols(my_S_fast,xmat_center)
  S_PMSE2$fast[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$fast[i] = sqrt(frobICA(M1=my_fastICA$A,M2=mmat,standardize = TRUE))
  TIME2$fast[i] = as.numeric(fast.end.time - fast.start.time)
  
  ######################################
  # Fast ICA - 40 restarts
  ######################################
  fast.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,restarts = 40,method = "C")
  fast.end.time = Sys.time()
  
  #my_S_fast = matchICA(my_fastICA$S,smat)
  #my_M_fast = est.M.ols(my_S_fast,xmat_center)
  S_PMSE2$fast40[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE2$fast40[i] = sqrt(frobICA(M1=my_fastICA$A,M2=mmat,standardize = TRUE))
  TIME2$fast40[i] = as.numeric(fast.end.time - fast.start.time)
  
  ######################################
  # Sparse ICA - EBM
  ######################################
  filename = paste0("../../Results/EBM/medium/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(rnorm(47*1089),nrow = 1089))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  S_PMSE2$ebm[i] = sqrt(frobICA(S1=my_S_EBM,S2=smat,standardize = TRUE))
  M_PMSE2$ebm[i] = sqrt(frobICA(M1=my_M_EBM,M2=mmat,standardize = TRUE))
  TIME2$ebm[i] = as.numeric(dat$tEnd)
  
  ######################################
  # SparseFastICA
  ######################################
  filename = paste0("../../Results/sparsefast/medium/estS_",i,".mat")
  dat = readMat(filename)
  #my_S_sfi = matchICA(t(dat$myS),smat)
  #my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_PMSE2$sfi[i] = sqrt(frobICA(S1=t(dat$myS),S2=smat,standardize = TRUE))
  M_PMSE2$sfi[i] = sqrt(frobICA(M1=t(dat$A),M2=mmat,standardize = TRUE))
  TIME2$sfi[i] = as.numeric(dat$tEnd)
}


###################################################################
# Simulations - High SNR
###################################################################
S_PMSE3 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))
M_PMSE3 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))
TIME3 = data.frame(sparse=numeric(rep),fast=numeric(rep),infomax=numeric(rep),ebm=numeric(rep),sfi=numeric(rep),sparse40=numeric(rep),fast40=numeric(rep),infomax40=numeric(rep))

load("../../Data/sim123_sparse_highSNR.RData")
for (i in 1:rep) {
  xmat = xmat_list[[i]]
  smat = smat_list[[i]]
  mmat = mmat_list[[i]]
  
  # standardize xmat
  xmat_center = scale(xmat,center = T,scale = F)
  
  ######################################
  # Sparse ICA - single restart
  ######################################
  select_sparseICA= BIC_sparseICA_Rcpp(xData = xmat,n.comp = 3,U.list = NULL,whiten = "eigenvec",
                                       eps = my_eps,maxit = my_maxit,method = "C",irlba = TRUE,
                                       verbose = FALSE,BIC_plot = FALSE,nu_list = seq(0.1,4,0.1))
  my_nu = select_sparseICA$best_nu
  
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat,n.comp = 3,nu = my_nu,whiten = "eigenvec",method = "C",
                           restarts = 1,irlba = TRUE,eps = my_eps,maxit = my_maxit,verbose = F)
  sparse.end.time = Sys.time()
  
  #my_S_sparse = matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE3$sparse[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE3$sparse[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  TIME3$sparse[i] = as.numeric(sparse.end.time - sparse.start.time)
  
  ######################################
  # Sparse ICA - 40 restarts
  ######################################
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat,n.comp = 3,nu = my_nu,whiten = "eigenvec",method = "C",
                           restarts = 40,irlba = TRUE,eps = my_eps,maxit = my_maxit,verbose = F)
  sparse.end.time = Sys.time()
  
  #my_S_sparse = matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(my_sparseICA$estS,xmat_center)
  S_PMSE3$sparse40[i] = sqrt(frobICA(S1=my_sparseICA$estS,S2=smat,standardize = TRUE))
  M_PMSE3$sparse40[i] = sqrt(frobICA(M1=my_M_sparse,M2=mmat,standardize = TRUE))
  TIME3$sparse40[i] = as.numeric(sparse.end.time - sparse.start.time)
  
  ######################################
  # Infomax ICA - single restart
  ######################################
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit)
  infomax.end.time = Sys.time()
  
  #my_S_infomax = matchICA(my_infomaxICA$S,smat)
  #my_M_infomax = est.M.ols(my_S_infomax,xmat_center)
  S_PMSE3$infomax[i] = sqrt(frobICA(S1=my_infomaxICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$infomax[i] = sqrt(frobICA(M1=my_infomaxICA$M,M2=mmat,standardize = TRUE))
  TIME3$infomax[i] = as.numeric(infomax.end.time - infomax.start.time)
  
  ######################################
  # Infomax ICA - 40 restarts
  ######################################
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 40)
  infomax.end.time = Sys.time()
  
  #my_S_infomax = matchICA(my_infomaxICA$S,smat)
  #my_M_infomax = est.M.ols(my_S_infomax,xmat_center)
  S_PMSE3$infomax40[i] = sqrt(frobICA(S1=my_infomaxICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$infomax40[i] = sqrt(frobICA(M1=my_infomaxICA$M,M2=mmat,standardize = TRUE))
  TIME3$infomax40[i] = as.numeric(infomax.end.time - infomax.start.time)
  
  ######################################
  # Fast ICA - single restart
  ######################################
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method = "C")
  fast.end.time = Sys.time()
  
  #my_S_fast = matchICA(my_fastICA$S,smat)
  #my_M_fast = est.M.ols(my_S_fast,xmat_center)
  S_PMSE3$fast[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$fast[i] = sqrt(frobICA(M1=my_fastICA$A,M2=mmat,standardize = TRUE))
  TIME3$fast[i] = as.numeric(fast.end.time - fast.start.time)
  
  ######################################
  # Fast ICA - 40 restarts
  ######################################
  fast.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,restarts = 40,method = "C")
  fast.end.time = Sys.time()
  
  #my_S_fast = matchICA(my_fastICA$S,smat)
  #my_M_fast = est.M.ols(my_S_fast,xmat_center)
  S_PMSE3$fast40[i] = sqrt(frobICA(S1=my_fastICA$S,S2=smat,standardize = TRUE))
  M_PMSE3$fast40[i] = sqrt(frobICA(M1=my_fastICA$A,M2=mmat,standardize = TRUE))
  TIME3$fast40[i] = as.numeric(fast.end.time - fast.start.time)
  
  ######################################
  # Sparse ICA - EBM
  ######################################
  filename = paste0("../../Results/EBM/high/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(rnorm(47*1089),nrow = 1089))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  S_PMSE3$ebm[i] = sqrt(frobICA(S1=my_S_EBM,S2=smat,standardize = TRUE))
  M_PMSE3$ebm[i] = sqrt(frobICA(M1=my_M_EBM,M2=mmat,standardize = TRUE))
  TIME3$ebm[i] = as.numeric(dat$tEnd)
  
  ######################################
  # SparseFastICA
  ######################################
  filename = paste0("../../Results/sparsefast/high/estS_",i,".mat")
  dat = readMat(filename)
  #my_S_sfi = matchICA(t(dat$myS),smat)
  #my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_PMSE3$sfi[i] = sqrt(frobICA(S1=t(dat$myS),S2=smat,standardize = TRUE))
  M_PMSE3$sfi[i] = sqrt(frobICA(M1=t(dat$A),M2=mmat,standardize = TRUE))
  TIME3$sfi[i] = as.numeric(dat$tEnd)
}

save(S_PMSE1,S_PMSE2,S_PMSE3,M_PMSE1,M_PMSE2,M_PMSE3,TIME1,TIME2,TIME3,file = "../../Results/02_PRMSE_sim123_sparse.RData")

