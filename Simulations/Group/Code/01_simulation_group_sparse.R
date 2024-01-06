###################################################################
# R codes for group-levels simulations
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
source("00_utils_group.R")

###################################################################
# General simulations settings
###################################################################

# Random seeds
seed1 = 79643
seed2 = 82569
seed3 = 49853
# the parameter for Gamma distribution
gamma.rate = 1e-4
gamma.shape = 0.02
# number of subjects
nsub = 20
# number of repetition
nrep = 50
# Full width at half maximum
FWHM = 6
# variance of inactive elements in components
var.inactive = 0
# whether it's under noisy ICA model
noisyICA=TRUE
# number of group components
n.group=3

###################################################################
# Group ICA in low SNR setting
###################################################################

SNR = 'low'
S_PRMSE1 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))
M_avg_PRMSE1 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))
TIME1 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))

set.seed(seed1)
for (k in 1:nrep) {
  # data generation 
  simData = list()
  for (i in 1:nsub){
    simData[[i]] = SimFMRI.ngca(snr = SNR, noisyICA=noisyICA,
                                gamma.rate = gamma.rate, gamma.shape = gamma.shape, 
                                FWHM = FWHM,var.inactive = var.inactive)
  }
  
  bic = BIC_group_sparse_ICA(simData,n.group = 3,BIC_plot = T)
  
  start.time = Sys.time()
  est_sparse = group_sparse_ICA(simData, n.group = n.group,restarts = 1,nu=bic$best_nu)
  end.time = Sys.time()
  TIME1$sparse[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_sparse40 = group_sparse_ICA(simData, n.group = n.group,restarts = 40,nu=bic$best_nu)
  end.time = Sys.time()
  TIME1$sparse40[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_fast = group_fast_ICA(simData, n.group = n.group,restarts = 1)
  end.time = Sys.time()
  TIME1$fast[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_fast40 = group_fast_ICA(simData, n.group = n.group,restarts = 40)
  end.time = Sys.time()
  TIME1$fast40[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_infomax = group_infomax_ICA(simData, n.group = n.group,restarts = 1)
  end.time = Sys.time()
  TIME1$infomax[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_infomax40 = group_infomax_ICA(simData, n.group = n.group,restarts = 40)
  end.time = Sys.time()
  TIME1$infomax40[k]=as.numeric(end.time - start.time)
  
  S_PRMSE1$sparse[k]=sqrt(frobICA(S1=est_sparse$estS,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE1$fast[k]=sqrt(frobICA(S1=est_fast$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE1$infomax[k]=sqrt(frobICA(S1=est_infomax$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE1$sparse40[k]=sqrt(frobICA(S1=est_sparse40$estS,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE1$fast40[k]=sqrt(frobICA(S1=est_fast40$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE1$infomax40[k]=sqrt(frobICA(S1=est_infomax40$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))

  a=0
  b=0
  c=0
  d=0
  e=0
  f=0
  for (j in 1:length(simData)) {
    xmat_center=scale(simData[[j]]$X,center = T,scale = F)
    M_sparse=est.M.ols(est_sparse$estS,xmat_center)
    M_fast=est.M.ols(est_fast$S,xmat_center)
    M_infomax=est.M.ols(est_infomax$S,xmat_center)
    M_sparse40=est.M.ols(est_sparse40$estS,xmat_center)
    M_fast40=est.M.ols(est_fast40$S,xmat_center)
    M_infomax40=est.M.ols(est_infomax40$S,xmat_center)
    
    a=a+sqrt(frobICA(M1=M_sparse,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    b=b+sqrt(frobICA(M1=M_fast,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    c=c+sqrt(frobICA(M1=M_infomax,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    d=d+sqrt(frobICA(M1=M_sparse40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    e=e+sqrt(frobICA(M1=M_fast40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    f=f+sqrt(frobICA(M1=M_infomax40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
  }
  M_avg_PRMSE1$sparse[k]=a/nsub
  M_avg_PRMSE1$fast[k]=b/nsub
  M_avg_PRMSE1$infomax[k]=c/nsub
  M_avg_PRMSE1$sparse40[k]=d/nsub
  M_avg_PRMSE1$fast40[k]=e/nsub
  M_avg_PRMSE1$infomax40[k]=f/nsub

  cat("The ",k,"th replication finished!\n")
}


###################################################################
# Group ICA in medium SNR setting
###################################################################

SNR = 'medium'
S_PRMSE2 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))
M_avg_PRMSE2 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))
TIME2 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))

set.seed(seed2)
for (k in 1:nrep) {
  # data generation 
  simData = list()
  for (i in 1:nsub){
    simData[[i]] = SimFMRI.ngca(snr = SNR, noisyICA=noisyICA,
                                gamma.rate = gamma.rate, gamma.shape = gamma.shape, 
                                FWHM = FWHM,var.inactive = var.inactive)
  }
  
  bic = BIC_group_sparse_ICA(simData,n.group = 3,BIC_plot = T)
  
  start.time = Sys.time()
  est_sparse = group_sparse_ICA(simData, n.group = n.group,restarts = 1,nu=bic$best_nu)
  end.time = Sys.time()
  TIME2$sparse[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_sparse40 = group_sparse_ICA(simData, n.group = n.group,restarts = 40,nu=bic$best_nu)
  end.time = Sys.time()
  TIME2$sparse40[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_fast = group_fast_ICA(simData, n.group = n.group,restarts = 1)
  end.time = Sys.time()
  TIME2$fast[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_fast40 = group_fast_ICA(simData, n.group = n.group,restarts = 40)
  end.time = Sys.time()
  TIME2$fast40[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_infomax = group_infomax_ICA(simData, n.group = n.group,restarts = 1)
  end.time = Sys.time()
  TIME2$infomax[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_infomax40 = group_infomax_ICA(simData, n.group = n.group,restarts = 40)
  end.time = Sys.time()
  TIME2$infomax40[k]=as.numeric(end.time - start.time)
  
  S_PRMSE2$sparse[k]=sqrt(frobICA(S1=est_sparse$estS,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE2$fast[k]=sqrt(frobICA(S1=est_fast$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE2$infomax[k]=sqrt(frobICA(S1=est_infomax$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE2$sparse40[k]=sqrt(frobICA(S1=est_sparse40$estS,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE2$fast40[k]=sqrt(frobICA(S1=est_fast40$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE2$infomax40[k]=sqrt(frobICA(S1=est_infomax40$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  
  a=0
  b=0
  c=0
  d=0
  e=0
  f=0
  for (j in 1:length(simData)) {
    xmat_center=scale(simData[[j]]$X,center = T,scale = F)
    M_sparse=est.M.ols(est_sparse$estS,xmat_center)
    M_fast=est.M.ols(est_fast$S,xmat_center)
    M_infomax=est.M.ols(est_infomax$S,xmat_center)
    M_sparse40=est.M.ols(est_sparse40$estS,xmat_center)
    M_fast40=est.M.ols(est_fast40$S,xmat_center)
    M_infomax40=est.M.ols(est_infomax40$S,xmat_center)
    
    a=a+sqrt(frobICA(M1=M_sparse,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    b=b+sqrt(frobICA(M1=M_fast,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    c=c+sqrt(frobICA(M1=M_infomax,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    d=d+sqrt(frobICA(M1=M_sparse40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    e=e+sqrt(frobICA(M1=M_fast40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    f=f+sqrt(frobICA(M1=M_infomax40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
  }
  M_avg_PRMSE2$sparse[k]=a/nsub
  M_avg_PRMSE2$fast[k]=b/nsub
  M_avg_PRMSE2$infomax[k]=c/nsub
  M_avg_PRMSE2$sparse40[k]=d/nsub
  M_avg_PRMSE2$fast40[k]=e/nsub
  M_avg_PRMSE2$infomax40[k]=f/nsub
  
  cat("The ",k,"th replication finished!\n")
}


###################################################################
# Group ICA in high SNR setting
###################################################################

SNR = 'high'
S_PRMSE3 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))
M_avg_PRMSE3 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))
TIME3 = data.frame(sparse=rep(NA,nrep),fast=rep(NA,nrep),infomax=rep(NA,nrep),sparse40=rep(NA,nrep),fast40=rep(NA,nrep),infomax40=rep(NA,nrep))

set.seed(seed3)
for (k in 1:nrep) {
  # data generation 
  simData = list()
  for (i in 1:nsub){
    simData[[i]] = SimFMRI.ngca(snr = SNR, noisyICA=noisyICA,
                                gamma.rate = gamma.rate, gamma.shape = gamma.shape, 
                                FWHM = FWHM,var.inactive = var.inactive)
  }
  
  bic = BIC_group_sparse_ICA(simData,n.group = 3,BIC_plot = T)
  
  start.time = Sys.time()
  est_sparse = group_sparse_ICA(simData, n.group = n.group,restarts = 1,nu=bic$best_nu)
  end.time = Sys.time()
  TIME3$sparse[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_sparse40 = group_sparse_ICA(simData, n.group = n.group,restarts = 40,nu=bic$best_nu)
  end.time = Sys.time()
  TIME3$sparse40[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_fast = group_fast_ICA(simData, n.group = n.group,restarts = 1)
  end.time = Sys.time()
  TIME3$fast[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_fast40 = group_fast_ICA(simData, n.group = n.group,restarts = 40)
  end.time = Sys.time()
  TIME3$fast40[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_infomax = group_infomax_ICA(simData, n.group = n.group,restarts = 1)
  end.time = Sys.time()
  TIME3$infomax[k]=as.numeric(end.time - start.time)
  
  start.time = Sys.time()
  est_infomax40 = group_infomax_ICA(simData, n.group = n.group,restarts = 40)
  end.time = Sys.time()
  TIME3$infomax40[k]=as.numeric(end.time - start.time)
  
  S_PRMSE3$sparse[k]=sqrt(frobICA(S1=est_sparse$estS,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE3$fast[k]=sqrt(frobICA(S1=est_fast$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE3$infomax[k]=sqrt(frobICA(S1=est_infomax$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE3$sparse40[k]=sqrt(frobICA(S1=est_sparse40$estS,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE3$fast40[k]=sqrt(frobICA(S1=est_fast40$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  S_PRMSE3$infomax40[k]=sqrt(frobICA(S1=est_infomax40$S,S2=simData[[1]]$S[,1:3],standardize = TRUE))
  
  a=0
  b=0
  c=0
  d=0
  e=0
  f=0
  for (j in 1:length(simData)) {
    xmat_center=scale(simData[[j]]$X,center = T,scale = F)
    M_sparse=est.M.ols(est_sparse$estS,xmat_center)
    M_fast=est.M.ols(est_fast$S,xmat_center)
    M_infomax=est.M.ols(est_infomax$S,xmat_center)
    M_sparse40=est.M.ols(est_sparse40$estS,xmat_center)
    M_fast40=est.M.ols(est_fast40$S,xmat_center)
    M_infomax40=est.M.ols(est_infomax40$S,xmat_center)
    
    a=a+sqrt(frobICA(M1=M_sparse,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    b=b+sqrt(frobICA(M1=M_fast,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    c=c+sqrt(frobICA(M1=M_infomax,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    d=d+sqrt(frobICA(M1=M_sparse40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    e=e+sqrt(frobICA(M1=M_fast40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
    f=f+sqrt(frobICA(M1=M_infomax40,M2=simData[[j]]$Ms[1:3,],standardize = TRUE))
  }
  M_avg_PRMSE3$sparse[k]=a/nsub
  M_avg_PRMSE3$fast[k]=b/nsub
  M_avg_PRMSE3$infomax[k]=c/nsub
  M_avg_PRMSE3$sparse40[k]=d/nsub
  M_avg_PRMSE3$fast40[k]=e/nsub
  M_avg_PRMSE3$infomax40[k]=f/nsub
  
  cat("The ",k,"th replication finished!\n")
}

save(S_PRMSE1,S_PRMSE2,S_PRMSE3,M_avg_PRMSE1,M_avg_PRMSE2,M_avg_PRMSE3,TIME1,TIME2,TIME3,file = "../Data/01_group_sparse.RData")

