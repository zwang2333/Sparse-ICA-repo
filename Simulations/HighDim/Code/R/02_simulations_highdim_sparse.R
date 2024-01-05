###################################################################
# R codes for simulation setting 2 (High-Dimensional) in single-subject simulations
# 02: run simulations
# Zihang Wang
# 1/5/2023
###################################################################
rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(R.matlab)
library(ciftiTools)
source("00_utils.R")

# ciftiTools package requires users using the following codes to point to the path to workbench
# For Mac
ciftiTools.setOption('wb_path', '/Applications/workbench')
# For Windows
#ciftiTools.setOption('wb_path', 'C:/Software/Workbench/workbench-windows64-v1.5.0/workbench') 
#ciftiTools.setOption('wb_path', 'D:/Softwares/workbench/workbench') 

###################################################################
# Load true components data
###################################################################

load("../../Data/group_SparseICA.RData")
load("../../Data/template_xifti.RData")
load("../../Data/template_timecourses2.RData")

###################################################################
# General simulation settings
###################################################################

# Simulate in low SNR level
my_SNR = 0.4
# Number of replications
rep = 100
# Level for variance of inactive points in components
var_level = 0
# Random seeds
seed1 = 159963
seed2 = 753741
seed3 = 852367
seed4 = 951753
seed5 = 364192
# number of time points in each simulation
my_nTR = 120
# AR structure in each simulation
my_phi = 0.47
# Convergence criterion
my_eps = 1e-6
# Maximum number of iterations
my_maxit = 1000
# Whether it's noisy ICA
my_noisyICA = TRUE
# group components to be used
my_components = c(17,6,23)
# the list of candidate nu
nu_list = seq(0.1,3,0.1)

###################################################################
# Start simulations
###################################################################

# Low SNR
S_cor1 = data.frame(sparse_IC1=numeric(rep),
                     sparse_IC2=numeric(rep),
                     sparse_IC3=numeric(rep),
                    sparse40_IC1=numeric(rep),
                    sparse40_IC2=numeric(rep),
                    sparse40_IC3=numeric(rep),
                     fast_IC1=numeric(rep),
                     fast_IC2=numeric(rep),
                     fast_IC3=numeric(rep),
                    fast40_IC1=numeric(rep),
                    fast40_IC2=numeric(rep),
                    fast40_IC3=numeric(rep),
                     infomax_IC1=numeric(rep),
                     infomax_IC2=numeric(rep),
                     infomax_IC3=numeric(rep),
                    infomax40_IC1=numeric(rep),
                    infomax40_IC2=numeric(rep),
                    infomax40_IC3=numeric(rep),
                    EBM_IC1=numeric(rep),
                    EBM_IC2=numeric(rep),
                    EBM_IC3=numeric(rep),
                    sfi_IC1=numeric(rep),
                    sfi_IC2=numeric(rep),
                    sfi_IC3=numeric(rep))
M_cor1 = data.frame(sparse_IC1=numeric(rep),
                     sparse_IC2=numeric(rep),
                     sparse_IC3=numeric(rep),
                    sparse40_IC1=numeric(rep),
                    sparse40_IC2=numeric(rep),
                    sparse40_IC3=numeric(rep),
                     fast_IC1=numeric(rep),
                     fast_IC2=numeric(rep),
                     fast_IC3=numeric(rep),
                    fast40_IC1=numeric(rep),
                    fast40_IC2=numeric(rep),
                    fast40_IC3=numeric(rep),
                     infomax_IC1=numeric(rep),
                     infomax_IC2=numeric(rep),
                     infomax_IC3=numeric(rep),
                    infomax40_IC1=numeric(rep),
                    infomax40_IC2=numeric(rep),
                    infomax40_IC3=numeric(rep),
                    EBM_IC1=numeric(rep),
                    EBM_IC2=numeric(rep),
                    EBM_IC3=numeric(rep),
                    sfi_IC1=numeric(rep),
                    sfi_IC2=numeric(rep),
                    sfi_IC3=numeric(rep))
TIME1 = data.frame(sparse=numeric(rep),sparse40=numeric(rep),fast=numeric(rep),fast40=numeric(rep),infomax=numeric(rep),infomax40=numeric(rep),EBM=numeric(rep),sfi=numeric(rep))

# rep 1-20
load("../../Data/highdim_sparse_1_20.RData")
for (i in 1:20) {
  xmat = xmat_list[[i]]
  smat = smat_list[[i]]
  mmat = mmat_list[[i]]
  
  # standardize xmat
  xmat_center = scale(xmat,center = T,scale = F)
  
  select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec",method = "C",
                                  lambda = sqrt(2)/2, eps = my_eps,maxit = my_maxit, 
                                  BIC_plot = F,nu_list = nu_list)
  my_nu = select_sparseICA$best_nu
  
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 1,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$sparse_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$sparse_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse_IC1[i] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse_IC2[i] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse_IC3[i] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse[i] = as.numeric(sparse.end.time - sparse.start.time)
  
  cat("#### Sparse ICA-1 finish! ####\n")
  
  sparse40.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 40,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse40.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse40_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$sparse40_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$sparse40_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse40_IC1[i] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse40_IC2[i] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse40_IC3[i] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse40[i] = as.numeric(sparse40.end.time - sparse40.start.time)
  
  cat("#### Sparse ICA-40 finish! ####\n")
  
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$infomax_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$infomax_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax_IC1[i] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax_IC2[i] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax_IC3[i] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax[i] = as.numeric(infomax.end.time - infomax.start.time)
  cat("#### Infomax ICA-1 finish! ####\n")
  
  infomax40.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax40.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax40_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$infomax40_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$infomax40_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax40_IC1[i] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax40_IC2[i] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax40_IC3[i] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax40[i] = as.numeric(infomax40.end.time - infomax40.start.time)
  cat("#### Infomax ICA-40 finish! ####\n")
  
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C")
  fast.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$fast_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$fast_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast_IC1[i] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast_IC2[i] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast_IC3[i] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast[i] = as.numeric(fast.end.time - fast.start.time)
  cat("#### Fast ICA-1 finish! ####\n")
  
  fast40.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C",restarts = 40)
  fast40.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast40_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$fast40_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$fast40_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast40_IC1[i] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast40_IC2[i] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast40_IC3[i] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast40[i] = as.numeric(fast40.end.time - fast40.start.time)
  cat("#### Fast ICA-40 finish! ####\n")
  
  ###############################################################################
  filename = paste0("../../Results/SICA_EBM/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(0,nrow = 59412,ncol = 117))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  
  S_cor1$EBM_IC1[i] = cor(my_S_EBM[,1],smat[,1])
  S_cor1$EBM_IC2[i] = cor(my_S_EBM[,2],smat[,2])
  S_cor1$EBM_IC3[i] = cor(my_S_EBM[,3],smat[,3])
  
  M_cor1$EBM_IC1[i] = cor(my_M_EBM[1,],mmat[1,])
  M_cor1$EBM_IC2[i] = cor(my_M_EBM[2,],mmat[2,])
  M_cor1$EBM_IC3[i] = cor(my_M_EBM[3,],mmat[3,])
  TIME1$EBM[i] = as.numeric(dat$tEnd)
  cat("#### SICA-EBM finish! ####\n")
  ###############################################################################
  filename = paste0("../../Results/sparsefast/estS_",i,".mat")
  dat = readMat(filename)
  my_S_sfi = matchICA(t(dat$myS),smat)
  my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_cor1$sfi_IC1[i] = cor(my_S_sfi[,1],smat[,1])
  S_cor1$sfi_IC2[i] = cor(my_S_sfi[,2],smat[,2])
  S_cor1$sfi_IC3[i] = cor(my_S_sfi[,3],smat[,3])
  
  M_cor1$sfi_IC1[i] = cor(my_M_sfi[1,],mmat[1,])
  M_cor1$sfi_IC2[i] = cor(my_M_sfi[2,],mmat[2,])
  M_cor1$sfi_IC3[i] = cor(my_M_sfi[3,],mmat[3,])
  TIME1$sfi[i] = as.numeric(dat$tEnd)
  cat("#### Sparse Fast ICA finish! ####\n")
  
  cat("############################ ",i," th finish! ##############################\n")
}


# rep 21-40
load("../../Data/highdim_sparse_21_40.RData")
gap = 20
for (i in 1:20) {
  xmat = xmat_list[[i]]
  smat = smat_list[[i]]
  mmat = mmat_list[[i]]
  
  # standardize xmat
  xmat_center = scale(xmat,center = T,scale = F)
  
  select_sparseICA= BIC_sparseICA(xData = xmat, n.comp = 3,whiten = "eigenvec",method = "C",
                                  lambda = sqrt(2)/2, eps = my_eps,maxit = my_maxit, 
                                  BIC_plot = F,nu_list = nu_list)
  my_nu = select_sparseICA$best_nu
  
  sparse.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 1,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$sparse_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$sparse_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse_IC1[i] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse_IC2[i] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse_IC3[i] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse[i] = as.numeric(sparse.end.time - sparse.start.time)
  
  cat("#### Sparse ICA-1 finish! ####\n")
  
  sparse40.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 40,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse40.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse40_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$sparse40_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$sparse40_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse40_IC1[i] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse40_IC2[i] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse40_IC3[i] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse40[i] = as.numeric(sparse40.end.time - sparse40.start.time)
  
  cat("#### Sparse ICA-40 finish! ####\n")
  
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$infomax_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$infomax_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax_IC1[i] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax_IC2[i] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax_IC3[i] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax[i] = as.numeric(infomax.end.time - infomax.start.time)
  cat("#### Infomax ICA-1 finish! ####\n")
  
  infomax40.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax40.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax40_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$infomax40_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$infomax40_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax40_IC1[i] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax40_IC2[i] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax40_IC3[i] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax40[i] = as.numeric(infomax40.end.time - infomax40.start.time)
  cat("#### Infomax ICA-40 finish! ####\n")
  
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C")
  fast.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$fast_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$fast_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast_IC1[i] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast_IC2[i] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast_IC3[i] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast[i] = as.numeric(fast.end.time - fast.start.time)
  cat("#### Fast ICA-1 finish! ####\n")
  
  fast40.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C",restarts = 40)
  fast40.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast40_IC1[i] = cor(S_match[,1],smat[,1])
  S_cor1$fast40_IC2[i] = cor(S_match[,2],smat[,2])
  S_cor1$fast40_IC3[i] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast40_IC1[i] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast40_IC2[i] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast40_IC3[i] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast40[i] = as.numeric(fast40.end.time - fast40.start.time)
  cat("#### Fast ICA-40 finish! ####\n")
  
  ###############################################################################
  filename = paste0("../../Results/SICA_EBM/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(0,nrow = 59412,ncol = 117))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  
  S_cor1$EBM_IC1[i] = cor(my_S_EBM[,1],smat[,1])
  S_cor1$EBM_IC2[i] = cor(my_S_EBM[,2],smat[,2])
  S_cor1$EBM_IC3[i] = cor(my_S_EBM[,3],smat[,3])
  
  M_cor1$EBM_IC1[i] = cor(my_M_EBM[1,],mmat[1,])
  M_cor1$EBM_IC2[i] = cor(my_M_EBM[2,],mmat[2,])
  M_cor1$EBM_IC3[i] = cor(my_M_EBM[3,],mmat[3,])
  TIME1$EBM[i] = as.numeric(dat$tEnd)
  cat("#### SICA-EBM finish! ####\n")
  ###############################################################################
  filename = paste0("../../Results/sparsefast/estS_",i,".mat")
  dat = readMat(filename)
  my_S_sfi = matchICA(t(dat$myS),smat)
  my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_cor1$sfi_IC1[i] = cor(my_S_sfi[,1],smat[,1])
  S_cor1$sfi_IC2[i] = cor(my_S_sfi[,2],smat[,2])
  S_cor1$sfi_IC3[i] = cor(my_S_sfi[,3],smat[,3])
  
  M_cor1$sfi_IC1[i] = cor(my_M_sfi[1,],mmat[1,])
  M_cor1$sfi_IC2[i] = cor(my_M_sfi[2,],mmat[2,])
  M_cor1$sfi_IC3[i] = cor(my_M_sfi[3,],mmat[3,])
  TIME1$sfi[i] = as.numeric(dat$tEnd)
  cat("#### Sparse Fast ICA finish! ####\n")
  
  cat("############################ ",i," th finish! ##############################\n")
}










save(S_cor1,M_cor1,TIME1,file = "../Results/13_simreal_all_BIC.RData")

#################################################################################################
load("../Results/13_simreal_all_BIC.RData")

library(ggplot2)
library(ggpubr)

#################################################################################
# plot single version first

# For S estimates
S_cor = data.frame(cor=c(S_cor1$sparse_IC1,S_cor1$fast_IC1,S_cor1$infomax_IC1,S_cor1$EBM_IC1,S_cor1$sfi_IC1,
                         S_cor1$sparse_IC2,S_cor1$fast_IC2,S_cor1$infomax_IC2,S_cor1$EBM_IC2,S_cor1$sfi_IC1,
                         S_cor1$sparse_IC3,S_cor1$fast_IC3,S_cor1$infomax_IC3,S_cor1$EBM_IC3,S_cor1$sfi_IC3),
                    Method=c(rep("Sparse ICA",length(S_cor1$sparse_IC1)),rep("Fast ICA",length(S_cor1$fast_IC1)),rep("Infomax ICA",length(S_cor1$infomax_IC1)),rep("SICA-EBM",length(S_cor1$EBM_IC1)),rep("Sparse Fast ICA",length(S_cor1$sfi_IC1)),
                             rep("Sparse ICA",length(S_cor1$sparse_IC2)),rep("Fast ICA",length(S_cor1$fast_IC2)),rep("Infomax ICA",length(S_cor1$infomax_IC2)),rep("SICA-EBM",length(S_cor1$EBM_IC2)),rep("Sparse Fast ICA",length(S_cor1$sfi_IC2)),
                             rep("Sparse ICA",length(S_cor1$sparse_IC3)),rep("Fast ICA",length(S_cor1$fast_IC3)),rep("Infomax ICA",length(S_cor1$infomax_IC3)),rep("SICA-EBM",length(S_cor1$EBM_IC3)),rep("Sparse Fast ICA",length(S_cor1$sfi_IC3))
                    ),
                    IC=c(rep("IC1",5*length(S_cor1$sparse_IC1)),
                          rep("IC2",5*length(S_cor1$sparse_IC2)),
                          rep("IC3",5*length(S_cor1$sparse_IC3)))
)

p1 = ggplot(S_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method)+ 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Source Signal")#+ylim(c(0.9,1))


# For M estimates
M_cor = data.frame(cor=c(M_cor1$sparse_IC1,M_cor1$fast_IC1,M_cor1$infomax_IC1,M_cor1$EBM_IC1,M_cor1$sfi_IC1,
                         M_cor1$sparse_IC2,M_cor1$fast_IC2,M_cor1$infomax_IC2,M_cor1$sfi_IC2,M_cor1$sfi_IC2,
                         M_cor1$sparse_IC3,M_cor1$fast_IC3,M_cor1$infomax_IC3,M_cor1$EBM_IC3,M_cor1$sfi_IC3),
                   Method=c(rep("Sparse ICA",length(M_cor1$sparse_IC1)),rep("Fast ICA",length(M_cor1$fast_IC1)),rep("Infomax ICA",length(M_cor1$infomax_IC1)),rep("SICA-EBM",length(M_cor1$EBM_IC1)),rep("Sparse Fast ICA",length(M_cor1$sfi_IC1)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC2)),rep("Fast ICA",length(M_cor1$fast_IC2)),rep("Infomax ICA",length(M_cor1$infomax_IC2)),rep("SICA-EBM",length(M_cor1$EBM_IC2)),rep("Sparse Fast ICA",length(M_cor1$sfi_IC2)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC3)),rep("Fast ICA",length(M_cor1$fast_IC3)),rep("Infomax ICA",length(M_cor1$infomax_IC3)),rep("SICA-EBM",length(M_cor1$EBM_IC3)),rep("Sparse Fast ICA",length(M_cor1$sfi_IC3))
                   ),
                   IC=c(rep("IC1",5*length(M_cor1$sparse_IC1)),
                        rep("IC2",5*length(M_cor1$sparse_IC2)),
                        rep("IC3",5*length(M_cor1$sparse_IC3)))
)


p2 = ggplot(M_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method) + 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Mixing Matrix")#+ylim(c(0.985,1))


ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../Figures/13_single_BIC.pdf",width = 14, height = 6)

#################################################################################
# plot no SICA-EBM and SFI version

# For S estimates
S_cor = data.frame(cor=c(S_cor1$sparse_IC1,S_cor1$fast_IC1,S_cor1$infomax_IC1,
                         S_cor1$sparse_IC2,S_cor1$fast_IC2,S_cor1$infomax_IC2,
                         S_cor1$sparse_IC3,S_cor1$fast_IC3,S_cor1$infomax_IC3),
                   Method=c(rep("Sparse ICA",length(S_cor1$sparse_IC1)),rep("Fast ICA",length(S_cor1$fast_IC1)),rep("Infomax ICA",length(S_cor1$infomax_IC1)),
                            rep("Sparse ICA",length(S_cor1$sparse_IC2)),rep("Fast ICA",length(S_cor1$fast_IC2)),rep("Infomax ICA",length(S_cor1$infomax_IC2)),
                            rep("Sparse ICA",length(S_cor1$sparse_IC3)),rep("Fast ICA",length(S_cor1$fast_IC3)),rep("Infomax ICA",length(S_cor1$infomax_IC3))
                   ),
                   IC=c(rep("IC1",3*length(S_cor1$sparse_IC1)),
                        rep("IC2",3*length(S_cor1$sparse_IC2)),
                        rep("IC3",3*length(S_cor1$sparse_IC3)))
)

p1 = ggplot(S_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method)+ 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Source Signal")#+ylim(c(0.9,1))


# For M estimates
M_cor = data.frame(cor=c(M_cor1$sparse_IC1,M_cor1$fast_IC1,M_cor1$infomax_IC1,
                         M_cor1$sparse_IC2,M_cor1$fast_IC2,M_cor1$infomax_IC2,
                         M_cor1$sparse_IC3,M_cor1$fast_IC3,M_cor1$infomax_IC3),
                   Method=c(rep("Sparse ICA",length(M_cor1$sparse_IC1)),rep("Fast ICA",length(M_cor1$fast_IC1)),rep("Infomax ICA",length(M_cor1$infomax_IC1)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC2)),rep("Fast ICA",length(M_cor1$fast_IC2)),rep("Infomax ICA",length(M_cor1$infomax_IC2)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC3)),rep("Fast ICA",length(M_cor1$fast_IC3)),rep("Infomax ICA",length(M_cor1$infomax_IC3))
                   ),
                   IC=c(rep("IC1",3*length(M_cor1$sparse_IC1)),
                        rep("IC2",3*length(M_cor1$sparse_IC2)),
                        rep("IC3",3*length(M_cor1$sparse_IC3)))
)


p2 = ggplot(M_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method) + 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Mixing Matrix")#+ylim(c(0.985,1))


ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../Figures/13_single_nomatlab_BIC.pdf",width = 14, height = 6)


#################################################################################
# plot single and 40 restarts

# For S estimates
S_cor = data.frame(cor=c(S_cor1$sparse_IC1,S_cor1$fast_IC1,S_cor1$infomax_IC1,S_cor1$sparse40_IC1,S_cor1$fast40_IC1,S_cor1$infomax40_IC1,
                         S_cor1$sparse_IC2,S_cor1$fast_IC2,S_cor1$infomax_IC2,S_cor1$sparse40_IC2,S_cor1$fast40_IC2,S_cor1$infomax40_IC2,
                         S_cor1$sparse_IC3,S_cor1$fast_IC3,S_cor1$infomax_IC3,S_cor1$sparse40_IC3,S_cor1$fast40_IC3,S_cor1$infomax40_IC3),
                   Method=c(rep("Sparse ICA-1",length(S_cor1$sparse_IC1)),rep("Fast ICA-1",length(S_cor1$fast_IC1)),rep("Infomax ICA-1",length(S_cor1$infomax_IC1)),rep("Sparse ICA-40",length(S_cor1$sparse40_IC1)),rep("Fast ICA-40",length(S_cor1$fast40_IC1)),rep("Infomax ICA-40",length(S_cor1$infomax40_IC1)),
                            rep("Sparse ICA-1",length(S_cor1$sparse_IC2)),rep("Fast ICA-1",length(S_cor1$fast_IC2)),rep("Infomax ICA-1",length(S_cor1$infomax_IC2)),rep("Sparse ICA-40",length(S_cor1$sparse40_IC2)),rep("Fast ICA-40",length(S_cor1$fast40_IC2)),rep("Infomax ICA-40",length(S_cor1$infomax40_IC2)),
                            rep("Sparse ICA-1",length(S_cor1$sparse_IC3)),rep("Fast ICA-1",length(S_cor1$fast_IC3)),rep("Infomax ICA-1",length(S_cor1$infomax_IC3)),rep("Sparse ICA-40",length(S_cor1$sparse40_IC3)),rep("Fast ICA-40",length(S_cor1$fast40_IC3)),rep("Infomax ICA-40",length(S_cor1$infomax40_IC3))
                   ),
                   IC=c(rep("IC1",6*length(S_cor1$sparse_IC1)),
                        rep("IC2",6*length(S_cor1$sparse_IC2)),
                        rep("IC3",6*length(S_cor1$sparse_IC3)))
)

p1 = ggplot(S_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method)+ 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Source Signal")#+ylim(c(0.9,1))


# For M estimates
M_cor = data.frame(cor=c(M_cor1$sparse_IC1,M_cor1$fast_IC1,M_cor1$infomax_IC1,M_cor1$sparse40_IC1,M_cor1$fast40_IC1,M_cor1$infomax40_IC1,
                         M_cor1$sparse_IC2,M_cor1$fast_IC2,M_cor1$infomax_IC2,M_cor1$sparse40_IC2,M_cor1$fast40_IC2,M_cor1$infomax40_IC2,
                         M_cor1$sparse_IC3,M_cor1$fast_IC3,M_cor1$infomax_IC3,M_cor1$sparse40_IC3,M_cor1$fast40_IC3,M_cor1$infomax40_IC3),
                   Method=c(rep("Sparse ICA-1",length(M_cor1$sparse_IC1)),rep("Fast ICA-1",length(M_cor1$fast_IC1)),rep("Infomax ICA-1",length(M_cor1$infomax_IC1)),rep("Sparse ICA-40",length(M_cor1$sparse40_IC1)),rep("Fast ICA-40",length(M_cor1$fast40_IC1)),rep("Infomax ICA-40",length(M_cor1$infomax40_IC1)),
                            rep("Sparse ICA-1",length(M_cor1$sparse_IC2)),rep("Fast ICA-1",length(M_cor1$fast_IC2)),rep("Infomax ICA-1",length(M_cor1$infomax_IC2)),rep("Sparse ICA-40",length(M_cor1$sparse40_IC2)),rep("Fast ICA-40",length(M_cor1$fast40_IC2)),rep("Infomax ICA-40",length(M_cor1$infomax40_IC2)),
                            rep("Sparse ICA-1",length(M_cor1$sparse_IC3)),rep("Fast ICA-1",length(M_cor1$fast_IC3)),rep("Infomax ICA-1",length(M_cor1$infomax_IC3)),rep("Sparse ICA-40",length(M_cor1$sparse40_IC3)),rep("Fast ICA-40",length(M_cor1$fast40_IC3)),rep("Infomax ICA-40",length(M_cor1$infomax40_IC3))
                   ),
                   IC=c(rep("IC1",6*length(M_cor1$sparse_IC1)),
                        rep("IC2",6*length(M_cor1$sparse_IC2)),
                        rep("IC3",6*length(M_cor1$sparse_IC3)))
)


p2 = ggplot(M_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method) + 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Mixing Matrix")#+ylim(c(0.985,1))


ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../Figures/13_single_40restarts_nomatlab_BIC.pdf",width = 14, height = 6)


apply(TIME1,2,mean)

mean(S_cor1$sparse_IC1)
sd(S_cor1$sparse_IC1)
mean(S_cor1$fast_IC1)
sd(S_cor1$fast_IC1)
mean(S_cor1$infomax_IC1)
sd(S_cor1$infomax_IC1)
mean(S_cor1$EBM_IC1)
sd(S_cor1$EBM_IC1)
mean(S_cor1$sfi_IC1)
sd(S_cor1$sfi_IC1)

mean(S_cor1$sparse_IC2)
sd(S_cor1$sparse_IC2)
mean(S_cor1$fast_IC2)
sd(S_cor1$fast_IC2)
mean(S_cor1$infomax_IC2)
sd(S_cor1$infomax_IC2)
mean(S_cor1$EBM_IC2)
sd(S_cor1$EBM_IC2)
mean(S_cor1$sfi_IC2)
sd(S_cor1$sfi_IC2)

mean(S_cor1$sparse_IC3)
sd(S_cor1$sparse_IC3)
mean(S_cor1$fast_IC3)
sd(S_cor1$fast_IC3)
mean(S_cor1$infomax_IC3)
sd(S_cor1$infomax_IC3)
mean(S_cor1$EBM_IC3)
sd(S_cor1$EBM_IC3)
mean(S_cor1$sfi_IC3)
sd(S_cor1$sfi_IC3)


mean(M_cor1$sparse_IC1)
sd(M_cor1$sparse_IC1)
mean(M_cor1$fast_IC1)
sd(M_cor1$fast_IC1)
mean(M_cor1$infomax_IC1)
sd(M_cor1$infomax_IC1)
mean(M_cor1$EBM_IC1)
sd(M_cor1$EBM_IC1)
mean(M_cor1$sfi_IC1)
sd(M_cor1$sfi_IC1)

mean(M_cor1$sparse_IC2)
sd(M_cor1$sparse_IC2)
mean(M_cor1$fast_IC2)
sd(M_cor1$fast_IC2)
mean(M_cor1$infomax_IC2)
sd(M_cor1$infomax_IC2)
mean(M_cor1$EBM_IC2)
sd(M_cor1$EBM_IC2)
mean(M_cor1$sfi_IC2)
sd(M_cor1$sfi_IC2)

mean(M_cor1$sparse_IC3)
sd(M_cor1$sparse_IC3)
mean(M_cor1$fast_IC3)
sd(M_cor1$fast_IC3)
mean(M_cor1$infomax_IC3)
sd(M_cor1$infomax_IC3)
mean(M_cor1$EBM_IC3)
sd(M_cor1$EBM_IC3)
mean(M_cor1$sfi_IC3)
sd(M_cor1$sfi_IC3)

boxplot(TIME1$sparse,TIME1$fast,TIME1$infomax,TIME1$sparse40,TIME1$fast40,TIME1$infomax40,
        names = c("Sparse ICA","Fast ICA","Infomax ICA","Sparse ICA-40","Fast ICA-40","Infomax ICA-40"),ylab="TIME")
summary(TIME1$sparse)
summary(TIME1$fast)
summary(TIME1$infomax)

hist(TIME1$sparse)
hist(TIME1$fast)
hist(TIME1$infomax)

sd(TIME1$sparse)
sd(TIME1$fast)
sd(TIME1$infomax)

boxplot(TIME1$sparse,TIME1$fast,TIME1$infomax)

apply(TIME1, 2, mean)
apply(TIME1, 2, median)

mean(TIME1$EBM)


###########################################################################################
# make BIC plot for realistic (high-dimensional) simulations
set.seed(2123)
simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_iter5,template_xifti = template_image,template_time = sample_time)
xmat = simData$X
smat = simData$S
mmat = simData$Ms

nu_list = seq(0.1,4,0.1)
pdf("../Figures/realistic_BIC_plot.pdf",width = 8,height = 4)
select_sparseICA= BIC_sparseICA_Rcpp(xData = xmat, n.comp = 3,U.list=NULL,whiten = "eigenvec", 
                                     restarts.pbyd = 30, lambda = sqrt(2)/2, eps = my_eps,maxit.laplace = my_maxit, 
                                     BIC_plot = T,nu_list = nu_list)
my_nu = select_sparseICA$best_nu
dev.off()


