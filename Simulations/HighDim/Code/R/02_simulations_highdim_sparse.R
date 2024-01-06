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
  
  M_cor1$sparse_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse[i+gap] = as.numeric(sparse.end.time - sparse.start.time)
  
  cat("#### Sparse ICA-1 finish! ####\n")
  
  sparse40.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 40,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse40.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$sparse40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$sparse40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse40_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse40_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse40_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse40[i+gap] = as.numeric(sparse40.end.time - sparse40.start.time)
  
  cat("#### Sparse ICA-40 finish! ####\n")
  
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax[i+gap] = as.numeric(infomax.end.time - infomax.start.time)
  cat("#### Infomax ICA-1 finish! ####\n")
  
  infomax40.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax40.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax40_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax40_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax40_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax40[i+gap] = as.numeric(infomax40.end.time - infomax40.start.time)
  cat("#### Infomax ICA-40 finish! ####\n")
  
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C")
  fast.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast[i+gap] = as.numeric(fast.end.time - fast.start.time)
  cat("#### Fast ICA-1 finish! ####\n")
  
  fast40.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C",restarts = 40)
  fast40.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast40_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast40_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast40_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast40[i+gap] = as.numeric(fast40.end.time - fast40.start.time)
  cat("#### Fast ICA-40 finish! ####\n")
  
  ###############################################################################
  filename = paste0("../../Results/SICA_EBM/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(0,nrow = 59412,ncol = 117))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  
  S_cor1$EBM_IC1[i+gap] = cor(my_S_EBM[,1],smat[,1])
  S_cor1$EBM_IC2[i+gap] = cor(my_S_EBM[,2],smat[,2])
  S_cor1$EBM_IC3[i+gap] = cor(my_S_EBM[,3],smat[,3])
  
  M_cor1$EBM_IC1[i+gap] = cor(my_M_EBM[1,],mmat[1,])
  M_cor1$EBM_IC2[i+gap] = cor(my_M_EBM[2,],mmat[2,])
  M_cor1$EBM_IC3[i+gap] = cor(my_M_EBM[3,],mmat[3,])
  TIME1$EBM[i+gap] = as.numeric(dat$tEnd)
  cat("#### SICA-EBM finish! ####\n")
  ###############################################################################
  filename = paste0("../../Results/sparsefast/estS_",i,".mat")
  dat = readMat(filename)
  my_S_sfi = matchICA(t(dat$myS),smat)
  my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_cor1$sfi_IC1[i+gap] = cor(my_S_sfi[,1],smat[,1])
  S_cor1$sfi_IC2[i+gap] = cor(my_S_sfi[,2],smat[,2])
  S_cor1$sfi_IC3[i+gap] = cor(my_S_sfi[,3],smat[,3])
  
  M_cor1$sfi_IC1[i+gap] = cor(my_M_sfi[1,],mmat[1,])
  M_cor1$sfi_IC2[i+gap] = cor(my_M_sfi[2,],mmat[2,])
  M_cor1$sfi_IC3[i+gap] = cor(my_M_sfi[3,],mmat[3,])
  TIME1$sfi[i+gap] = as.numeric(dat$tEnd)
  cat("#### Sparse Fast ICA finish! ####\n")
  
  cat("############################ ",i," th finish! ##############################\n")
}


# rep 41 - 60
load("../../Data/highdim_sparse_41_60.RData")
gap = 40
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
  
  M_cor1$sparse_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse[i+gap] = as.numeric(sparse.end.time - sparse.start.time)
  
  cat("#### Sparse ICA-1 finish! ####\n")
  
  sparse40.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 40,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse40.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$sparse40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$sparse40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse40_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse40_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse40_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse40[i+gap] = as.numeric(sparse40.end.time - sparse40.start.time)
  
  cat("#### Sparse ICA-40 finish! ####\n")
  
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax[i+gap] = as.numeric(infomax.end.time - infomax.start.time)
  cat("#### Infomax ICA-1 finish! ####\n")
  
  infomax40.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax40.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax40_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax40_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax40_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax40[i+gap] = as.numeric(infomax40.end.time - infomax40.start.time)
  cat("#### Infomax ICA-40 finish! ####\n")
  
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C")
  fast.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast[i+gap] = as.numeric(fast.end.time - fast.start.time)
  cat("#### Fast ICA-1 finish! ####\n")
  
  fast40.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C",restarts = 40)
  fast40.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast40_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast40_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast40_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast40[i+gap] = as.numeric(fast40.end.time - fast40.start.time)
  cat("#### Fast ICA-40 finish! ####\n")
  
  ###############################################################################
  filename = paste0("../../Results/SICA_EBM/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(0,nrow = 59412,ncol = 117))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  
  S_cor1$EBM_IC1[i+gap] = cor(my_S_EBM[,1],smat[,1])
  S_cor1$EBM_IC2[i+gap] = cor(my_S_EBM[,2],smat[,2])
  S_cor1$EBM_IC3[i+gap] = cor(my_S_EBM[,3],smat[,3])
  
  M_cor1$EBM_IC1[i+gap] = cor(my_M_EBM[1,],mmat[1,])
  M_cor1$EBM_IC2[i+gap] = cor(my_M_EBM[2,],mmat[2,])
  M_cor1$EBM_IC3[i+gap] = cor(my_M_EBM[3,],mmat[3,])
  TIME1$EBM[i+gap] = as.numeric(dat$tEnd)
  cat("#### SICA-EBM finish! ####\n")
  ###############################################################################
  filename = paste0("../../Results/sparsefast/estS_",i,".mat")
  dat = readMat(filename)
  my_S_sfi = matchICA(t(dat$myS),smat)
  my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_cor1$sfi_IC1[i+gap] = cor(my_S_sfi[,1],smat[,1])
  S_cor1$sfi_IC2[i+gap] = cor(my_S_sfi[,2],smat[,2])
  S_cor1$sfi_IC3[i+gap] = cor(my_S_sfi[,3],smat[,3])
  
  M_cor1$sfi_IC1[i+gap] = cor(my_M_sfi[1,],mmat[1,])
  M_cor1$sfi_IC2[i+gap] = cor(my_M_sfi[2,],mmat[2,])
  M_cor1$sfi_IC3[i+gap] = cor(my_M_sfi[3,],mmat[3,])
  TIME1$sfi[i+gap] = as.numeric(dat$tEnd)
  cat("#### Sparse Fast ICA finish! ####\n")
  
  cat("############################ ",i," th finish! ##############################\n")
}


# rep 61 - 80
load("../../Data/highdim_sparse_61_80.RData")
gap = 60
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
  
  M_cor1$sparse_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse[i+gap] = as.numeric(sparse.end.time - sparse.start.time)
  
  cat("#### Sparse ICA-1 finish! ####\n")
  
  sparse40.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 40,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse40.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$sparse40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$sparse40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse40_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse40_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse40_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse40[i+gap] = as.numeric(sparse40.end.time - sparse40.start.time)
  
  cat("#### Sparse ICA-40 finish! ####\n")
  
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax[i+gap] = as.numeric(infomax.end.time - infomax.start.time)
  cat("#### Infomax ICA-1 finish! ####\n")
  
  infomax40.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax40.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax40_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax40_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax40_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax40[i+gap] = as.numeric(infomax40.end.time - infomax40.start.time)
  cat("#### Infomax ICA-40 finish! ####\n")
  
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C")
  fast.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast[i+gap] = as.numeric(fast.end.time - fast.start.time)
  cat("#### Fast ICA-1 finish! ####\n")
  
  fast40.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C",restarts = 40)
  fast40.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast40_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast40_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast40_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast40[i+gap] = as.numeric(fast40.end.time - fast40.start.time)
  cat("#### Fast ICA-40 finish! ####\n")
  
  ###############################################################################
  filename = paste0("../../Results/SICA_EBM/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(0,nrow = 59412,ncol = 117))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  
  S_cor1$EBM_IC1[i+gap] = cor(my_S_EBM[,1],smat[,1])
  S_cor1$EBM_IC2[i+gap] = cor(my_S_EBM[,2],smat[,2])
  S_cor1$EBM_IC3[i+gap] = cor(my_S_EBM[,3],smat[,3])
  
  M_cor1$EBM_IC1[i+gap] = cor(my_M_EBM[1,],mmat[1,])
  M_cor1$EBM_IC2[i+gap] = cor(my_M_EBM[2,],mmat[2,])
  M_cor1$EBM_IC3[i+gap] = cor(my_M_EBM[3,],mmat[3,])
  TIME1$EBM[i+gap] = as.numeric(dat$tEnd)
  cat("#### SICA-EBM finish! ####\n")
  ###############################################################################
  filename = paste0("../../Results/sparsefast/estS_",i,".mat")
  dat = readMat(filename)
  my_S_sfi = matchICA(t(dat$myS),smat)
  my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_cor1$sfi_IC1[i+gap] = cor(my_S_sfi[,1],smat[,1])
  S_cor1$sfi_IC2[i+gap] = cor(my_S_sfi[,2],smat[,2])
  S_cor1$sfi_IC3[i+gap] = cor(my_S_sfi[,3],smat[,3])
  
  M_cor1$sfi_IC1[i+gap] = cor(my_M_sfi[1,],mmat[1,])
  M_cor1$sfi_IC2[i+gap] = cor(my_M_sfi[2,],mmat[2,])
  M_cor1$sfi_IC3[i+gap] = cor(my_M_sfi[3,],mmat[3,])
  TIME1$sfi[i+gap] = as.numeric(dat$tEnd)
  cat("#### Sparse Fast ICA finish! ####\n")
  
  cat("############################ ",i," th finish! ##############################\n")
}


# rep 81 - 100
load("../../Data/highdim_sparse_81_100.RData")
gap = 80
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
  
  M_cor1$sparse_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse[i+gap] = as.numeric(sparse.end.time - sparse.start.time)
  
  cat("#### Sparse ICA-1 finish! ####\n")
  
  sparse40.start.time = Sys.time()
  my_sparseICA = sparseICA(xData = xmat, n.comp = 3,method = "C",
                           whiten = "eigenvec", restarts = 40,nu = my_nu,
                           eps = my_eps,maxit = my_maxit, converge_plot = F,verbose = F)
  sparse40.end.time = Sys.time()
  S_match=matchICA(my_sparseICA$estS,smat)
  my_M_sparse = est.M.ols(S_match,xmat_center)
  
  S_cor1$sparse40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$sparse40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$sparse40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$sparse40_IC1[i+gap] = cor(my_M_sparse[1,],mmat[1,])
  M_cor1$sparse40_IC2[i+gap] = cor(my_M_sparse[2,],mmat[2,])
  M_cor1$sparse40_IC3[i+gap] = cor(my_M_sparse[3,],mmat[3,])
  TIME1$sparse40[i+gap] = as.numeric(sparse40.end.time - sparse40.start.time)
  
  cat("#### Sparse ICA-40 finish! ####\n")
  
  infomax.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax[i+gap] = as.numeric(infomax.end.time - infomax.start.time)
  cat("#### Infomax ICA-1 finish! ####\n")
  
  infomax40.start.time = Sys.time()
  my_infomaxICA = infomaxICA(X=xmat,n.comp = 3,whiten = TRUE,eps = my_eps,maxit = my_maxit,restarts = 1)
  infomax40.end.time = Sys.time()
  S_match=matchICA(my_infomaxICA$S,smat)
  my_M_infomax = est.M.ols(S_match,xmat_center)
  
  S_cor1$infomax40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$infomax40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$infomax40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$infomax40_IC1[i+gap] = cor(my_M_infomax[1,],mmat[1,])
  M_cor1$infomax40_IC2[i+gap] = cor(my_M_infomax[2,],mmat[2,])
  M_cor1$infomax40_IC3[i+gap] = cor(my_M_infomax[3,],mmat[3,])
  TIME1$infomax40[i+gap] = as.numeric(infomax40.end.time - infomax40.start.time)
  cat("#### Infomax ICA-40 finish! ####\n")
  
  fast.start.time = Sys.time()
  my_fastICA = fastICA(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C")
  fast.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast[i+gap] = as.numeric(fast.end.time - fast.start.time)
  cat("#### Fast ICA-1 finish! ####\n")
  
  fast40.start.time = Sys.time()
  my_fastICA = fastICA_restarts(xmat,n.comp = 3,maxit = my_maxit,tol=my_eps,method="C",restarts = 40)
  fast40.end.time = Sys.time()
  S_match=matchICA(my_fastICA$S,smat)
  my_M_fast = est.M.ols(S_match,xmat_center)
  
  S_cor1$fast40_IC1[i+gap] = cor(S_match[,1],smat[,1])
  S_cor1$fast40_IC2[i+gap] = cor(S_match[,2],smat[,2])
  S_cor1$fast40_IC3[i+gap] = cor(S_match[,3],smat[,3])
  
  M_cor1$fast40_IC1[i+gap] = cor(my_M_fast[1,],mmat[1,])
  M_cor1$fast40_IC2[i+gap] = cor(my_M_fast[2,],mmat[2,])
  M_cor1$fast40_IC3[i+gap] = cor(my_M_fast[3,],mmat[3,])
  TIME1$fast40[i+gap] = as.numeric(fast40.end.time - fast40.start.time)
  cat("#### Fast ICA-40 finish! ####\n")
  
  ###############################################################################
  filename = paste0("../../Results/SICA_EBM/estS_",i,".mat")
  dat = readMat(filename)
  temp = cbind(smat,matrix(0,nrow = 59412,ncol = 117))
  my_S_EBM = matchICA(t(dat$myS),temp)
  my_S_EBM = my_S_EBM[,1:3]
  my_M_EBM = est.M.ols(my_S_EBM,xmat_center)
  
  S_cor1$EBM_IC1[i+gap] = cor(my_S_EBM[,1],smat[,1])
  S_cor1$EBM_IC2[i+gap] = cor(my_S_EBM[,2],smat[,2])
  S_cor1$EBM_IC3[i+gap] = cor(my_S_EBM[,3],smat[,3])
  
  M_cor1$EBM_IC1[i+gap] = cor(my_M_EBM[1,],mmat[1,])
  M_cor1$EBM_IC2[i+gap] = cor(my_M_EBM[2,],mmat[2,])
  M_cor1$EBM_IC3[i+gap] = cor(my_M_EBM[3,],mmat[3,])
  TIME1$EBM[i+gap] = as.numeric(dat$tEnd)
  cat("#### SICA-EBM finish! ####\n")
  ###############################################################################
  filename = paste0("../../Results/sparsefast/estS_",i,".mat")
  dat = readMat(filename)
  my_S_sfi = matchICA(t(dat$myS),smat)
  my_M_sfi = est.M.ols(my_S_sfi,xmat_center)
  S_cor1$sfi_IC1[i+gap] = cor(my_S_sfi[,1],smat[,1])
  S_cor1$sfi_IC2[i+gap] = cor(my_S_sfi[,2],smat[,2])
  S_cor1$sfi_IC3[i+gap] = cor(my_S_sfi[,3],smat[,3])
  
  M_cor1$sfi_IC1[i+gap] = cor(my_M_sfi[1,],mmat[1,])
  M_cor1$sfi_IC2[i+gap] = cor(my_M_sfi[2,],mmat[2,])
  M_cor1$sfi_IC3[i+gap] = cor(my_M_sfi[3,],mmat[3,])
  TIME1$sfi[i+gap] = as.numeric(dat$tEnd)
  cat("#### Sparse Fast ICA finish! ####\n")
  
  cat("############################ ",i," th finish! ##############################\n")
}

save(S_cor1,M_cor1,TIME1,file = "../../Results/02_highdim_sparse.RData")

