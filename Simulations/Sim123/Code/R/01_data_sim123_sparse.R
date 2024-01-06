###################################################################
# R codes for simulation setting 1 in single-subject simulations
# 01: generate simulation data
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

# library(SparseICA)
# library(steadyICA)
# library(fastICA)
# library(irlba)
library(R.matlab)

source("00_utils.R")

###################################################################
# General simulations settings
###################################################################

# The three given SNR levels
my_SNR = c(0.4,1.5,3)
# Number of replications
rep = 100
# Variance level of inactive elements in components
var_level = 0
# Random seeds
seed1 = 1470
seed2 = 2580
seed3 = 3690
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
# Generate simulation data for setting 1, 
# saved in RData and Matlab data format 
# can be used for R functions (fastICA, InfomaxICA, SparseICA) 
# and Matlab functions (SICE-EBM, SparseFastICA)
###################################################################

# low SNR
xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed1)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = var_level,noisyICA = my_noisyICA, snr=my_SNR[1],nTR = my_nTR,phi = my_phi)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  xmat_list[[i]] = xmat
  smat_list[[i]] = smat
  mmat_list[[i]] = mmat
  xmat_matlab_lowSNR = rbind(xmat_matlab_lowSNR,simData$X)
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/sim123_sparse_lowSNR.RData")
writeMat("../../Data/sim123_sparse_lowSNR.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

# medium SNR
xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_mediumSNR=c()
set.seed(seed2)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = var_level,noisyICA = my_noisyICA, snr=my_SNR[2],nTR = my_nTR,phi = my_phi)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  xmat_list[[i]] = xmat
  smat_list[[i]] = smat
  mmat_list[[i]] = mmat
  xmat_matlab_mediumSNR = rbind(xmat_matlab_mediumSNR,simData$X)
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/sim123_sparse_mediumSNR.RData")
writeMat("../../Data/sim123_sparse_mediumSNR.mat",xmat_matlab_mediumSNR=xmat_matlab_mediumSNR)

#high SNR
xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_highSNR=c()
set.seed(seed3)
for (i in 1:rep) {
  simData = SimFMRI123(var.inactive = var_level,noisyICA = my_noisyICA, snr=my_SNR[3],nTR = my_nTR,phi = my_phi)
  xmat = simData$X
  smat = simData$S
  mmat = simData$Ms
  xmat_list[[i]] = xmat
  smat_list[[i]] = smat
  mmat_list[[i]] = mmat
  xmat_matlab_highSNR = rbind(xmat_matlab_highSNR,simData$X)
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/sim123_sparse_highSNR.RData")
writeMat("../../Data/sim123_sparse_highSNR.mat",xmat_matlab_highSNR=xmat_matlab_highSNR)

