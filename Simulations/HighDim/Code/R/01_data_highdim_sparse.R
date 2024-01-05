###################################################################
# R codes for simulation setting 2 (High-Dimensional) in single-subject simulations
# 01: generate simulated high-dim data
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

# true_image=read_cifti("../../Data/group_sparseICA.nii",brainstructures = "all")
# # We picked three components:
# IC17: Dorsal Attention 
# view_cifti_surface(xifti = select_xifti(true_image,17),zlim = c(-2.43,2.82))
# IC6: Default
# view_cifti_surface(xifti = select_xifti(true_image,6),zlim = c(-2.43,2.82))
# IC23: Visual 
# view_cifti_surface(xifti = select_xifti(true_image,23),zlim = c(-2.43,2.82))

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


#######################################################################################################
#generate data for matlab
#split into 5 parts

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed1)
for (i in 1:20) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i]] = simData$X
  smat_list[[i]] = simData$S
  mmat_list[[i]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_sparse_1_20.RData")
writeMat("../../Data/highdim_sparse_1_20.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed2)
for (i in 21:40) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-20]] = simData$X
  smat_list[[i-20]] = simData$S
  mmat_list[[i-20]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_sparse_21_40.RData")
writeMat("../../Data/highdim_sparse_21_40.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed3)
for (i in 41:60) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-40]] = simData$X
  smat_list[[i-40]] = simData$S
  mmat_list[[i-40]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_sparse_41_60.RData")
writeMat("../../Data/highdim_sparse_41_60.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed4)
for (i in 61:80) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-60]] = simData$X
  smat_list[[i-60]] = simData$S
  mmat_list[[i-60]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_sparse_61_80.RData")
writeMat("../../Data/highdim_sparse_61_80.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)

xmat_list = list()
smat_list = list()
mmat_list = list()
xmat_matlab_lowSNR=c()
set.seed(seed5)
for (i in 81:100) {
  simData = SimFMRIreal(snr=my_SNR, noisyICA = my_noisyICA, nTR=my_nTR, phi=my_phi, var.inactive = var_level, components = my_components, true_data = dat_fastICA,template_xifti = template_image,template_time = sample_time)
  xmat_list[[i-80]] = simData$X
  smat_list[[i-80]] = simData$S
  mmat_list[[i-80]] = simData$Ms
  xmat_matlab_lowSNR = cbind(xmat_matlab_lowSNR,simData$X)
  cat(i," finished!\n")
}
save(xmat_list,mmat_list,smat_list,file = "../../Data/highdim_sparse_81_100.RData")
writeMat("../../Data/highdim_sparse_81_100.mat",xmat_matlab_lowSNR=xmat_matlab_lowSNR)
