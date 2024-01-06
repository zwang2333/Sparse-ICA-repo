###################################################################
# R codes for real data application
# 04: Run group ICA using Sparse ICA and Fast ICA
# Note: This script was run on a personal computer.
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(irlba)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench') 

load("../Data/PC30_subPC85.RData")
load("../Results/nu_selection_groupPC_BIC.RData")


# run group Sparse ICA
best_nu=nu_selection$best_nu
my_sparseICA_all = sparseICA(xData = PC30,n.comp = 30,nu = best_nu,method = "C",
                              whiten = "none",restarts = 40,irlba = T,converge_plot = T)

# sign change
a=signchange(t(my_sparseICA_all$estS))
my_sparseICA_all$estS_sign=t(a$S)

save(my_sparseICA_all,file="../Results/group_sparseICA.RData")

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")
cifti_image$data$cortex_left=my_sparseICA_all1$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_sparseICA_all1$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = 30)
write_cifti(cifti_image,"../Results/group_sparseICA")


# Run group Fast ICA
my_fastICA = fastICA_restarts(PC30,n.comp = 30,restarts = 40,method = "C",maxit = 500,tol = 1e-6,verbose = T)

# match with sparse ICA
my_fastICA$estS_sign=matchICA(my_fastICA$S,my_sparseICA_all$estS_sign)

save(my_fastICA,file="../Results/group_fastICA.RData")

cifti_image=read_cifti("../Data/cifti_template.dtseries.nii",brainstructures = "all")
cifti_image$data$cortex_left=my_fastICA2$estS_sign[1:29696,]
cifti_image$data$cortex_right=my_fastICA2$estS_sign[29697:59412,]
cifti_image$data$subcort=matrix(0,nrow = 31870,ncol = 30)
write_cifti(cifti_image,"../Results/group_fastICA")

