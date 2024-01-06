###################################################################
# R codes for real data application
# 05: Estimate mixing matrix for each subject in ABIDEI-NYU
# Note: This script was run on the group cluster server.
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(oro.nifti)
library(neurobase)
library(ciftiTools)
ciftiTools.setOption('wb_path', '/usr/local/workbench')
source("00_utils.R")

load("../Results/group_fastICA.RData")

subj_list = read.csv("../Data/ICA_subjlist_combine.csv")
table(subj_list$SITE_ID)


# Estimate M on ABIDEI-NYU subjects
cat("###########################################################################\n")
cat("################## Estimate M matrix on ABIDEI-NYU ########################\n")
cat("###########################################################################\n")

#npc=85
iter=5
my_list = subj_list[which(subj_list$SITE_ID=="ABIDEI-NYU"),]
setwd("/data/home4/risk_share/ImproveFConnASD/ABIDE/fmriprep_preprocessed/36p/ABIDEI-NYU/")
for (i in 1:dim(my_list)[1]) {
  file_list=list.files(paste0(my_list$SUB_ID[i],"/cifti_based/"))
  nii_index=grep("_9p.dtseries.nii",file_list)
  nii_name=file_list[nii_index]
  # read cifti
  cifti_image=read_cifti(cifti_fname = paste0(my_list$SUB_ID[i],"/cifti_based/",nii_name),brainstructures = "all")
  # prepare data 
  dat_left = cifti_image$data$cortex_left
  dat_right = cifti_image$data$cortex_right
  #dat_subcort = cifti_image$data$subcort
  dat_all = rbind(dat_left,dat_right)
  #dat_all = rbind(dat_all,dat_subcort)
  
  dat_all = scale(dat_all,center = T,scale = F)
  
  my_M = est.M.ols(my_fastICA$estS_sign,dat_all)
  
  save(my_M,file=paste0("/home/zwan873/Real_Data_Application/est_M_fastICA/ABIDEI-NYU/",my_list$SUB_ID[i],"_estM_subPC85.RData"))
  
  cat(my_list$SUB_ID[i]," finished!\n")
}


