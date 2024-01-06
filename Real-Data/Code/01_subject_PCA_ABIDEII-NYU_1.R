###################################################################
# R codes for real data application
# 01: Perform subject-level PCA on ABIDEII-NYU data
# Note: This script was run on our group cluster server.
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(ciftiTools)
ciftiTools.setOption('wb_path', '/usr/local/workbench')
library(irlba)

# Read the list of participants who passed QC procedure
subj_list = read.csv("../Data/ICA_subjlist_combine.csv")
table(subj_list$SITE_ID)

# Perform subject-level PCA on ABIDEII-NYU_1 subjects
cat("###########################################################################\n")
cat("#################### Start PCA on ABIDEII-NYU_1 ##############################\n")
cat("###########################################################################\n")

# number of PCs to be retained for each subject
npc = 85
# number of iterative standardizations
iter = 5
my_list = subj_list[which(subj_list$SITE_ID=="ABIDEII-NYU_1"),]
setwd("/data/home4/risk_share/ImproveFConnASD/ABIDE/fmriprep_preprocessed/36p/ABIDEII-NYU_1/")
for (i in 1:dim(my_list)[1]) {
  file_list=list.files(paste0(my_list$SUB_ID[i],"/cifti_based/"))
  nii_index=grep("_9p.dtseries.nii",file_list)
  nii_name=file_list[nii_index]
  # read cifti files
  cifti_image=read_cifti(cifti_fname = paste0(my_list$SUB_ID[i],"/cifti_based/",nii_name),brainstructures = "all")
  # prepare data 
  dat_left = cifti_image$data$cortex_left
  dat_right = cifti_image$data$cortex_right
  #dat_subcort = cifti_image$data$subcort
  dat_all = rbind(dat_left,dat_right)
  #dat_all = rbind(dat_all,dat_subcort)
  
  # iterative standardization
  for (j in 1:iter) {
    # standardization at each voxel
    dat_all[which(dat_all[,1]!=0),] = t(scale(t(dat_all[which(dat_all[,1]!=0),])))
    # standardization at each time point
    dat_all = scale(dat_all)
  }
  
  # perform PCA
  subj_PCA=prcomp_irlba(dat_all,npc)
  PC=subj_PCA$x
  dimnames(PC)=NULL
  # cifti_image$data$cortex_left=PC30[1:29696,]
  # cifti_image$data$cortex_right=PC30[29697:59412,]
  # cifti_image$data$subcort=PC30[59413:91282,]

  #write_cifti(cifti_image,paste0("/home/zwan873/group_ICA/sub_PC/ABIDEI/",file_list$V1[i],"_PC30"))
  save(PC,file=paste0("/home/zwan873/Real_Data_Application/sub_PC/ABIDEII-NYU_1/",my_list$SUB_ID[i],"_PC",npc,"_cortex_9p_iter.RData"))
  
  cat(my_list$SUB_ID[i]," finished!\n")
}



