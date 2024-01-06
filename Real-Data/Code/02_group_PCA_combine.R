###################################################################
# R codes for real data application
# 02: Concatenate all subject data together and extract group PCs
# Note: This script was run on our group cluster server.
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(ciftiTools)
ciftiTools.setOption('wb_path', '/usr/local/workbench')
library(irlba)
source("00_utils.R")

subj_list = read.csv("../Data/ICA_subjlist_combine.csv")

setwd("/home/zwan873/Real_Data_Application/sub_PC")

study_list = c("ABIDEI-KKI","ABIDEI-NYU","ABIDEII-KKI","ABIDEII-NYU_1")

dat_whole=NULL
npc=85
for (study in study_list) {
  cat("###########################################################################\n")
  cat("############################ Start ",study, "##############################\n")
  cat("###########################################################################\n")
  my_list = subj_list[which(subj_list$SITE_ID==study),]
  for (subj in my_list$SUB_ID) {
    load(paste0(study,"/",subj,"_PC",npc,"_cortex_9p_iter.RData"))
    dat_whole = cbind(dat_whole,PC)
    cat("The subject, ",subj," in study ",study, "finished!\n")
  }
  cat("###########################################################################\n")
  cat("############################## End ",study, "##############################\n")
  cat("###########################################################################\n")
}

cat("The whole data matrix is generated! The dimension is ",dim(dat_whole),"\n")
#save(dat_whole,file=paste0("dat_whole_subPC",npc,".RData"))

# perform group PCA
temp = whitener(X = dat_whole,n.comp = 30,irlba=T)
PC30 = temp$Z

save(PC30,file=paste0("/home/zwan873/Real_Data_Application/sub_PC/PC30_subPC",npc,".RData"))

