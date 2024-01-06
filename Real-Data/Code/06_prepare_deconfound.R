###################################################################
# R codes for real data application
# 06: Prepare data for de-biased (deconfounded) estimation for group difference of functional connectivity
# Note: This script was run on a personal computer.
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

my_list = read.csv("../Data/ICA_subjlist_combine.csv")
table(my_list$DX_GROUP)

# create functional connectivity matrix
cor_all=array(dim = c(30,30,312))
for (i in 1:dim(my_list)[1]) {
  load(paste0("../Data/est_M_sparseICA/",my_list$SITE_ID[i],"/",my_list$SUB_ID[i],"_estM_subPC85.RData"))
  my_cor = cor(t(my_M))
  cor_all[,,i]=my_cor
  cat(my_list$SITE_ID[i],", ",my_list$SUB_ID[i]," finished!\n")
}

####################################################
# site harmonization
####################################################

# install combat r package
#library(devtools)
#install_github("jfortin1/neuroCombat_Rpackage")
#https://github.com/Jfortin1/neuroCombat_Rpackage/tree/2ccab2b42c59163717f012c70a54d20e8757b1c7

library(neuroCombat)

# receiving coil and slice thickness
receiving_coil = my_list$RECEIVING_COIL
slice_thickness = ifelse(my_list$SITE_ID %in% c("ABIDEI-NYU", "ABIDEII-NYU_1"), 4, 3)

# batch
batch = numeric(312)
batch[grepl("KKI", my_list$SITE_ID)] = ifelse(receiving_coil[grepl("KKI", my_list$SITE_ID)] == "8 channel", 1, 2)
batch[grepl("NYU", my_list$SITE_ID)] = ifelse(slice_thickness[grepl("NYU", my_list$SITE_ID)] == 4, 3, 4)
# batch only has three categories 

# # create a corresponding column
my_list$site_feature = NA
my_list$site_feature[batch == 1] = "KKI-8 channel"
my_list$site_feature[batch == 2] = "KKI-32 channel"
my_list$site_feature[batch == 3] = "NYU-4 thickness"
# # write the csv
# write.csv(abide, file = "./data/delta_all_pheno.csv")

# KKI
# ABIDEI_KKI     8 channel, 3 thickness
# ABIDEII_KKI_1  8 channel, 3 thickness
# ABIDEII_KKI_2  32 channel, 3 thickness 
# NYU
# ABIDEI_NYU     8 channel, 4 thickness
# ABIDEII_NYU_1  8 channel, 4 thickness

table(my_list$site_feature)
# 
# KKI-32 channel   KKI-8 channel NYU-4 thickness 
#             45             154             113 

# input:
# dat (p*n). p is the number of features, and n is the number of participants. Here p = 435, n = 312
# batch (a vector of length n), indicates the site
# mod: a matrix containing biological covariates
#      to ensure that biological variability is preserved in the harmonization process

#####################
# create dat by z-transforming cor mat
Y = matrix(nrow = 435,ncol = 312)
for (i in 1:312) {
  my_cor = cor_all[,,i]
  #Y[,i] = FisherZ(my_cor[lower.tri(my_cor)])
  Y[,i] = atanh(my_cor[lower.tri(my_cor)])
}

# mod
# the Np. 157 subject has missing FIQ, impute by nearby values
sum(is.na(my_list$FIQ))
my_list$FIQ[157]
my_list$FIQ[157]=124

mod = model.matrix(~ my_list$DX_GROUP + my_list$Mean_FD + my_list$AGE_AT_SCAN + my_list$SEX + my_list$HANDEDNESS_LR +
                     my_list$Stimulant + my_list$NonStimulant + my_list$ADOS_combine + my_list$FIQ + my_list$Propfd02 + my_list$PropRMSD025) 

# harmonization with parametric adjustment 
# parametric priors are used in the EB estimation
data.harmonized = neuroCombat(dat = Y, batch = batch, mod = mod)
#save(data.harmonized, file = "corvec_harmonized.RData")

# results
Y_harmonized = data.harmonized$dat.combat

# save result
save(Y_harmonized, file = "../Data/Y_harmonized_sparseICA.RData")

####################################################
# adjusted residuals
####################################################

covariates_adj = my_list[,c(1,10,12,14,16,17,18,24,25,33,34,35)]

# center covariates
covariates_adj$Mean_FD=as.vector(scale(covariates_adj$Mean_FD,scale = F))
covariates_adj$SEX[which(covariates_adj$SEX==2)]=-1
covariates_adj$AGE_AT_SCAN=as.vector(scale(covariates_adj$AGE_AT_SCAN,scale = F))
covariates_adj$FIQ=as.vector(scale(covariates_adj$FIQ,scale = F))
covariates_adj$Propfd02=as.vector(scale(covariates_adj$Propfd02,scale = F))
covariates_adj$PropRMSD025=as.vector(scale(covariates_adj$PropRMSD025,scale = F))
covariates_adj$HANDEDNESS_LR=relevel(as.factor(covariates_adj$HANDEDNESS_LR),ref = "R")

# calculate adjusted residuals
res_adj = matrix(nrow = 435,ncol=312)
for (i in 1:435) {
  covariates_adj$y = Y_harmonized[i,]
  my_lm = lm(y ~ DX_GROUP + Mean_FD + Propfd02 + PropRMSD025 + SEX + AGE_AT_SCAN + HANDEDNESS_LR,data = covariates_adj)
  res_adj[i,]= my_lm$residuals + my_lm$coefficients[1] + my_lm$coefficients[2]*covariates_adj$DX_GROUP
}

dim(res_adj)

save(res_adj,file = "../Results/adjusted_residuals_sparseICA.RData")


####################################################
# naive z-statistics
####################################################

res_adj_ASD = res_adj[,which(my_list$DX_GROUP==1)]
res_adj_TD = res_adj[,which(my_list$DX_GROUP==2)]

z_stat = matrix(nrow = 435,ncol = 3)
for (i in 1:435) {
  asd = res_adj_ASD[i,]
  td = res_adj_TD[i,]
  
  diff = mean(asd) - mean(td)
  diff_se = sqrt(var(asd)/length(asd)+var(td)/length(td))
  z = diff/diff_se
  z_p = 2*pnorm(abs(z),lower.tail = F)
  
  z_stat[i,]=c(z,z_p,NA)
}

hist(z_stat[,1])
hist(z_stat[,2])

# FDR control
z_stat[,3] = p.adjust(z_stat[,2],method = "BH")

length(which(z_stat[,2]<0.05/435))
length(which(z_stat[,3]<0.2))
length(which(z_stat[,3]<0.05))

save(z_stat,res_adj_ASD,res_adj_TD,res_adj,file = "../Results/naive_z_stat_sparseICA.RData")

####################################################
# Make figures for naive z-stat results
####################################################

library(corrplot)
library(ppcor)

# Bonferroni
full_cor = diag(30)
full_cor[lower.tri(full_cor)] = z_stat[,1]
full_cor[upper.tri(full_cor)] = t(full_cor)[upper.tri(full_cor)]
colnames(full_cor)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                     "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                     "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")
rownames(full_cor)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                     "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                     "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")

full_cor_p = diag(x=1,30)
full_cor_p[lower.tri(full_cor_p)] = z_stat[,2]
full_cor_p[upper.tri(full_cor_p)] = t(full_cor_p)[upper.tri(full_cor_p)]

colnames(full_cor_p)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                       "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                       "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")
rownames(full_cor_p)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                       "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                       "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")

corrplot(full_cor,method = "color",tl.col="black",is.corr = F,p.mat = full_cor_p,insig = "label_sig",sig.level = 0.05/435)


# FDR=0.2
full_cor = diag(30)
full_cor[lower.tri(full_cor)] = z_stat[,1]
full_cor[upper.tri(full_cor)] = t(full_cor)[upper.tri(full_cor)]
colnames(full_cor)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                     "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                     "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")
rownames(full_cor)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                     "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                     "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")

full_cor_p = diag(x=1,30)
full_cor_p[lower.tri(full_cor_p)] = z_stat[,3]
full_cor_p[upper.tri(full_cor_p)] = t(full_cor_p)[upper.tri(full_cor_p)]

colnames(full_cor_p)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                       "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                       "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")
rownames(full_cor_p)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                       "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                       "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")

corrplot(full_cor,method = "color",tl.col="black",is.corr = F,p.mat = full_cor_p,insig = "label_sig",sig.level = 0.2)



####################################################
# Prepare data for de-biased (deconfounded) group difference - DRTMLE/AIPWE
####################################################

library(SuperLearner)
library(drtmle)

# clean data
dat_used = my_list[,c(1,10,12,14,16,17,18,24,25,33,34,35)]
dat_used$delta = 1
dat_used[,14:448] = t(res_adj)

load("../Data/abide_Y_9p_all_419.Rda")
abide$cor_mat=NULL

# 396 subjects in total

# 396-312 = 84 unused data
dat_unused = abide[-which(abide$SUB_ID%in%dat_used$SUB_ID),]
dat_unused = dat_unused[,c(1,10,12,14,16,17,18,24,25,33,34,35)]
dat_unused$delta = 0
dat_unused[,14:448]=NA

dat_all = rbind(dat_used,dat_unused)

save(dat_all,dat_unused,dat_used, file = "../Data/dat_for_deconfound_sparseICA.RData")

