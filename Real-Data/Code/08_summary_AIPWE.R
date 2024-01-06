###################################################################
# R codes for real data application
# 08: summarize AIPWE estimation results and make figures
# Note: This script was run on a personal computer
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

###################################################################
# extract AIPWE results of 400 different seeds
###################################################################
library(corrplot)
library(ppcor)

# load DRTMLE and AIPTW data
z_drtmle_stat_all = matrix(nrow = 435,ncol = 400)
z_aiptw_stat_all = matrix(nrow = 435,ncol = 400)

for (i in 1:400) {
  load(paste0("../Data/deconfound_sparseICA/AIPWE/AIPWE_ASD_TD_z_stat_",i,".RData"))
  
  z_drtmle_stat_all[,i]=z_drtmle_stat[,3]
  z_aiptw_stat_all[,i]=z_aiptw_stat[,3]
}

###################################################################
# make a plot about the relationship between naive and DRTMLE/AIPTW z stat
###################################################################
load("../Data/naive_z_stat_sparseICA.RData")

dat_long = data.frame(x=numeric(174000),y1=numeric(174000),y2=numeric(174000))
my_x = c()
for (i in 1:435) {
  temp_x = rep(z_stat[i,1],400)
  my_x = c(my_x,temp_x)
}
dat_long$x=my_x
my_y2 = c()
for (i in 1:435) {
  my_y2 = c(my_y2,z_aiptw_stat_all[i,])
}
dat_long$y2=my_y2

pdf("../Figures/AIPWE_vs_naive_sparseICA.pdf",width = 10, height = 6)
plot(dat_long$x,dat_long$y2,xlab = "Naive Z",ylab = "AIPWE Z",main = "AIPWE vs Naive")
abline(a=0,b=1,col="red")
dev.off()


###################################################################
# calculate mean AIPWE
###################################################################
z_aiptw_stat_mean = apply(z_aiptw_stat_all,1,mean)

# add p-value
z_aiptw = matrix(nrow = 435,ncol = 3)

z_aiptw[,1]=z_aiptw_stat_mean
z_aiptw[,2]=2*pnorm(abs(z_aiptw_stat_mean),lower.tail = F)
z_aiptw[,3] = p.adjust(z_aiptw[,2],method = "BH")

length(which(z_stat[,3]<0.05))
length(which(z_aiptw[,3]<0.05))
cor(z_aiptw[,1],z_stat[,1])


###################################################################
# Make figures for AIPWE-adjusted results
###################################################################
full_cor = diag(30)
full_cor[lower.tri(full_cor)] = z_aiptw[,1]
full_cor[upper.tri(full_cor)] = t(full_cor)[upper.tri(full_cor)]
colnames(full_cor)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                     "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                     "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")
rownames(full_cor)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                     "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                     "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")

full_cor_p = diag(x=1,30)
full_cor_p[lower.tri(full_cor_p)] = z_aiptw[,3]
full_cor_p[upper.tri(full_cor_p)] = t(full_cor_p)[upper.tri(full_cor_p)]

colnames(full_cor_p)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                       "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                       "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")
rownames(full_cor_p)=c("IC1","IC2","IC3","IC4","IC5","IC6","IC7","IC8","IC9","IC10",
                       "IC11","IC12","IC13","IC14","IC15","IC16","IC17","IC18","IC19","IC20",
                       "IC21","IC22","IC23","IC24","IC25","IC26","IC27","IC28","IC29","IC30")

pdf("../Figures/AIPWE_adjusted_sparseICA.pdf",width = 7,height = 7)
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F)
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F,order = "hclust")
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F,order = "hclust",p.mat = full_cor_p,insig = "label_sig",sig.level = 0.05/435)
corrplot(full_cor,method = "color",tl.col="black",is.corr = F,p.mat = full_cor_p,insig = "label_sig",sig.level = 0.05)
dev.off()

which(z_aiptw[,3]<0.05)
z_aiptw[which(z_aiptw[,3]<0.05),]

# naive
# FDR=0.05
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

pdf("../Figures/Naive_sparseICA.pdf",width = 7,height = 7)
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F)
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F,order = "hclust")
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F,order = "hclust",hclust.method = "ward.D")
#corrplot(full_cor,method = "color",tl.col="black",is.corr = F,order = "hclust",p.mat = full_cor_p,insig = "label_sig",sig.level = 0.05/435)
corrplot(full_cor,method = "color",tl.col="black",is.corr = F,p.mat = full_cor_p,insig = "label_sig",sig.level = 0.05)
dev.off()

save(z_stat,z_drtmle,z_aiptw,z_drtmle_stat_all,z_aiptw_stat_all,file = "../Results/deconfounded_results.RData")


###################################################################
# make group comparison plot for z-transformed correlation
###################################################################

library(ggplot2)

# read outcome model
out_asd=0
out_td=0

for (i in 1:400) {
  load(paste0("../Data/deconfound_sparseICA/outcome_model/outcome_predict_",i,".RData"))
  out_asd = out_asd + Qbar.SL.asd_mat[118,]
  out_td = out_td + Qbar.SL.td_mat[118,]
}
out_asd = out_asd/400
out_td = out_td/400

dat_jitter = data.frame(group = c(rep("ASD",144),rep("TD",252)),z = c(out_asd,out_td))

ggplot(dat_jitter,aes(group,z))+geom_jitter()

plot(out_asd)


# load AIPWE data
asd_aiptw_all = matrix(nrow = 435,ncol = 400)
td_aiptw_all = matrix(nrow = 435,ncol = 400)

for (i in 1:400) {
  load(paste0("../Data/deconfound_sparseICA/AIPWE/AIPWE_ASD_TD_z_stat_",i,".RData"))
  asd_aiptw_all[,i]=z_aiptw_stat[,1]
  td_aiptw_all[,i]=z_aiptw_stat[,2]
}

dat_jitter = data.frame(group = c(rep("ASD",99),rep("TD",213)),z = c(res_adj_ASD[118,],res_adj_TD[118,]))
pdf("../Figures/IC18_20_Jitter_sparseICA.pdf",width = 4,height = 4)
ggplot(dat_jitter,aes(group,z))+
  geom_jitter(color="#666666",width = 0.2) + 
  labs(x = "Primary Diagnosis", y = "IC18 - IC20") +
  geom_point(aes(x="ASD",y=mean(asd_aiptw_all[118,])),colour="red",size=3,shape=15)+
  geom_point(aes(x="TD",y=mean(td_aiptw_all[118,])),colour="red",size=3,shape=15)
dev.off()

dat_jitter = data.frame(group = c(rep("ASD",99),rep("TD",213)),z = c(res_adj_ASD[359,],res_adj_TD[359,]))
pdf("../Figures/IC5_13_Jitter_sparseICA.pdf",width = 4,height = 4)
ggplot(dat_jitter,aes(group,z))+
  geom_jitter(color="#666666",width = 0.2) + 
  labs(x = "Primary Diagnosis", y = "IC5 - IC13") +
  geom_point(aes(x="ASD",y=mean(asd_aiptw_all[359,])),colour="red",size=3,shape=15)+
  geom_point(aes(x="TD",y=mean(td_aiptw_all[359,])),colour="red",size=3,shape=15)
dev.off()


