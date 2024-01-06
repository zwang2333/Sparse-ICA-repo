###################################################################
# R codes for real data application
# 07: Build propensity and outcome model, estimate de-biased group difference for a single seed
# Note: This script was run on the group cluster with other bash files
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

local=FALSE
# Change this for cluster or local:

if(local){
  setwd('D:/Files/Real_Data_Application/Deconfound')
  save.input.data = FALSE
  getOption("mc.cores")
  options(mc.cores=8)
  seed=1
} else {
  setwd('/home/zwan873/Real_Data_Application/Deconfound')
  save.input.data = FALSE
  options(mc.cores=1)
  seed = seedID
}
.libPaths('/home/zwan873/R/x86_64-redhat-linux-gnu-library/4.1')
.libPaths()

getOption("mc.cores")
set.seed(seed, "L'Ecuyer-CMRG")

tic = proc.time()
library(drtmle)
library(SuperLearner)
library(earth)
library(quadprog)
library(gam)
library(nloptr)
library(future)
library(future.apply)
library(xgboost)
library(ranger)
library(visdat)
library(ggplot2)
library(gridExtra)
#library(tidyr)
library(e1071)
library(glmnet)
library(readxl)
library(ROCit)

####################################################
# Deconfounded group difference - AIPWE
####################################################

load("../Data/dat_for_deconfound_sparseICA.RData")

####################################################
# fit propensity model
####################################################

gn.variables=c("delta","SEX","AGE_AT_SCAN","FIQ","HANDEDNESS_LR",
               "DX_GROUP","Stimulant","NonStimulant","ADOS_combine")

# extract design matrix for outcome model fitting
temp.data = dat_all[gn.variables]
gn.xmat = data.frame(model.matrix(delta~.,data=temp.data))[,-1]

# fit with super learner
my.SL.libs.gn = c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.step","SL.step.interaction","SL.xgboost","SL.mean")
gn.est=mcSuperLearner(Y = dat_all$delta, X = gn.xmat, family=binomial(),SL.library = my.SL.libs.gn, cvControl = list(V = 10), method='method.CC_nloglik')

gn.est.predict = as.vector(gn.est$SL.predict)
summary(rocit(score=gn.est.predict,class=dat_all$delta,method='nonparametric'))

save(gn.est.predict,file=paste0("../Data/deconfound_sparseICA/propensity_model/gn_est_predict_",seed,".RData"))

####################################################
# fit outcome model
####################################################

Qn.variables=c("delta","SEX","AGE_AT_SCAN","FIQ","HANDEDNESS_LR",
               "DX_GROUP","Stimulant","NonStimulant","ADOS_combine")

# extract design matrix for outcome model fitting
temp.data = dat_used[Qn.variables]
Qn.xmat.fit = data.frame(model.matrix(delta~.,data=temp.data))[,-1]

# design matrix for outcome prediction
temp.data = dat_all[Qn.variables]
Qn.xmat.predict = data.frame(model.matrix(delta~.,data=temp.data))[,-1]

# Separate ASD and TD datasets are necessary to obtain the DRTMLE estimates:
temp.data = Qn.xmat.predict[Qn.xmat.predict$DX_GROUP==1,]
Qn.xmat.predict.asd = data.frame(model.matrix(numeric(nrow(temp.data))~.,data=temp.data)[,-1])

temp.data = Qn.xmat.predict[Qn.xmat.predict$DX_GROUP==2,]
Qn.xmat.predict.td = data.frame(model.matrix(numeric(nrow(temp.data))~.,data=temp.data))[,-1]

my.SL.libs.Qbar= c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.ridge","SL.step","SL.step.interaction","SL.svm","SL.xgboost","SL.mean")

# fit outcome model for each edge
Qbar.SL.asd_mat = matrix(nrow = 435,ncol = 144)
Qbar.SL.td_mat = matrix(nrow = 435,ncol = 252)
for (i in 1:435) {
  outcome.SL = mcSuperLearner(Y = dat_used[,13+i],X=Qn.xmat.fit,family=gaussian(), SL.library = my.SL.libs.Qbar,cvControl = list(V = 10), method = drtmle:::tmp_method.CC_LS)
  Qbar.SL.asd_mat[i,] = predict(outcome.SL, newdata = Qn.xmat.predict.asd)[[1]]
  Qbar.SL.td_mat[i,] = predict(outcome.SL, newdata = Qn.xmat.predict.td)[[1]]
  cat(i," th edge finished!\n")
}

save(Qbar.SL.asd_mat,Qbar.SL.td_mat,file=paste0("../Data/deconfound_sparseICA/outcome_model/outcome_predict_",seed,".RData"))

####################################################
# Final AIPWE estimation
####################################################

#### AIPTW function for one Y vector input.
AIPTW <- function(Y,Delta,Qn,gn){ #Y has NA when \delta = 0.
  
  Y[Delta==0]=0 # transfer the NA in Y to 0.
  n=length(Y)
  mean_fconn_AIP <- mean(Delta/gn*Y+ (1-Delta/gn)*Qn) #\psi
  Zi=Delta/gn*(Y-Qn)+Qn-mean_fconn_AIP  #Zi
  Zbar=sum(Zi)/n  #Zbar=1/n*sum(Zi)
  var_fconn_AIP <- sum((Zi-Zbar)^2)/(n-1)/n #sum((Zi-Zbar)^2)/(n-1)/n
  AIP=list(mean=mean_fconn_AIP,var=var_fconn_AIP)
  return(AIP)
}

z_aiptw_stat = matrix(nrow = 435,ncol = 6)
for (i in 1:435) {
  my_col = i+13
  
  Qbar.SL.asd = Qbar.SL.asd_mat[i,]
  Qbar.SL.td = Qbar.SL.td_mat[i,]
  
  mean_fconn_asd <- AIPTW(Y = dat_all[,my_col][which(dat_all$DX_GROUP==1)],Delta = dat_all$delta[which(dat_all$DX_GROUP==1)],Qn = Qbar.SL.asd,gn = gn.est.predict[which(dat_all$DX_GROUP==1)])
  mean_fconn_td <- AIPTW(Y = dat_all[,my_col][which(dat_all$DX_GROUP==2)],Delta = dat_all$delta[which(dat_all$DX_GROUP==2)],Qn = Qbar.SL.asd,gn = gn.est.predict[which(dat_all$DX_GROUP==2)])
    
  diff=mean_fconn_asd$mean-mean_fconn_td$mean
  se.diff=sqrt(mean_fconn_asd$var+mean_fconn_td$var)
  z.diff=diff/se.diff
    
  z_p = 2*pnorm(abs(z.diff),lower.tail = F)
    
  z_aiptw_stat[i,]=c(mean_fconn_asd$mean,mean_fconn_td$mean,z.diff,z_p,diff,se.diff)
    
  cat(i,"th AIPWE connectivity finished!\n")
  
}

save(z_aiptw_stat,file = paste0("../Data/deconfound_sparseICA/AIPWE/AIPWE_ASD_TD_z_stat_",seed,".RData"))

proc.time()-tic
