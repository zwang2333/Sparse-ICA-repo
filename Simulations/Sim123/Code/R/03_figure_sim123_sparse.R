###################################################################
# R codes for simulation setting 1 in single-subject simulations
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(SparseICA)
library(steadyICA)
library(fastICA)
library(irlba)
library(R.matlab)
library(ggplot2)
library(ggpubr)

source("00_utils.R")

###################################################################
# General simulations settings
###################################################################

# The three given SNR levels
my_SNR = c(0.4,1.5,3)
# Number of replications
rep = 100
# Variance level of inactive points in components
var_level = 0
# Overall seed
seed = 2023
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

load("../../Results/02_PRMSE_sim123_sparse.RData")

###################################################################
# Figures for single and multiple restarts
###################################################################

# For S estimates
S_PMSE = data.frame(PMSE=c(S_PMSE1$sparse,S_PMSE1$fast,S_PMSE1$infomax,S_PMSE1$sparse40,S_PMSE1$fast40,S_PMSE1$infomax40,S_PMSE1$ebm,S_PMSE1$sfi,
                           S_PMSE2$sparse,S_PMSE2$fast,S_PMSE2$infomax,S_PMSE2$sparse40,S_PMSE2$fast40,S_PMSE2$infomax40,S_PMSE2$ebm,S_PMSE2$sfi,
                           S_PMSE3$sparse,S_PMSE3$fast,S_PMSE3$infomax,S_PMSE3$sparse40,S_PMSE3$fast40,S_PMSE3$infomax40,S_PMSE3$ebm,S_PMSE3$sfi),
                    Method=c(rep("Sparse ICA-1",length(S_PMSE1$sparse)),rep("Fast ICA-1",length(S_PMSE1$fast)),rep("Infomax ICA-1",length(S_PMSE1$infomax)),rep("Sparse ICA-40",length(S_PMSE1$sparse)),rep("Fast ICA-40",length(S_PMSE1$fast)),rep("Infomax ICA-40",length(S_PMSE1$infomax)),rep("SICA-EBM",length(S_PMSE1$ebm)),rep("Sparse Fast ICA",length(S_PMSE1$sfi)),
                             rep("Sparse ICA-1",length(S_PMSE2$sparse)),rep("Fast ICA-1",length(S_PMSE2$fast)),rep("Infomax ICA-1",length(S_PMSE2$infomax)),rep("Sparse ICA-40",length(S_PMSE2$sparse)),rep("Fast ICA-40",length(S_PMSE2$fast)),rep("Infomax ICA-40",length(S_PMSE2$infomax)),rep("SICA-EBM",length(S_PMSE2$ebm)),rep("Sparse Fast ICA",length(S_PMSE2$sfi)),
                             rep("Sparse ICA-1",length(S_PMSE3$sparse)),rep("Fast ICA-1",length(S_PMSE3$fast)),rep("Infomax ICA-1",length(S_PMSE3$infomax)),rep("Sparse ICA-40",length(S_PMSE3$sparse)),rep("Fast ICA-40",length(S_PMSE3$fast)),rep("Infomax ICA-40",length(S_PMSE3$infomax)),rep("SICA-EBM",length(S_PMSE3$ebm)),rep("Sparse Fast ICA",length(S_PMSE3$sfi))
                             ),
                    SNR=c(rep("Low",8*length(S_PMSE1$sparse)),
                          rep("Medium",8*length(S_PMSE1$sparse)),
                          rep("High",8*length(S_PMSE1$sparse)))
                    )


p1 = ggplot(S_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) +
geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Independent Components")

# For M estimates
M_PMSE = data.frame(PMSE=c(M_PMSE1$sparse,M_PMSE1$fast,M_PMSE1$infomax,M_PMSE1$sparse40,M_PMSE1$fast40,M_PMSE1$infomax40,M_PMSE1$ebm,M_PMSE1$sfi,
                           M_PMSE2$sparse,M_PMSE2$fast,M_PMSE2$infomax,M_PMSE2$sparse40,M_PMSE2$fast40,M_PMSE2$infomax40,M_PMSE2$ebm,M_PMSE2$sfi,
                           M_PMSE3$sparse,M_PMSE3$fast,M_PMSE3$infomax,M_PMSE3$sparse40,M_PMSE3$fast40,M_PMSE3$infomax40,M_PMSE3$ebm,M_PMSE3$sfi),
                    Method=c(rep("Sparse ICA",length(M_PMSE1$sparse)),rep("Fast ICA",length(M_PMSE1$fast)),rep("Infomax ICA",length(M_PMSE1$infomax)),rep("Sparse ICA-40",length(M_PMSE1$sparse)),rep("Fast ICA-40",length(M_PMSE1$fast)),rep("Infomax ICA-40",length(M_PMSE1$infomax)),rep("SICA-EBM",length(M_PMSE1$ebm)),rep("Sparse Fast ICA",length(M_PMSE1$sfi)),
                             rep("Sparse ICA",length(M_PMSE2$sparse)),rep("Fast ICA",length(M_PMSE2$fast)),rep("Infomax ICA",length(M_PMSE2$infomax)),rep("Sparse ICA-40",length(M_PMSE2$sparse)),rep("Fast ICA-40",length(M_PMSE2$fast)),rep("Infomax ICA-40",length(M_PMSE2$infomax)),rep("SICA-EBM",length(M_PMSE1$ebm)),rep("Sparse Fast ICA",length(M_PMSE1$sfi)),
                             rep("Sparse ICA",length(M_PMSE3$sparse)),rep("Fast ICA",length(M_PMSE3$fast)),rep("Infomax ICA",length(M_PMSE3$infomax)),rep("Sparse ICA-40",length(M_PMSE3$sparse)),rep("Fast ICA-40",length(M_PMSE3$fast)),rep("Infomax ICA-40",length(M_PMSE3$infomax)),rep("SICA-EBM",length(M_PMSE1$ebm)),rep("Sparse Fast ICA",length(M_PMSE1$sfi))
                    ),
                    SNR=c(rep("Low",8*length(M_PMSE1$sparse)),
                          rep("Medium",8*length(M_PMSE1$sparse)),
                          rep("High",8*length(M_PMSE1$sparse)))
)

p2 = ggplot(M_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) +
 geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Mixing Matrix")


# For computation time
TIME = data.frame(TIME=c(TIME1$sparse,TIME1$fast,TIME1$infomax,TIME1$sparse40,TIME1$fast40,TIME1$infomax40,TIME1$ebm,TIME1$sfi,
                         TIME2$sparse,TIME2$fast,TIME2$infomax,TIME2$sparse40,TIME2$fast40,TIME2$infomax40,TIME2$ebm,TIME2$sfi,
                         TIME3$sparse,TIME3$fast,TIME3$infomax,TIME3$sparse40,TIME3$fast40,TIME3$infomax40,TIME3$ebm,TIME3$sfi),
                    Method=c(rep("Sparse ICA",length(TIME1$sparse)),rep("Fast ICA",length(TIME1$fast)),rep("Infomax ICA",length(TIME1$infomax)),rep("Sparse ICA-40",length(TIME1$sparse)),rep("Fast ICA-40",length(TIME1$fast)),rep("Infomax ICA-40",length(TIME1$infomax)),rep("SICA-EBM",length(TIME1$ebm)),rep("Sparse Fast ICA",length(TIME1$sfi)),
                             rep("Sparse ICA",length(TIME2$sparse)),rep("Fast ICA",length(TIME2$fast)),rep("Infomax ICA",length(TIME2$infomax)),rep("Sparse ICA-40",length(TIME2$sparse)),rep("Fast ICA-40",length(TIME2$fast)),rep("Infomax ICA-40",length(TIME2$infomax)),rep("SICA-EBM",length(TIME2$ebm)),rep("Sparse Fast ICA",length(TIME2$sfi)),
                             rep("Sparse ICA",length(TIME3$sparse)),rep("Fast ICA",length(TIME3$fast)),rep("Infomax ICA",length(TIME3$infomax)),rep("Sparse ICA-40",length(TIME3$sparse)),rep("Fast ICA-40",length(TIME3$fast)),rep("Infomax ICA-40",length(TIME3$infomax)),rep("SICA-EBM",length(TIME3$ebm)),rep("Sparse Fast ICA",length(TIME3$sfi))),
                    SNR=c(rep("Low",8*length(TIME1$sparse)),
                          rep("Medium",8*length(TIME2$sparse)),
                          rep("High",8*length(TIME3$sparse)))
)

p3 = ggplot(TIME, aes(factor(SNR,levels = c("Low","Medium","High")), TIME)) +
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="TIME",title = "Computation Time")#+ylim(c(0,0.25))


ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../../Figures/03_sim123_PRMSE_single_multi.pdf",width = 10, height = 6)


###################################################################
# Figures for single restart
###################################################################
# For S estimates
S_PMSE = data.frame(PMSE=c(S_PMSE1$sparse,S_PMSE1$fast,S_PMSE1$infomax,S_PMSE1$ebm,S_PMSE1$sfi,
                           S_PMSE2$sparse,S_PMSE2$fast,S_PMSE2$infomax,S_PMSE2$ebm,S_PMSE2$sfi,
                           S_PMSE3$sparse,S_PMSE3$fast,S_PMSE3$infomax,S_PMSE3$ebm,S_PMSE3$sfi),
                    Method=c(rep("Sparse ICA",length(S_PMSE1$sparse)),rep("Fast ICA",length(S_PMSE1$fast)),rep("Infomax ICA",length(S_PMSE1$infomax)),rep("SICA-EBM",length(S_PMSE1$ebm)),rep("Sparse Fast ICA",length(S_PMSE1$sfi)),
                             rep("Sparse ICA",length(S_PMSE2$sparse)),rep("Fast ICA",length(S_PMSE2$fast)),rep("Infomax ICA",length(S_PMSE2$infomax)),rep("SICA-EBM",length(S_PMSE2$ebm)),rep("Sparse Fast ICA",length(S_PMSE2$sfi)),
                             rep("Sparse ICA",length(S_PMSE3$sparse)),rep("Fast ICA",length(S_PMSE3$fast)),rep("Infomax ICA",length(S_PMSE3$infomax)),rep("SICA-EBM",length(S_PMSE3$ebm)),rep("Sparse Fast ICA",length(S_PMSE3$sfi))
                    ),
                    SNR=c(rep("Low",5*length(S_PMSE1$sparse)),
                          rep("Medium",5*length(S_PMSE1$sparse)),
                          rep("High",5*length(S_PMSE1$sparse)))
)


p1 = ggplot(S_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) +
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Independent Components")
#theme(legend.position  = "none")


# For M estimates
M_PMSE = data.frame(PMSE=c(M_PMSE1$sparse,M_PMSE1$fast,M_PMSE1$infomax,M_PMSE1$ebm,M_PMSE1$sfi,
                           M_PMSE2$sparse,M_PMSE2$fast,M_PMSE2$infomax,M_PMSE2$ebm,M_PMSE2$sfi,
                           M_PMSE3$sparse,M_PMSE3$fast,M_PMSE3$infomax,M_PMSE3$ebm,M_PMSE3$sfi),
                    Method=c(rep("Sparse ICA",length(M_PMSE1$sparse)),rep("Fast ICA",length(M_PMSE1$fast)),rep("Infomax ICA",length(M_PMSE1$infomax)),rep("SICA-EBM",length(M_PMSE1$ebm)),rep("Sparse Fast ICA",length(M_PMSE1$sfi)),
                             rep("Sparse ICA",length(M_PMSE2$sparse)),rep("Fast ICA",length(M_PMSE2$fast)),rep("Infomax ICA",length(M_PMSE2$infomax)),rep("SICA-EBM",length(M_PMSE1$ebm)),rep("Sparse Fast ICA",length(M_PMSE1$sfi)),
                             rep("Sparse ICA",length(M_PMSE3$sparse)),rep("Fast ICA",length(M_PMSE3$fast)),rep("Infomax ICA",length(M_PMSE3$infomax)),rep("SICA-EBM",length(M_PMSE1$ebm)),rep("Sparse Fast ICA",length(M_PMSE1$sfi))
                    ),
                    SNR=c(rep("Low",5*length(M_PMSE1$sparse)),
                          rep("Medium",5*length(M_PMSE1$sparse)),
                          rep("High",5*length(M_PMSE1$sparse)))
)


p2 = ggplot(M_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) +
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Mixing Matrix")


# For computation time
TIME = data.frame(TIME=c(TIME1$sparse,TIME1$fast,TIME1$infomax,TIME1$ebm,TIME1$sfi,
                         TIME2$sparse,TIME2$fast,TIME2$infomax,TIME2$ebm,TIME2$sfi,
                         TIME3$sparse,TIME3$fast,TIME3$infomax,TIME3$ebm,TIME3$sfi),
                  Method=c(rep("Sparse ICA",length(TIME1$sparse)),rep("Fast ICA",length(TIME1$fast)),rep("Infomax ICA",length(TIME1$infomax)),rep("SICA-EBM",length(TIME1$ebm)),rep("Sparse Fast ICA",length(TIME1$sfi)),
                           rep("Sparse ICA",length(TIME2$sparse)),rep("Fast ICA",length(TIME2$fast)),rep("Infomax ICA",length(TIME2$infomax)),rep("SICA-EBM",length(TIME2$ebm)),rep("Sparse Fast ICA",length(TIME2$sfi)),
                           rep("Sparse ICA",length(TIME3$sparse)),rep("Fast ICA",length(TIME3$fast)),rep("Infomax ICA",length(TIME3$infomax)),rep("SICA-EBM",length(TIME3$ebm)),rep("Sparse Fast ICA",length(TIME3$sfi))),
                  SNR=c(rep("Low",5*length(TIME1$sparse)),
                        rep("Medium",5*length(TIME2$sparse)),
                        rep("High",5*length(TIME3$sparse)))
)


p3 = ggplot(TIME, aes(factor(SNR,levels = c("Low","Medium","High")), TIME)) +
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="TIME",title = "Computation Time")#+ylim(c(0,0.25))


ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../../Figures/03_sim123_PRMSE_single.pdf",width = 10, height = 6)


apply(TIME1, 2, mean)
apply(TIME2, 2, mean)
apply(TIME3, 2, mean)


