###################################################################
# R codes for group-levels simulations
# Zihang Wang
# 1/5/2023
###################################################################

rm(list = ls())

library(ggplot2)
library(ggpubr)

load("../Data/01_group_sparse.RData")

###################################################################
# Figures for 40 restarts S and M, single restart TIME
###################################################################

# For S estimates
S_PMSE = data.frame(PMSE=c(S_PRMSE1$sparse40,S_PRMSE1$fast40,S_PRMSE1$infomax40,
                           S_PRMSE2$sparse40,S_PRMSE2$fast40,S_PRMSE2$infomax40,
                           S_PRMSE3$sparse40,S_PRMSE3$fast40,S_PRMSE3$infomax40),
                    Method=c(rep("Group Sparse ICA",length(S_PRMSE1$sparse40)),rep("Group Fast ICA",length(S_PRMSE1$fast40)),rep("Group Infomax ICA",length(S_PRMSE1$infomax40)),
                             rep("Group Sparse ICA",length(S_PRMSE2$sparse40)),rep("Group Fast ICA",length(S_PRMSE2$fast40)),rep("Group Infomax ICA",length(S_PRMSE2$infomax40)),
                             rep("Group Sparse ICA",length(S_PRMSE3$sparse40)),rep("Group Fast ICA",length(S_PRMSE3$fast40)),rep("Group Infomax ICA",length(S_PRMSE3$infomax40))
                    ),
                    SNR=c(rep("Low",3*length(S_PRMSE1$sparse40)),
                          rep("Medium",3*length(S_PRMSE2$sparse40)),
                          rep("High",3*length(S_PRMSE3$sparse40)))
)


p1 = ggplot(S_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) + 
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Source Signal")#+scale_fill_brewer(palette="BuPu")
#theme(legend.position  = "none")


# For M estimates
M_PMSE = data.frame(PMSE=c(M_avg_PRMSE1$sparse40,M_avg_PRMSE1$fast40,M_avg_PRMSE1$infomax40,
                           M_avg_PRMSE2$sparse40,M_avg_PRMSE2$fast40,M_avg_PRMSE2$infomax40,
                           M_avg_PRMSE3$sparse40,M_avg_PRMSE3$fast40,M_avg_PRMSE3$infomax40),
                    Method=c(rep("Group Sparse ICA",length(M_avg_PRMSE1$sparse40)),rep("Group Fast ICA",length(M_avg_PRMSE1$fast40)),rep("Group Infomax ICA",length(M_avg_PRMSE1$infomax40)),
                             rep("Group Sparse ICA",length(M_avg_PRMSE2$sparse40)),rep("Group Fast ICA",length(M_avg_PRMSE2$fast40)),rep("Group Infomax ICA",length(M_avg_PRMSE2$infomax40)),
                             rep("Group Sparse ICA",length(M_avg_PRMSE3$sparse40)),rep("Group Fast ICA",length(M_avg_PRMSE3$fast40)),rep("Group Infomax ICA",length(M_avg_PRMSE3$infomax40))
                    ),
                    SNR=c(rep("Low",3*length(M_avg_PRMSE1$sparse40)),
                          rep("Medium",3*length(M_avg_PRMSE2$sparse40)),
                          rep("High",3*length(M_avg_PRMSE3$sparse40)))
)


p2 = ggplot(M_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) + 
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Mixing Matrix")#+scale_fill_brewer(palette="BuPu")


# For computation time, we display the single version for fail comparisons
TIME = data.frame(TIME=c(TIME1$sparse,TIME1$fast,TIME1$infomax,
                         TIME2$sparse,TIME2$fast,TIME2$infomax,
                         TIME3$sparse,TIME3$fast,TIME3$infomax),
                  Method=c(rep("Sparse ICA",length(TIME1$sparse)),rep("Fast ICA",length(TIME1$fast)),rep("Infomax ICA",length(TIME1$infomax)),
                           rep("Sparse ICA",length(TIME2$sparse)),rep("Fast ICA",length(TIME2$fast)),rep("Infomax ICA",length(TIME2$infomax)),
                           rep("Sparse ICA",length(TIME3$sparse)),rep("Fast ICA",length(TIME3$fast)),rep("Infomax ICA",length(TIME3$infomax))),
                  SNR=c(rep("Low",3*length(TIME1$sparse)),
                        rep("Medium",3*length(TIME2$sparse)),
                        rep("High",3*length(TIME3$sparse)))
)


p3 = ggplot(TIME, aes(factor(SNR,levels = c("Low","Medium","High")), TIME)) + 
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="Seconds",title = "Computation Time")#+scale_fill_brewer(palette="BuPu")


ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../Figures/02_group_sparse_40restarts.png",width = 10, height = 6)


###################################################################
# Figure for single restart S, M, and TIME
###################################################################

# For S estimates
S_PMSE = data.frame(PMSE=c(S_PRMSE1$sparse,S_PRMSE1$fast,S_PRMSE1$infomax,
                           S_PRMSE2$sparse,S_PRMSE2$fast,S_PRMSE2$infomax,
                           S_PRMSE3$sparse,S_PRMSE3$fast,S_PRMSE3$infomax),
                    Method=c(rep("Group Sparse ICA",length(S_PRMSE1$sparse)),rep("Group Fast ICA",length(S_PRMSE1$fast)),rep("Group Infomax ICA",length(S_PRMSE1$infomax)),
                             rep("Group Sparse ICA",length(S_PRMSE2$sparse)),rep("Group Fast ICA",length(S_PRMSE2$fast)),rep("Group Infomax ICA",length(S_PRMSE2$infomax)),
                             rep("Group Sparse ICA",length(S_PRMSE3$sparse)),rep("Group Fast ICA",length(S_PRMSE3$fast)),rep("Group Infomax ICA",length(S_PRMSE3$infomax))
                    ),
                    SNR=c(rep("Low",3*length(S_PRMSE1$sparse)),
                          rep("Medium",3*length(S_PRMSE2$sparse)),
                          rep("High",3*length(S_PRMSE3$sparse)))
)


p1 = ggplot(S_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) + 
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Source Signal")#+scale_fill_brewer(palette="BuPu")
#theme(legend.position  = "none")


# For M estimates
M_PMSE = data.frame(PMSE=c(M_avg_PRMSE1$sparse,M_avg_PRMSE1$fast,M_avg_PRMSE1$infomax,
                           M_avg_PRMSE2$sparse,M_avg_PRMSE2$fast,M_avg_PRMSE2$infomax,
                           M_avg_PRMSE3$sparse,M_avg_PRMSE3$fast,M_avg_PRMSE3$infomax),
                    Method=c(rep("Group Sparse ICA",length(M_avg_PRMSE1$sparse)),rep("Group Fast ICA",length(M_avg_PRMSE1$fast)),rep("Group Infomax ICA",length(M_avg_PRMSE1$infomax)),
                             rep("Group Sparse ICA",length(M_avg_PRMSE2$sparse)),rep("Group Fast ICA",length(M_avg_PRMSE2$fast)),rep("Group Infomax ICA",length(M_avg_PRMSE2$infomax)),
                             rep("Group Sparse ICA",length(M_avg_PRMSE3$sparse)),rep("Group Fast ICA",length(M_avg_PRMSE3$fast)),rep("Group Infomax ICA",length(M_avg_PRMSE3$infomax))
                    ),
                    SNR=c(rep("Low",3*length(M_avg_PRMSE1$sparse)),
                          rep("Medium",3*length(M_avg_PRMSE2$sparse)),
                          rep("High",3*length(M_avg_PRMSE3$sparse)))
)


p2 = ggplot(M_PMSE, aes(factor(SNR,levels = c("Low","Medium","High")), PMSE)) + 
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="PRMSE",title = "Mixing Matrix")#+scale_fill_brewer(palette="BuPu")


# For computation time, we display the single version for fail comparisons
TIME = data.frame(TIME=c(TIME1$sparse,TIME1$fast,TIME1$infomax,
                         TIME2$sparse,TIME2$fast,TIME2$infomax,
                         TIME3$sparse,TIME3$fast,TIME3$infomax),
                  Method=c(rep("Sparse ICA",length(TIME1$sparse)),rep("Fast ICA",length(TIME1$fast)),rep("Infomax ICA",length(TIME1$infomax)),
                           rep("Sparse ICA",length(TIME2$sparse)),rep("Fast ICA",length(TIME2$fast)),rep("Infomax ICA",length(TIME2$infomax)),
                           rep("Sparse ICA",length(TIME3$sparse)),rep("Fast ICA",length(TIME3$fast)),rep("Infomax ICA",length(TIME3$infomax))),
                  SNR=c(rep("Low",3*length(TIME1$sparse)),
                        rep("Medium",3*length(TIME2$sparse)),
                        rep("High",3*length(TIME3$sparse)))
)


p3 = ggplot(TIME, aes(factor(SNR,levels = c("Low","Medium","High")), TIME)) + 
  geom_boxplot(aes(fill = Method))+ theme_bw() +labs(x="SNR Levels",y="Seconds",title = "Computation Time")#+scale_fill_brewer(palette="BuPu")


ggarrange(p1, p2, p3,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../Figures/02_group_sparse_single.png",width = 10, height = 6)

