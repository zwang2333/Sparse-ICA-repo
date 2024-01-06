###################################################################
# R codes for simulation setting 2 (High-Dimensional) in single-subject simulations
# 03: organize results and make figures
# Zihang Wang
# 1/5/2023
###################################################################
rm(list = ls())

library(ggplot2)
library(ggpubr)

load("../../Results/02_highdim_sparse.RData")

###################################################################
# Figures for single restart version
###################################################################

# For S estimates
S_cor = data.frame(cor=c(S_cor1$sparse_IC1,S_cor1$fast_IC1,S_cor1$infomax_IC1,S_cor1$EBM_IC1,S_cor1$sfi_IC1,
                         S_cor1$sparse_IC2,S_cor1$fast_IC2,S_cor1$infomax_IC2,S_cor1$EBM_IC2,S_cor1$sfi_IC1,
                         S_cor1$sparse_IC3,S_cor1$fast_IC3,S_cor1$infomax_IC3,S_cor1$EBM_IC3,S_cor1$sfi_IC3),
                    Method=c(rep("Sparse ICA",length(S_cor1$sparse_IC1)),rep("Fast ICA",length(S_cor1$fast_IC1)),rep("Infomax ICA",length(S_cor1$infomax_IC1)),rep("SICA-EBM",length(S_cor1$EBM_IC1)),rep("Sparse Fast ICA",length(S_cor1$sfi_IC1)),
                             rep("Sparse ICA",length(S_cor1$sparse_IC2)),rep("Fast ICA",length(S_cor1$fast_IC2)),rep("Infomax ICA",length(S_cor1$infomax_IC2)),rep("SICA-EBM",length(S_cor1$EBM_IC2)),rep("Sparse Fast ICA",length(S_cor1$sfi_IC2)),
                             rep("Sparse ICA",length(S_cor1$sparse_IC3)),rep("Fast ICA",length(S_cor1$fast_IC3)),rep("Infomax ICA",length(S_cor1$infomax_IC3)),rep("SICA-EBM",length(S_cor1$EBM_IC3)),rep("Sparse Fast ICA",length(S_cor1$sfi_IC3))
                    ),
                    IC=c(rep("IC1",5*length(S_cor1$sparse_IC1)),
                          rep("IC2",5*length(S_cor1$sparse_IC2)),
                          rep("IC3",5*length(S_cor1$sparse_IC3)))
)

p1 = ggplot(S_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method)+ 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Source Signal")#+ylim(c(0.9,1))


# For M estimates
M_cor = data.frame(cor=c(M_cor1$sparse_IC1,M_cor1$fast_IC1,M_cor1$infomax_IC1,M_cor1$EBM_IC1,M_cor1$sfi_IC1,
                         M_cor1$sparse_IC2,M_cor1$fast_IC2,M_cor1$infomax_IC2,M_cor1$sfi_IC2,M_cor1$sfi_IC2,
                         M_cor1$sparse_IC3,M_cor1$fast_IC3,M_cor1$infomax_IC3,M_cor1$EBM_IC3,M_cor1$sfi_IC3),
                   Method=c(rep("Sparse ICA",length(M_cor1$sparse_IC1)),rep("Fast ICA",length(M_cor1$fast_IC1)),rep("Infomax ICA",length(M_cor1$infomax_IC1)),rep("SICA-EBM",length(M_cor1$EBM_IC1)),rep("Sparse Fast ICA",length(M_cor1$sfi_IC1)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC2)),rep("Fast ICA",length(M_cor1$fast_IC2)),rep("Infomax ICA",length(M_cor1$infomax_IC2)),rep("SICA-EBM",length(M_cor1$EBM_IC2)),rep("Sparse Fast ICA",length(M_cor1$sfi_IC2)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC3)),rep("Fast ICA",length(M_cor1$fast_IC3)),rep("Infomax ICA",length(M_cor1$infomax_IC3)),rep("SICA-EBM",length(M_cor1$EBM_IC3)),rep("Sparse Fast ICA",length(M_cor1$sfi_IC3))
                   ),
                   IC=c(rep("IC1",5*length(M_cor1$sparse_IC1)),
                        rep("IC2",5*length(M_cor1$sparse_IC2)),
                        rep("IC3",5*length(M_cor1$sparse_IC3)))
)


p2 = ggplot(M_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method) + 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Mixing Matrix")#+ylim(c(0.985,1))


ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../../Figures/03_highdim_sparse_single.pdf",width = 14, height = 6)

###################################################################
# Figures for single restart version, without SICA-EBM and SparseFastICA
###################################################################

# For S estimates
S_cor = data.frame(cor=c(S_cor1$sparse_IC1,S_cor1$fast_IC1,S_cor1$infomax_IC1,
                         S_cor1$sparse_IC2,S_cor1$fast_IC2,S_cor1$infomax_IC2,
                         S_cor1$sparse_IC3,S_cor1$fast_IC3,S_cor1$infomax_IC3),
                   Method=c(rep("Sparse ICA",length(S_cor1$sparse_IC1)),rep("Fast ICA",length(S_cor1$fast_IC1)),rep("Infomax ICA",length(S_cor1$infomax_IC1)),
                            rep("Sparse ICA",length(S_cor1$sparse_IC2)),rep("Fast ICA",length(S_cor1$fast_IC2)),rep("Infomax ICA",length(S_cor1$infomax_IC2)),
                            rep("Sparse ICA",length(S_cor1$sparse_IC3)),rep("Fast ICA",length(S_cor1$fast_IC3)),rep("Infomax ICA",length(S_cor1$infomax_IC3))
                   ),
                   IC=c(rep("IC1",3*length(S_cor1$sparse_IC1)),
                        rep("IC2",3*length(S_cor1$sparse_IC2)),
                        rep("IC3",3*length(S_cor1$sparse_IC3)))
)

p1 = ggplot(S_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method)+ 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Source Signal")#+ylim(c(0.9,1))


# For M estimates
M_cor = data.frame(cor=c(M_cor1$sparse_IC1,M_cor1$fast_IC1,M_cor1$infomax_IC1,
                         M_cor1$sparse_IC2,M_cor1$fast_IC2,M_cor1$infomax_IC2,
                         M_cor1$sparse_IC3,M_cor1$fast_IC3,M_cor1$infomax_IC3),
                   Method=c(rep("Sparse ICA",length(M_cor1$sparse_IC1)),rep("Fast ICA",length(M_cor1$fast_IC1)),rep("Infomax ICA",length(M_cor1$infomax_IC1)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC2)),rep("Fast ICA",length(M_cor1$fast_IC2)),rep("Infomax ICA",length(M_cor1$infomax_IC2)),
                            rep("Sparse ICA",length(M_cor1$sparse_IC3)),rep("Fast ICA",length(M_cor1$fast_IC3)),rep("Infomax ICA",length(M_cor1$infomax_IC3))
                   ),
                   IC=c(rep("IC1",3*length(M_cor1$sparse_IC1)),
                        rep("IC2",3*length(M_cor1$sparse_IC2)),
                        rep("IC3",3*length(M_cor1$sparse_IC3)))
)


p2 = ggplot(M_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method) + 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Mixing Matrix")#+ylim(c(0.985,1))


ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../../Figures/03_highdim_sparse_single_nomatlab.pdf",width = 14, height = 6)


###################################################################
# Figures for comparing single restart version and 40 restarts version
###################################################################

# For S estimates
S_cor = data.frame(cor=c(S_cor1$sparse_IC1,S_cor1$fast_IC1,S_cor1$infomax_IC1,S_cor1$sparse40_IC1,S_cor1$fast40_IC1,S_cor1$infomax40_IC1,
                         S_cor1$sparse_IC2,S_cor1$fast_IC2,S_cor1$infomax_IC2,S_cor1$sparse40_IC2,S_cor1$fast40_IC2,S_cor1$infomax40_IC2,
                         S_cor1$sparse_IC3,S_cor1$fast_IC3,S_cor1$infomax_IC3,S_cor1$sparse40_IC3,S_cor1$fast40_IC3,S_cor1$infomax40_IC3),
                   Method=c(rep("Sparse ICA-1",length(S_cor1$sparse_IC1)),rep("Fast ICA-1",length(S_cor1$fast_IC1)),rep("Infomax ICA-1",length(S_cor1$infomax_IC1)),rep("Sparse ICA-40",length(S_cor1$sparse40_IC1)),rep("Fast ICA-40",length(S_cor1$fast40_IC1)),rep("Infomax ICA-40",length(S_cor1$infomax40_IC1)),
                            rep("Sparse ICA-1",length(S_cor1$sparse_IC2)),rep("Fast ICA-1",length(S_cor1$fast_IC2)),rep("Infomax ICA-1",length(S_cor1$infomax_IC2)),rep("Sparse ICA-40",length(S_cor1$sparse40_IC2)),rep("Fast ICA-40",length(S_cor1$fast40_IC2)),rep("Infomax ICA-40",length(S_cor1$infomax40_IC2)),
                            rep("Sparse ICA-1",length(S_cor1$sparse_IC3)),rep("Fast ICA-1",length(S_cor1$fast_IC3)),rep("Infomax ICA-1",length(S_cor1$infomax_IC3)),rep("Sparse ICA-40",length(S_cor1$sparse40_IC3)),rep("Fast ICA-40",length(S_cor1$fast40_IC3)),rep("Infomax ICA-40",length(S_cor1$infomax40_IC3))
                   ),
                   IC=c(rep("IC1",6*length(S_cor1$sparse_IC1)),
                        rep("IC2",6*length(S_cor1$sparse_IC2)),
                        rep("IC3",6*length(S_cor1$sparse_IC3)))
)

p1 = ggplot(S_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method)+ 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Source Signal")#+ylim(c(0.9,1))


# For M estimates
M_cor = data.frame(cor=c(M_cor1$sparse_IC1,M_cor1$fast_IC1,M_cor1$infomax_IC1,M_cor1$sparse40_IC1,M_cor1$fast40_IC1,M_cor1$infomax40_IC1,
                         M_cor1$sparse_IC2,M_cor1$fast_IC2,M_cor1$infomax_IC2,M_cor1$sparse40_IC2,M_cor1$fast40_IC2,M_cor1$infomax40_IC2,
                         M_cor1$sparse_IC3,M_cor1$fast_IC3,M_cor1$infomax_IC3,M_cor1$sparse40_IC3,M_cor1$fast40_IC3,M_cor1$infomax40_IC3),
                   Method=c(rep("Sparse ICA-1",length(M_cor1$sparse_IC1)),rep("Fast ICA-1",length(M_cor1$fast_IC1)),rep("Infomax ICA-1",length(M_cor1$infomax_IC1)),rep("Sparse ICA-40",length(M_cor1$sparse40_IC1)),rep("Fast ICA-40",length(M_cor1$fast40_IC1)),rep("Infomax ICA-40",length(M_cor1$infomax40_IC1)),
                            rep("Sparse ICA-1",length(M_cor1$sparse_IC2)),rep("Fast ICA-1",length(M_cor1$fast_IC2)),rep("Infomax ICA-1",length(M_cor1$infomax_IC2)),rep("Sparse ICA-40",length(M_cor1$sparse40_IC2)),rep("Fast ICA-40",length(M_cor1$fast40_IC2)),rep("Infomax ICA-40",length(M_cor1$infomax40_IC2)),
                            rep("Sparse ICA-1",length(M_cor1$sparse_IC3)),rep("Fast ICA-1",length(M_cor1$fast_IC3)),rep("Infomax ICA-1",length(M_cor1$infomax_IC3)),rep("Sparse ICA-40",length(M_cor1$sparse40_IC3)),rep("Fast ICA-40",length(M_cor1$fast40_IC3)),rep("Infomax ICA-40",length(M_cor1$infomax40_IC3))
                   ),
                   IC=c(rep("IC1",6*length(M_cor1$sparse_IC1)),
                        rep("IC2",6*length(M_cor1$sparse_IC2)),
                        rep("IC3",6*length(M_cor1$sparse_IC3)))
)


p2 = ggplot(M_cor) + geom_boxplot(aes(x=IC, y=cor)) + facet_grid(~Method) + 
  theme_bw()+labs(x="Independent Components",y="Correlation",title = "Mixing Matrix")#+ylim(c(0.985,1))


ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("../../Figures/03_highdim_sparse_single_40restarts_nomatlab.pdf",width = 14, height = 6)


###################################################################
# Results for constructing the computation time table
###################################################################

apply(TIME1,2,mean)

mean(S_cor1$sparse_IC1)
sd(S_cor1$sparse_IC1)
mean(S_cor1$fast_IC1)
sd(S_cor1$fast_IC1)
mean(S_cor1$infomax_IC1)
sd(S_cor1$infomax_IC1)
mean(S_cor1$EBM_IC1)
sd(S_cor1$EBM_IC1)
mean(S_cor1$sfi_IC1)
sd(S_cor1$sfi_IC1)

mean(S_cor1$sparse_IC2)
sd(S_cor1$sparse_IC2)
mean(S_cor1$fast_IC2)
sd(S_cor1$fast_IC2)
mean(S_cor1$infomax_IC2)
sd(S_cor1$infomax_IC2)
mean(S_cor1$EBM_IC2)
sd(S_cor1$EBM_IC2)
mean(S_cor1$sfi_IC2)
sd(S_cor1$sfi_IC2)

mean(S_cor1$sparse_IC3)
sd(S_cor1$sparse_IC3)
mean(S_cor1$fast_IC3)
sd(S_cor1$fast_IC3)
mean(S_cor1$infomax_IC3)
sd(S_cor1$infomax_IC3)
mean(S_cor1$EBM_IC3)
sd(S_cor1$EBM_IC3)
mean(S_cor1$sfi_IC3)
sd(S_cor1$sfi_IC3)


mean(M_cor1$sparse_IC1)
sd(M_cor1$sparse_IC1)
mean(M_cor1$fast_IC1)
sd(M_cor1$fast_IC1)
mean(M_cor1$infomax_IC1)
sd(M_cor1$infomax_IC1)
mean(M_cor1$EBM_IC1)
sd(M_cor1$EBM_IC1)
mean(M_cor1$sfi_IC1)
sd(M_cor1$sfi_IC1)

mean(M_cor1$sparse_IC2)
sd(M_cor1$sparse_IC2)
mean(M_cor1$fast_IC2)
sd(M_cor1$fast_IC2)
mean(M_cor1$infomax_IC2)
sd(M_cor1$infomax_IC2)
mean(M_cor1$EBM_IC2)
sd(M_cor1$EBM_IC2)
mean(M_cor1$sfi_IC2)
sd(M_cor1$sfi_IC2)

mean(M_cor1$sparse_IC3)
sd(M_cor1$sparse_IC3)
mean(M_cor1$fast_IC3)
sd(M_cor1$fast_IC3)
mean(M_cor1$infomax_IC3)
sd(M_cor1$infomax_IC3)
mean(M_cor1$EBM_IC3)
sd(M_cor1$EBM_IC3)
mean(M_cor1$sfi_IC3)
sd(M_cor1$sfi_IC3)

boxplot(TIME1$sparse,TIME1$fast,TIME1$infomax,TIME1$sparse40,TIME1$fast40,TIME1$infomax40,
        names = c("Sparse ICA","Fast ICA","Infomax ICA","Sparse ICA-40","Fast ICA-40","Infomax ICA-40"),ylab="TIME")
summary(TIME1$sparse)
summary(TIME1$fast)
summary(TIME1$infomax)

hist(TIME1$sparse)
hist(TIME1$fast)
hist(TIME1$infomax)

sd(TIME1$sparse)
sd(TIME1$fast)
sd(TIME1$infomax)

boxplot(TIME1$sparse,TIME1$fast,TIME1$infomax)



