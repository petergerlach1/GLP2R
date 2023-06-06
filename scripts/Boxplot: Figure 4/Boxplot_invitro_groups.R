library(ComplexHeatmap)
library("ggplot2")
library("reshape2")
library("viridis")
library(tidyverse)
library(NbClust)
library(tidyverse)
library(readxl)
library(factoextra)
library(NbClust)
library(circlize)
library(ggbeeswarm)
library(ggrepel)


setwd("~/GLP2R_final/scripts/Boxplot: Figure 4")


# loading dataframe 
df <-read.csv("Data/RaSP_GEMME_cluster.csv")





ggplot(df, aes(x=RaSP_pred, y=GEMME_pred, color = Cluster, label = mutant)) + 
  geom_point() + 
  geom_vline(xintercept = 2, linetype="dashed") + 
  geom_hline(yintercept = -3, linetype="dashed") +
  geom_text_repel() +
  scale_fill_discrete(guide = "none") + 
  guides(fill=guide_legend(override.aes = list(color = NA)), color = FALSE, shape = FALSE) +
  ylab("GEMME \u0394E [a.u.]") +
  #xlab("RaSP \u0394\u0394G [kcal*mol-1]") + 
  xlab(expression(RaSP~Delta*Delta*G~"["*kcal %*% mol^{-1}*"]")) +
  theme_bw()





