library(ComplexHeatmap)
library("ggplot2")
library("reshape2")
library("viridis")
library(tidyverse)
library("ggrepel")
library(ggpmisc)

set.seed(42)




setwd("~/GLP2R_final/scripts/VEP comparison: Figure 3")


# loading dataframe 
df <-read.csv("Data/cAMP_LogEC50_GEMME.csv", row.names = "mutant")

cor.test(df$cAMP_LogEC50, df$GEMME_pred, method = "spearman")


ggplot(df, aes(x= cAMP_LogEC50, y = GEMME_pred)) + 
  geom_point() + 
  geom_text_repel(aes(label = rownames(df)), size = 3) +
  stat_poly_line() +
  #stat_poly_eq() +
  xlab("cAMP Log(EC50)")+
  ylab("GEMME prediction score") +
  annotate(geom = "text", x = -9.54, y = 0, label = "rho", parse=TRUE, size = 5) +
  annotate(geom = "text", x = -9.4, y = 0, label = "= -0.38", size = 5) +
  theme_bw(base_size = 18) 
    