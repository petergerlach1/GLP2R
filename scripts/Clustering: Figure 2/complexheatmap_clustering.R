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
#set.seed(21)

setwd("~/GLP2R_project/final_heatmap_clustering")

# loading dataframe
df <-read.csv("Data/variant_zscore.csv", row.names = "mutant")



df$cAMP_LogEC50_zscore <- -1*df$cAMP_LogEC50_zscore
df$arrestin_LogEC50_zscore <- -1*df$arrestin_LogEC50_zscore



matrix_heatmap <- data.matrix(df)
newOrder <- c("cAMP_Emax_zscore", "cAMP_LogEC50_zscore", "arrestin_Emax_zscore", "arrestin_LogEC50_zscore")  # Specify the new order of columns
matrix_heatmap <- matrix_heatmap[, newOrder]

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

#row_km = 5
HM <- Heatmap(matrix_heatmap, column_names_rot = 45, row_km = 6,
              column_names_gp = gpar(fontsize = 10),
              row_names_gp = gpar(fontsize = 6), name = "Z-score",
              show_column_dend = FALSE, row_km_repeats = 100,
              col = col_fun)

draw(HM)

fviz_nbclust(df, kmeans, method = "wss") +
  #geom_vline(xintercept = 4, linetype = 2) +
  labs(subtitle = "Elbow method")

for (i in 1:length(row_order(HM))){
  if (i == 1) {
    clu <- t(t(row.names(matrix_heatmap[row_order(HM)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  }
  else {
    clu <- t(t(row.names(matrix_heatmap[row_order(HM)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
out

write.csv(out, "Data/regenie_clustering_groups.csv", row.names=FALSE)
