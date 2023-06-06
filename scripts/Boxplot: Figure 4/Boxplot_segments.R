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


setwd("~/GLP2R_final/scripts/Boxplot: Figure 4")


# loading dataframe 
df <-read_excel("Data/all_variants_segmentconverted.xlsx")


ggplot(df, aes(x=Segment, y=GEMME_pred)) + 
  geom_boxplot(outlier.shape = NA, fill = "cornflowerblue") + 
  coord_flip() + 
  ylab("GEMME \u0394E [a.u.]") +
  xlab("") +
  theme_classic(base_size = 18)
