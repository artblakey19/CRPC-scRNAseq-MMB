#Add necessary tools to library

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)

DefaultAssay(ARCre2) <- "RNA"

ARCre2_Krt8 <- subset(x=ARCre2, subset = Krt8 > 0)
Idents(object = ARCre2_Krt8) <- "Krt8Pos"
ARCre2_Krt8[["Krt8Exp"]] <- Idents(object = ARCre2_Krt8)
Idents(object = ARCre2_Krt8) <- "Krt8Exp"
ARCre2$Krt8Exp <- Idents(object = ARCre2_Krt8)

Idents(object = ARCre2) <- "Krt8Exp"
ARCre2_Krt8 <- subset(ARCre2, idents = c("Krt8Pos"))
DimPlot(ARCre2_Krt8, reduction = "umap", pt.size = 0.3)
table(Idents(ARCre2_Krt8))