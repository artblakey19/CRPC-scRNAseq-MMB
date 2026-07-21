  library(Seurat)
  raw_CRPC1 <- Read10X(data.dir = "Raw_data/CRPC1")
  raw_CRPC1 <- CreateSeuratObject(counts = raw_CRPC1, project = "copykat.hCRPC", min.cells = 0, min.features = 0)
  exp.hCRPC1 <- as.matrix(raw_CRPC1@assays$RNA@counts)

write.table(exp.hCRPC1, file="exp.hCRPC1.txt", sep="\t", quote = FALSE, row.names = TRUE)