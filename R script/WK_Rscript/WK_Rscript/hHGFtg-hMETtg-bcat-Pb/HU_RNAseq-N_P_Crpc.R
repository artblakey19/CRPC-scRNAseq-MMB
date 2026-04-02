##Env
library(tidyverse)
library(readr)
library(dplyr)

# Make count table from count files
#FPKM/TCGA
ls()
df <- list.files(path='//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNA-Seq_analysis_Hu/All_TCGA_Pros',
                 pattern='tsv$',recursive = TRUE, full.names = TRUE)

df.list <- sapply(df, read.delim, simplify=FALSE,comment.char="#")
df2 <- dplyr::bind_cols (df.list)
#Curate
df3 <- df2 %>% select(!starts_with(c("gene_id","gene_type","unstranded",
                                     "stranded_first","stranded_second",
                                     "fpkm_uq_unstranded")))
df3 <- as.data.frame(df3)
write.csv(df3, "samples_TCGA_pros_fpkm.csv")
#Metadata
sample_cols <- df3 %>% select(starts_with("gene_name"))
sample_cols <- colnames(sample_cols)
sample_cols <- as.data.frame(sample_cols)
df <- as.data.frame(df)
data_meta <- cbind(df, sample_cols)
write.csv(df, "samples_TCGA_meta.csv")
write.csv(data_meta, "sample_cols_meta.csv")

#Sequence of joining#metadata
head(df)
df <- as.data.frame(df)
write.csv(df, "samples_TCGA_meta.csv")
#Join The Metadata Files
samples_TCGA_meta <- read_csv("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNA-Seq_analysis_Hu/samples_TCGA_meta.csv")
gdc_sample_metadata<- read.delim("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNA-Seq_analysis_Hu/gdc_sample_sheet.2022-10-28.tsv")
colnames(samples_TCGA_meta)[2] <- "File_ID"
tcga_metadat_curated <- right_join(samples_TCGA_meta, gdc_sample_metadata, by = "File_ID",sort = F)
sample_cols_curated <- full_join(sample_cols_meta, gdc_sample_metadata, by = "File_ID",sort = F)
colnames(gdc_sample_metadata)
write.csv(tcga_metadat_curated, "tcga_naive_metadata.csv")
write.csv(sample_cols_curated, "tcga_naive_metadata_master.csv")
          
#Curate the data set for only naive TCGA
naive_tcga <- samples_TCGA_pros_naive_fpkm 
colnames(naive_tcga)[1] = "Genes"
naive_tcga <- naive_tcga %>% select(!starts_with(c("gene_name","tpm")))
write.csv(naive_tcga, "samples_TCGA_pros_naive_fpkm.csv")
