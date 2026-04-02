setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/RNAseq/TumorvsPIN/Clipper-1")

####Method1####
wk_RAW1 <- cbind(genes_df$Sample_31888_s7112, genes_df$Sample_31889_s7181)
wk_RAW1 <- as.data.frame(wk_RAW1)
colnames(wk_RAW1) <- c("31888_s7112", "31889_s7181")
rownames(wk_RAW1) <- row.names(genes_df)

wk_RAW2 <- cbind(genes_df$Sample_31885_s9364, genes_df$Sample_31886_s6562, genes_df$Sample_31892_s9369)
wk_RAW2 <- as.data.frame(wk_RAW2)
colnames(wk_RAW2) <- c("31885_s9364", "31886_s6562", "31892_s9369")
rownames(wk_RAW2) <- row.names(genes_df)

r1 = ncol(wk_RAW1)
r2 = ncol(wk_RAW2)
cond_idx = rep(2, r1 + r2)
cond_idx[1:r1] = 1
dat = cbind(wk_RAW1, wk_RAW2)
cond_idx = factor(cond_idx)
y <- DGEList(counts=dat,group=cond_idx)
y <- calcNormFactors(y)
count_norm = cpm(y)

re_clipper = Clipper(score.exp = log(base = 2, count_norm[,1:r1] + 1), 
                     score.back = log(base = 2, count_norm[,-(1:r1)] + 1), 
                     FDR = 0.25,
                     analysis = "differential")

head(re_clipper$contrast.score.value)
write.csv(re_clipper$contrast.score.value, "TumorvPIN.clipper.0.25.contrast.score.value.csv")

head(re_clipper$q)
write.csv(re_clipper$q, "TumorvPIN.clipper.0.25.q.csv")

head(re_clipper$discoveries)
write.csv(re_clipper$discoveries, "TumorvPIN.clipper.0.25.discoveries.csv")

####Method2####
wk_RAW <- cbind(genes_df$Sample_31885_s9364, genes_df$Sample_31886_s6562, genes_df$Sample_31892_s9369, genes_df$Sample_31888_s7112, genes_df$Sample_31889_s7181)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("31885_s9364", "31886_s6562", "31892_s9369", "31888_s7112", "31889_s7181")
rownames(wk_RAW) <- row.names(genes_df)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("PIN", "PIN", "PIN", "Tumor", "Tumor"), genes=data.frame(Chr=as.character(genes_df$seqnames), Start=genes_df$start, End=genes_df$end, Strand=as.character(genes_df$strand), ID=genes_df$gene_id, Symbol=genes_df$symbol, Length=genes_df$total_exon_length))
colnames(wk_RAW_DGEList) <- c("31885_s9364", "31886_s6562", "31892_s9369", "31888_s7112", "31889_s7181")
names(wk_RAW_DGEList)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "upperquartile")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

re_clipper = Clipper(score.exp = log(base = 2, CPM_w[,4:5] + 0.5), 
                     score.back = log(base = 2, CPM_w[,1:3] + 0.5), 
                     FDR = 0.05,
                     analysis = "differential")

head(re_clipper$contrast.score.value)
write.csv(re_clipper$contrast.score.value, "TumorvPIN.clipper.contrast.score.value.csv")

head(re_clipper$q)
write.csv(re_clipper$q, "TumorvPIN.clipper.q.csv")

head(re_clipper$discoveries)
write.csv(re_clipper$discoveries, "TumorvPIN.clipper.discoveries.csv")
