library(SingleCellExperiment)
library(Clipper)
library(edgeR)
library(MAST)
library(SingleCellExperiment)
library(monocle3)

source("my_functions.R")
protocol_list <- c('drop-seq', '10x')
DE_pair = c('CD4+ T cell','Cytotoxic T cell')
q_tot = seq(0.01, 0.1, 0.01)
ncores = 8
for(iter1 in 1:length(protocol_list)){
  re = mclapply(1:15, function(iter2){
    sce <- readRDS(paste0('../scdesign2_simu_dataset/sce/sce_', protocol_list[iter1], '_', iter2, '_pbmc2_hca_new.rds'))
    expression_data <- assay(sce, "counts_nozi")
    DE_truth <- sce@rowRanges@elementMetadata$de_CD4..T.cell_Cytotoxic.T.cell
    trueidx <- which(DE_truth==1)
    score_exp <- expression_data[, colData(sce)$cell_type %in% DE_pair[1]]
    score_back <- expression_data[, colData(sce)$cell_type %in% DE_pair[2]]
    score_exp_log <- log(expression_data[, colData(sce)$cell_type %in% DE_pair[1]]+1)
    score_back_log <- log(expression_data[, colData(sce)$cell_type %in% DE_pair[2]]+1)
    
    
    # Clipper
    re_clipper <- Clipper(score_exp_log, score_back_log, 'd', FDR = q_tot)
    fdppow_clipper = compute_fdppow(discovery_ls = re_clipper$discoveries, trueidx = trueidx)
    
    # re_clipper_diff <- Clipper(score_exp, score_back, 'd', contrast.score = 'diff',  n.permutation = 50, FDR = q_tot)
    # fdppow_clipper_diff = compute_fdppow(discovery_ls = re_clipper_diff$discoveries, trueidx = trueidx)
    
    # ttest
    re_ttest <- my.ttest(score_exp_log, score_back_log)
    discovery_ls = find_discovery_wpval(pval = re_ttest, q = q_tot)
    fdppow_t = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    # wilcoxon test
    re_wilcoxon <- my.wilcoxon(score_exp, score_back)
    discovery_ls = find_discovery_wpval(pval = re_wilcoxon, q = q_tot)
    fdppow_wilcoxon = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    # MAST
    discovery_ls <- mymast2(score_exp, score_back, FDR = q_tot, genename = rownames(score_exp))
    fdppow_mast = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    # monocle3
    re_monocle3 <- my_monocle3Pip(score_exp, score_back, ncores = 1)
    discovery_monocle3 = find_discovery_wpval(pval = re_monocle3, q = q_tot)
    fdppow_monocle3 = compute_fdppow(discovery_ls = discovery_monocle3, trueidx = trueidx)
    
    # edgeR
    re_edger = myedger(count1 = score_exp, count2 = score_back, trueDE = rownames(score_exp)[trueidx], FDR = q_tot)
    fdppow_edger = re_edger$fdppow_ls
    
    # re_edger_det = myedger(count1 = score_exp, count2 = score_back, trueDE = rownames(score_exp)[trueidx], det = T, FDR = q_tot)
    # fdppow_edger_det = re_edger_det$fdppow_ls
    
    # re_clipper_edger <- clipper_imp(imp_ls = -log(edger_ls), FDR = q_tot)
    # fdppow_clipper_edger = compute_fdppow(discovery_ls = re_clipper_edger$discoveries, trueidx = trueidx)
    fdppow_all <- t(cbind(fdppow_clipper, fdppow_t, fdppow_wilcoxon, fdppow_edger, fdppow_mast, fdppow_monocle3))
    fdppow_all <- as.data.frame(fdppow_all)
    fdppow_all$q <- rep(q_tot, 6)
    fdppow_all$methods <- rep(c("Clipper", "t-test", "Wilcoxon", "edgeR", "MAST", "Monocle3"), each = 10)
    return(fdppow_all)
  }, mc.cores = ncores)
  saveRDS(re, paste0('scdesign2_', protocol_list[iter1],'_fdppow_0619.rds'))
}