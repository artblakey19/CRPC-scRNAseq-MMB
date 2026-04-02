yeast_perm = readRDS("pval_NDEGifuniform_yeast_perm.rds")
pval_yeast_perm = lapply(yeast_perm$pval,function(x){x[!is.na(x)]})
yeast_2019 = readRDS("pval_NDEGifuniform_yeast_2019.rds")
pval_yeast_2019 = lapply(yeast_2019$pval,function(x){x[!is.na(x)]})

human_perm = readRDS("pval_NDEGifuniform_human_perm.rds")
pval_human_perm = lapply(human_perm$pval,function(x){x[!is.na(x)]})
human_2019 = readRDS("pval_NDEGifuniform_human_2019.rds")
pval_human_2019 = lapply(human_2019$pval,function(x){x[!is.na(x)]})

naprop_yeast_perm = lapply(yeast_perm$pval, function(x){mean(is.na(x))})
naprop_yeast_2019 = lapply(yeast_2019$pval, function(x){mean(is.na(x))})
naprop_human_perm = lapply(human_perm$pval, function(x){mean(is.na(x))})
naprop_human_2019 = lapply(human_2019$pval, function(x){mean(is.na(x))})



palatte = c(alpha(2, 0.3), alpha(4, 0.3))
n_breaks = 30
pdf(file = 'pvalues_ifpvalofNDEGfollowUnif.pdf', height = 8, width = 8)
par(mfrow = c(2,2))
hist(pval_yeast_perm$deseq2_norm,
     freq = F, 
     col =palatte[1],
     main = 'DESeq2 on yeast semi-synthetic datasets',
     xlab = 'p-values',border = NA, breaks = n_breaks, xlim = c(0,1))
hist(pval_yeast_2019$deseq2_norm, freq = F, add = T, col = palatte[2],border = NA, breaks = n_breaks)
legend('topright',legend = c('Strategy 1','Strategy 2'), fill = palatte,  bty = 'n')
# text(x = 0.8,y = 5, labels = paste('Prop. of NA p-values: ', paste(paste0("Strategy 1: ", round(naprop_yeast_perm$deseq2_norm,2)),
#                                                                    paste0("Strategy 2: ", round(naprop_yeast_2019$deseq2_norm,2)), sep ='\n'), sep ='\n'))

hist(pval_yeast_perm$edger_norm,
     freq = F, 
     col =palatte[1],
     main = 'edgeR on yeast semi-synthetic datasets',
     xlab = 'p-values',border = NA, breaks = n_breaks)
hist(pval_yeast_2019$edger_norm, freq = F, add = T, col = palatte[2],border = NA, breaks = n_breaks)
legend('topright',legend = c('Strategy 1','Strategy 2'), fill = palatte,  bty = 'n')
# text(x = 0.8,y = 3, labels = paste('Prop. of NA p-values: ', paste(paste0("Strategy 1: ", round(naprop_yeast_perm$edger_norm,2)),
#                                                                   paste0("Strategy 2: ", round(naprop_yeast_2019$edger_norm,2)), sep ='\n'), sep ='\n'))

hist(pval_human_perm$deseq2_norm,
     freq = F, 
     col =palatte[1],
     main = 'DESeq2 on human semi-synthetic datasets',
     xlab = 'p-values',border = NA, breaks = n_breaks)
hist(pval_human_2019$deseq2_norm, freq = F, add = T, col = palatte[2],border = NA, breaks = n_breaks)
legend('topright',legend = c('Strategy 1','Strategy 2'), fill = palatte,  bty = 'n')
# text(x = 0.8,y = 5, labels = paste('Prop. of NA p-values: ', paste(paste0("Strategy 1: ", round(naprop_human_perm$deseq2_norm,2)),
#                                                                   paste0("Strategy 2: ", round(naprop_human_2019$deseq2_norm,2)), sep ='\n'), sep ='\n'))


hist(pval_human_perm$edger_norm,
     freq = F, 
     col =palatte[1],
     main = 'edgeR on human semi-synthetic datasets',
     xlab = 'p-values',border = NA, breaks = n_breaks)
hist(pval_human_2019$edger_norm, freq = F, add = T, col = palatte[2],border = NA, breaks = n_breaks)
legend('topright',legend = c('Strategy 1','Strategy 2'), fill = palatte,  bty = 'n')
# text(x = 0.8,y = 1.5, labels = paste('Prop. of NA p-values: ', paste(paste0("Strategy 1: ", round(naprop_human_perm$edger_norm,2)),
#                                                                   paste0("Strategy 2: ", round(naprop_human_2019$edger_norm,2)), sep ='\n'), sep ='\n'))


dev.off()
