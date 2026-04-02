remove(list = ls())
library(parallel)
library(DESeq2)
library(edgeR)
library(IHW)
library(preprocessCore)
source('clipper042520_w2sided.R')
source('auxiliaryfunctions_070321.R', echo=TRUE)
filename = "DEG_spikeinFollow2019_wbatcheffect_48vs48data_070321"
FDR = (1:10)/100
r1 = 3
r2 = 3
count = readRDS('count_matrices_gierlinski.rds')
count1 = count$count1
count2 = count$count2
prop0_1 = apply(count1, 1, function(x){
  mean(x==0)
})
prop0_2 = apply(count2, 1, function(x){
  mean(x==0)
})
hist(prop0_1)
r1_tot = ncol(count1)
r2_tot = ncol(count2)
geneid = rownames(count1)
# ############ use edgeR to normalize ############
# dat = cbind(count1, count2)
# cond_idx = rep(c(1,2), c(r1_tot, r2_tot))
# y <- DGEList(counts = dat, group = cond_idx)
# y <- calcNormFactors(y)
# count_norm_edger = cpm(y) #### normalized matrix
# count1 = count_norm_edger[, (1:r1_tot)]
# count2 = count_norm_edger[,-(1:r1_tot)]
# # ############ quantile normalize ############
# count1 = normalize.quantiles(as.matrix(count1),copy=TRUE)
# count2 = normalize.quantiles(as.matrix(count2),copy = T)
# rownames(count1) = geneid
# rownames(count2) = geneid
idx0_1 = apply(count1, 1, function(x){
  any(x == 0)
})
idx0_2 = apply(count2, 1, function(x){
  any(x == 0)
})

idx2fd = which( (!idx0_1) & (!idx0_2))

########
set.seed(12345)
m1 = rowMeans(count1)
m2 = rowMeans(count2)
fd =  (m1 + 1)/(m2+1)
fd = fd[(fd) >= 1.5 ] ### 335
# sd_ls = apply(count2, 1, sd)
# length(fd)
# saveRDS(fd, file = 'foldchange_fromGierlinski2015.rds')
# hist(fd, breaks = 1000)
# mu = mean(logfd)
# sd = sd(logfd)
trueDE_idx = sample(idx2fd, 0.3*nrow(count1))
# trueDE_idx = which(0 > 1)
trueDE = rownames(count1)[trueDE_idx]
fd_trueDE = sample_fd(fd = fd, n = length(trueDE_idx))
hist(fd)
hist(fd_trueDE)
fd_tot = rep(1, nrow(count1))
fd_tot[trueDE_idx] = fd_trueDE
saveRDS(list(fd_tot = fd_tot, trueDE_idx = trueDE_idx), file = paste0(filename, '_trueFD&DEidx.rds'))
count2_wfd = apply(count2, 2, function(x){x*fd_tot})
re = list(count_wofd = count2, count_wfd = count2_wfd, trueDE_idx = trueDE_idx)
saveRDS(re, file = 'realdata_yeast_2019.rds')


# *(as.numeric(runif(length(trueDE_idx)) < 0.5)*2 - 1)

# fd = rep(1, nrow(count1))
# fd[trueDE_idx] = fd_trueDE
hist(fd_trueDE, breaks = 100)
# generate_subcount= function(count, r1, r2, trueDE_idx, fd_trueDE){
#   
#   lgfd = log(fd_trueDE)
#   n_degene = length(trueDE_idx)
#   idx_cond1 = as.numeric(runif(n_degene) <= 1/2)
#   
#   scalar1 = exp(idx_cond1*lgfd)
#   scalar2 =  fd_trueDE/scalar1
#   
#   scaling_m1 = matrix(scalar1, ncol = 1) %*% matrix(rep(1, r1), nrow =1, byrow = T)
#   scaling_m2 = matrix(scalar2, ncol = 1) %*% matrix(rep(1, r2), nrow =1, byrow = T)
#   
#   r1_tot = ncol(count)
#   subcount = count[, sample(r1_tot, r1 + r2)]
#   subcount1 = subcount[, 1:r1]
#   subcount2 = subcount[, -(1:r1)]
#   
#   subcount1[trueDE_idx, ] = subcount1[trueDE_idx, ]*scaling_m1
#   subcount2[trueDE_idx, ] = subcount2[trueDE_idx, ]*scaling_m2
#   
#   
#   re = list(subcount1 = as.matrix(subcount1), subcount2 = as.matrix(subcount2))
#   
#   return(re)
#   
# }

ncores = 33
iter = 100
########
set.seed(1)

re = Compare(spikeInTarget = 'DEG')

saveRDS(re, file = paste0(filename,'_fd3divby2.rds'))



############ coarse plotting
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
re =  readRDS(file = paste0(filename,'_fd3divby2.rds'))

err_idx = (sapply(re, class) == 'try-error')

re = re[!err_idx]

re_test_uniformpval = test_uniformpval(re, trueDE_idx)
pval_test_uniformpval = lapply(re_test_uniformpval, function(x){x[2,]})
naprop_test_uniformpval = lapply(re_test_uniformpval, function(x){x[1,]})
saveRDS(list(pval = pval_test_uniformpval, 
             naprop = naprop_test_uniformpval), 'pval_NDEGifuniform_yeast_2019.rds')

# pdf(file = paste0('boxplot_',filename,'_nDEG_uniformtest.pdf'), height = 5, width = 7)
# re_test_uniformpval = lapply(pval_test_uniformpval, function(x){
#    x[!is.na(x)]
# })
# dev.off()
# saveRDS(re_test_uniformpval, file = 'pval_NDEGifuniform_yeast_2019.rds')

# vioplot( re_test_uniformpval,
#          # as.matrix(as.data.frame(re_test_uniformpval)),
#          main = filename)
# dev.off()


p_ls = coarseplot(re)
p1 = p_ls$p1
p2 = p_ls$p2
pdf(file = paste0(filename,'_fd3divby2.pdf'),  width = 10, height = 5)
p = cbind(p1, p2)
grid.arrange(p)
dev.off()
# 
# re = readRDS('DEG_spikeinFollow2019_wbatcheffect_48vs48data_070321.rds')
# dat = re
# print(length(re[[1]]$fdp))
# fdp = lapply(dat, '[[','fdp')
# fdp = lapply(1:length(re[[1]]$fdp), function(i){
#   rowMeans(sapply(fdp, '[[',i))
# })
# 
# pow = lapply(dat, '[[','pow')
# pow = lapply(1:length(re[[1]]$fdp), function(i){
#   rowMeans(sapply(pow, '[[',i))
# })
# 
# fdp
# pow
# ################### lineplots ###################
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(egg)
# library(gtable)
# 
# method_names = names(re[[1]]$fdp)
# dat = cbind.data.frame(fdr = unlist(fdp),
#                        pow = unlist(pow),
#                        methods = rep(method_names, each = 10),
#                        q = rep(FDR, length(method_names)))
# p1 = ggplot(dat,aes(x = q, y = fdr, group = methods)) +
#   geom_line(aes(q, y = fdr, col = methods), lwd = 1) + 
#   geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5)  
#   
# p2 = ggplot(dat,aes(x = q, y = pow, group = methods)) +
#   geom_line(aes(q, y = pow, col = methods), lwd = 1) 
# 
# pdf(file = 'fdr_DEG_spikeinFollow2019_wbatcheffect_48vs48data_070321.pdf',  width = 5, height = 5)
# p1
# dev.off()
# pdf(file = 'pow_DEG_spikeinFollow2019_wbatcheffect_48vs48data_070321.pdf',  width = 5, height = 5)
# p2
# dev.off()
# methods_sub = c('DESeq2','DESeq2 (IHW)','Clipper (DESeq2 normalized)','edgeR','edgeR (IHW)','Clipper (edgeR normalized)','Clipper')
# 
# 
# transformy = function(x){
#   (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
# }
# dat = cbind.data.frame(fdr = unlist(fdp),
#                        pow = unlist(pow),
#                        methods = rep(methods_sub, each = 10), 
#                        q = rep(FDR, length(methods_sub)))
# dat$fdr_tr = transformy(dat$fdr)
# 
# 
# textsize = 16
# pt_shape = c(1:9)
# 
# 
# ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
# yticks = transformy(ylabels)
# 
# col_clipper = rgb(0, 114, 178, max = 255)
# col_edger = "#D69840"  # yellow
# col_edger_ihw = "#FDB552"
# col_clipper_edger =  "#3C8A88"    # light blue
# col_deseq2 = "#B03A3D" # red
# col_deseq2_ihw = "#E48068"
# col_clipper_deseq2 = "#8A91C4"
# col_palette = c(
#   deseq2 = col_deseq2,
#   deseq2_ihw = col_deseq2_ihw,
#   clipper_deseq2 = col_clipper_deseq2,
#   edger = col_edger,
#   edger_ihw = col_edger_ihw,
#   clipper_edger = col_clipper_edger,
#   clipper = col_clipper
# )
# names(col_palette) = c('DESeq2','DESeq2 (IHW)','Clipper (DESeq2 normalized)','edgeR','edgeR (IHW)','Clipper (edgeR normalized)','Clipper')
# lwd = 0.75
# pt_size = 2
# 
# dat = dat[dat$methods %in% methods_sub,]
# 
# pdf(file = paste0('lineplot_DEanalysis_fdr_spikeinFollow2019_070321.pdf'), height = 5, width = 5)
# p = ggplot(dat, 
#            aes(x = q, y = fdr_tr, group = methods)) +
#   theme(text = element_text(size=textsize, color = 'black'),
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.text = element_text( colour = 1, size = textsize),
#         axis.line.y = element_line(size = 0.5),
#         axis.line.x = element_line(size = 0.5),
#         # axis.title.x = element_blank(),
#         legend.title=element_blank(),
#         legend.key = element_rect(colour = NA, fill = NA)) +
#   scale_x_continuous(limits = c(0, 0.11),
#                      breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
#                      expand = c(0,0),
#                      labels = function(x) paste0(x*100)) +
#   scale_y_continuous(limits = c(0, transformy(1.05)),
#                      breaks = yticks,
#                      expand = c(0,0),
#                      labels = ylabels*100) + 
#   # coord_cartesian( ylim=c(0, 1.1)) +
#   geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
#   geom_line(aes(q, y = fdr_tr, col = methods), lwd = lwd) +
#   scale_color_manual(values = col_palette) +
#   geom_point(aes(q, y = fdr_tr, shape = methods, col = methods), stroke = 1,size = pt_size) +
#   scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Actual FDR (%)')
# # saveRDS(p, file='lineplot_DEanalysis_fdr.rds')
# 
# grid.newpage()
# grid.draw(set_panel_size(p, width  = unit(2, "in"),
#                          height = unit(2/0.11*transformy(1), "in")))
# dev.off()
# 
# pdf(file = paste0('lineplot_DEanalysis_pow_spikeinFollow2019_070321.pdf'), height = 5, width = 5)
# p = ggplot(dat, 
#            aes(x = q, y = pow, group = methods)) +
#   theme(text = element_text(size=textsize, color = 'black'),
#         panel.grid = element_blank(),
#         panel.background = element_blank(),
#         axis.line.y = element_line(size = 0.5),
#         axis.line.x = element_line(size = 0.5),
#         # axis.title.x = element_blank(),
#         axis.text = element_text( colour = 1, size = textsize),
#         legend.title=element_blank(),
#         legend.key = element_rect(colour = NA, fill = NA),
#         axis.title.y = element_blank()) +
#   scale_x_continuous(limits = c(0, 0.11),
#                      breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
#                      expand = c(0,0),
#                      labels = function(x) paste0(x*100)) +
#   scale_y_continuous(limits = c(0, 0.75),
#                      breaks = c(0.25, 0.5, 0.75),
#                      expand = c(0,0),
#                      labels = c(0.25, 0.5, 0.75)*100) + 
#   coord_cartesian( ylim=c(0, 0.70)) +
#   # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
#   geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
#   scale_color_manual(values = col_palette) +
#   geom_point(aes(q, y = pow, shape = methods, col = methods),stroke = 1,size = pt_size) +
#   scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)")
# # saveRDS(p, file='lineplot_DEanalysis_pow.rds')
# grid.newpage()
# grid.draw(set_panel_size(p, width  = unit(2, "in"),
#                          height = unit(2/0.11*transformy(1), "in")))
# dev.off()
# 
# 
# 
