remove(list = ls())
library(parallel)
library(DESeq2)
library(edgeR)
library(IHW)
source('clipper042520_w2sided.R')
source('auxiliaryfunctions_070321.R', echo=TRUE)
filename = "DEG_spikeinNDEGbypermute_17vs17data_070321"
FDR = (1:10)/100
# dat_tot = readRDS(file = paste0('trueDE_FDR', 100*FDR,'.rds'))
r1 = 3
r2 = 3
r = 17
count1 = readRDS('count1.rds')
count2 = readRDS('count2.rds')
##########
############ use edgeR to normalize ############
dat = cbind(count1, count2)
cond_idx = rep(c(1,2), each = r)
y <- DGEList(counts = dat, group = cond_idx)
y <- calcNormFactors(y)
count_norm_edger = cpm(y) #### normalized matrix
count1 = count_norm_edger[, (1:r)]
count2 = count_norm_edger[,-(1:r)]



m1 = rowMeans(count1)
m2 = rowMeans(count2)
logfd = log(base = 2, (m1 + 1)/(m2+1))

# idx_m1_na = m1 == 0.0
# idx_m2_na = m2 == 0.0
# hist(logfd)
# sum(logfd > 2, na.rm = T)
trueDE_idx = which(abs(logfd) > 4)
trueDE = rownames(count1)[trueDE_idx] ## 191 trueDE
count1_trueDE = count1[trueDE_idx,]
count2_trueDE = count2[trueDE_idx,]
fd_tot = rep(1, nrow(count1))
fd_tot[trueDE_idx] = (2^logfd)[trueDE_idx]
saveRDS(list(fd_tot = fd_tot, trueDE_idx = trueDE_idx), file = paste0(filename, '_trueFD&DEidx.rds'))
# count2_wfd = apply(count2, 2, function(x){x*fd_tot})
re = list(count1 = count2, count2 = count2, trueDE_idx = trueDE_idx)
saveRDS(re, file = 'realdata_human_perm.rds')
# saveRDS(list(count1 = count1, count2 = count2, trueDE = trueDE ), file = 'new.rds')
# hist(log(base = 2, (m1[trueDE_idx] + 1)/(m2[trueDE_idx]+ 1)), breaks = 100)
# plot(log(base = 2, m1+1), log(base = 2, m2 + 1), cex =0.3)
# points(log(base = 2, m1[trueDE_idx]+1), log(base = 2, m2[trueDE_idx] + 1), col = 2, cex = 0.3)
# points(log(base = 2, m1[-trueDE_idx]+1), log(base = 2, m2[-trueDE_idx] + 1), col = 3, cex = 0.3)

ncores = 34
iter = 100
########
set.seed(1)





re = Compare(spikeInTarget = 'NDEG')

# rowMeans(sapply(re, '[[','fdp'))
# rowMeans(sapply(re, '[[','pow'))

saveRDS(re, file = 'DEG_spikeinNDEGbypermute_17vs17data_070321.rds')


################### lineplots ###################
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
re = readRDS( 'DEG_spikeinNDEGbypermute_17vs17data_070321.rds')


re_test_uniformpval = test_uniformpval(re, trueDE_idx)
pval_test_uniformpval = lapply(re_test_uniformpval, function(x){x[2,]})
naprop_test_uniformpval = lapply(re_test_uniformpval, function(x){x[1,]})
saveRDS(list(pval = pval_test_uniformpval, 
             naprop = naprop_test_uniformpval), 'pval_NDEGifuniform_human_perm.rds')
# pdf(file = paste0('boxplot_',filename,'_nDEG_uniformtest.pdf'), height = 5, width = 7)
# re_test_uniformpval = lapply(re_test_uniformpval, function(x){x[!is.na(x)]})
# vioplot( re_test_uniformpval,
#          # as.matrix(as.data.frame(re_test_uniformpval)),
#          main = filename)
# dev.off()
# saveRDS(re_test_uniformpval, file = 'pval_NDEGifuniform_human_perm.rds')


p_ls = coarseplot(re)

p1 = p_ls$p1
p2 = p_ls$p2
pdf(file = 'DEG_spikeinNDEGbypermute_17vs17data_070321.pdf',  width = 10, height = 5)
p = cbind(p1, p2)
grid.arrange(p)
dev.off()

# ######## drawing pvalues
# pval = lapply(dat, '[[','pval')
# pval_deseq = sapply(pval,'[[','deseq2')
# # nullpval_edger = sapply(nullpval, '[[','edger')
# pvaladj_deseq = apply(pval_deseq, 2, function(x){
#   p.adjust(x, method = 'BH')
# })
# FDR = 0.05
# false_discovery = apply(pvaladj_deseq, 2, function(x){
#   disc = rownames(pval_deseq)[which(x <= FDR)]
#   disc[! (disc %in% trueDE)]
# })
# genes2examine = names(which(table(unlist(false_discovery))>= 9))
# # par(mfrow = c(3,3))
# textsize = 12
# p_ls = vector('list', 25)
# for(i in 1:25){
#   # hist(pval_deseq[rownames(pval_deseq) %in% genes2examine[i], ],
#   #      xlim = c(0,1), freq = F, main = genes2examine[i], xlab = '')
#   value = pval_deseq[rownames(pval_deseq) %in% genes2examine[i], ]
#   value = value[!is.na(value)]
#   df = data.frame(value = value)
#   p = ggplot(df, aes(x=value)) +
#     geom_histogram(aes(y =..density..),bins = 20, colour = 1, fill = NA, size = 0.5) +
#     theme(text = element_text(size=textsize, color = "#525252"),
#           # panel.grid = element_rect(),
#           panel.background = element_rect(fill = 'white', colour = "#525252", size = 0.3),
#           axis.text = element_text( colour = "#525252", size = textsize),
#           axis.line.y = element_line(size = 0.3, colour = "#525252"),
#           axis.line.x = element_line(size = 0.3, colour = "#525252"),
#           axis.title.x = element_blank(),
#           legend.title = element_blank(),
#           plot.title = element_text(color = 1),
#           # title = element_text(genes2examine[i]),
#           legend.key = element_rect(colour = NA, fill = NA)) +
#     ggtitle(genes2examine[i])+
#     geom_abline(slope = 0, intercept = 1, lty = 'dashed', color ="#525252", size = 0.5) 
#     
#   if(i > 4){
#     p = p + theme(axis.title.y = element_blank())
#   }
#   
#   p = set_panel_size(p, width = unit(2.5, 'in'), height = unit(2, 'in'))
#   p_ls[[i]] = p
# 
# }
# 
# blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
#   cowplot::theme_nothing()
# 
# # p_tot = arrangeGrob(nrow = 7, ncol = 7,
# #             p_ls[[1]],blankPlot, p_ls[[2]],blankPlot,  p_ls[[3]], blankPlot, p_ls[[4]],
# #             blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, 
# #             p_ls[[5]],blankPlot,  p_ls[[6]], blankPlot,  p_ls[[7]],blankPlot,  p_ls[[8]], 
# #             blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, 
# #             p_ls[[9]],blankPlot,  p_ls[[10]], blankPlot, p_ls[[11]], blankPlot, p_ls[[12]], 
# #             blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, blankPlot, 
# #             p_ls[[13]],blankPlot,  p_ls[[14]],blankPlot,  p_ls[[15]],blankPlot,  p_ls[[16]], 
# #             # p_ls[[17]], p_ls[[18]],  p_ls[[19]], p_ls[[20]],
# #             # p_ls[[21]], p_ls[[22]],  p_ls[[23]], p_ls[[24]], p_ls[[25]],  
# #             widths = c(3, 0.2, 3, 0.2, 3, 0.2, 3), heights = c(3, 0.3, 3, 0.3, 3, 0.3, 3))
# 
# 
# g1 <- rbind(p_ls[[1]], p_ls[[2]], p_ls[[3]], p_ls[[4]])
# g2 <- rbind(p_ls[[5]], p_ls[[6]], p_ls[[7]], p_ls[[8]])
# g3 <- rbind(p_ls[[9]], p_ls[[10]], p_ls[[11]], p_ls[[12]])
# g4 <- rbind(p_ls[[13]], p_ls[[14]], p_ls[[15]], p_ls[[16]])
# g <- cbind(g1, g2, g3, g4)
# # g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
# 
# 
# 
# pdf(file =  'hist_pvalue_deseq2.pdf', height = 11, width = 12)
# grid.newpage()
# grid.draw(g)
# dev.off()

# naidx = apply(pval_deseq, 1, function(x){
#   all(is.na(x)) | names(pval_deseq[[1]])
# })

# nullpval_deseq = nullpval_deseq[!naidx, ]
# 
# idx_sp = sample(nrow(nullpval_deseq), 25)
# pdf(file = 'nullpval_deseq.pdf', width = 10, height = 10)
# par(mfrow = c(5,5))
# for(i in 1:25){
#   hist(nullpval_deseq[idx_sp[i],], xlim = c(0,1))
# }
# dev.off()

# pow
# transformy = function(x){
#   (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
# }
# dat = cbind.data.frame(fdr = unlist(fdp),
#                        pow = unlist(pow),
#                        methods = rep(c('DESeq2','DESeq2 (IHW)','edgeR','edgeR (IHW)','Clipper'), each = 10), 
#                        q = rep(FDR, 5))
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
# methods_sub = c('DESeq2','DESeq2 (IHW)','edgeR','edgeR (IHW)','Clipper')
# dat = dat[dat$methods %in% methods_sub,]
# 
# pdf(file = paste0('lineplot_DEanalysis_fdr_DEG_spikeinNDEGbypermute_17vs17data_070321.pdf'), height = 5, width = 5)
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
# pdf(file = paste0('lineplot_DEanalysis_pow_DEG_spikeinNDEGbypermute_17vs17data_070321.pdf'), height = 5, width = 5)
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

# ######################## using all 17 replicates########################
# r1 = r2 = 17
# re_r17 = mclapply(1:iter, function(it){
#   set.seed(it)
#   dat_i = generate_subcount(r1, r2)
#   # subcount1 = dat_i$subcount1
#   # subcount2 = dat_i$subcount2
#   # plot(rowMeans(subcount1), rowMeans(subcount2), cex = 0.3)
#   # points(rowMeans(subcount1)[trueDE_idx],rowMeans(subcount2)[trueDE_idx], cex = 0.3, col =2)
#   # abline(a = 0, b = 1)
#   ##### deseq2 
#   re_deseq2 = suppressMessages(mydeseq2(count1 = round(dat_i$subcount1), count2 = round(dat_i$subcount2), FDR = FDR, trueDE = trueDE))
#   fdppow_deseq2 = re_deseq2$fdppow_ls
#   count_norm_deseq2 = re_deseq2$count_norm
#   
#   ##### edgeR
#   re_edger = myedger(count1 = dat_i$subcount1, count2 = dat_i$subcount2 , FDR = FDR, trueDE = trueDE)
#   fdppow_edger = re_edger$fdppow_ls
#   count_norm_edger = re_edger$count_norm
#   
#   ##### clipper 
#   re_clipper = clipper2sided(score_exp =  dat_i$subcount1 , 
#                              score_back =  dat_i$subcount2 , 
#                              FDR = FDR,
#                              importanceScore_method = 'diff',
#                              contrastScore_method = 'max',
#                              # nknockoff = 4,
#                              FDR_control_method = 'GZ',
#                              ifpowerful = F)
#   discovery = rownames(dat_i$subcount1)[re_clipper$results[[1]]$discovery]
#   fdp_clipper = sum(!discovery %in% trueDE)/max(1, length(discovery))
#   pow_clipper = sum(discovery %in% trueDE)/length(trueDE)
#   
#   re_clipper_log = clipper2sided(score_exp = log(base = 2, dat_i$subcount1 + 1), 
#                                  score_back = log(base = 2, dat_i$subcount2 + 1), 
#                                  FDR = FDR,
#                                  importanceScore_method = 'diff',
#                                  contrastScore_method = 'max',
#                                  # nknockoff = 4,
#                                  FDR_control_method = 'GZ',
#                                  ifpowerful = F)
#   discovery = rownames(dat_i$subcount1)[re_clipper_log$results[[1]]$discovery]
#   fdp_clipper_log = sum(!discovery %in% trueDE)/max(1, length(discovery))
#   pow_clipper_log = sum(discovery %in% trueDE)/length(trueDE)
#   
#   ##### clipper + deseq2
#   re_clipper_deseq2 = clipper2sided(score_exp =  count_norm_deseq2[,1:r1] , 
#                                     score_back =  count_norm_deseq2[,-(1:r1)], 
#                                     FDR = FDR,
#                                     importanceScore_method = 'diff',
#                                     contrastScore_method = 'max',
#                                     # nknockoff = 4,
#                                     FDR_control_method = 'GZ',
#                                     ifpowerful = F)
#   discovery = rownames(count_norm_deseq2)[re_clipper_deseq2$results[[1]]$discovery]
#   fdp_clipper_deseq2 = sum(!discovery %in% trueDE)/max(1, length(discovery))
#   pow_clipper_deseq2 = sum(discovery %in% trueDE)/length(trueDE)
#   
#   re_clipper_log_deseq2 = clipper2sided(score_exp = log(base = 2, count_norm_deseq2[,1:r1] + 1), 
#                                         score_back = log(base = 2, count_norm_deseq2[,-(1:r1)] + 1), 
#                                         FDR = FDR,
#                                         importanceScore_method = 'diff',
#                                         contrastScore_method = 'max',
#                                         # nknockoff = 4,
#                                         FDR_control_method = 'GZ',
#                                         ifpowerful = F)
#   discovery = rownames(count_norm_deseq2)[re_clipper_log_deseq2$results[[1]]$discovery]
#   fdp_clipper_log_deseq2 = sum(!discovery %in% trueDE)/max(1, length(discovery))
#   pow_clipper_log_deseq2 = sum(discovery %in% trueDE)/length(trueDE)
#   
#   ##### clipper + edgeR
#   re_clipper_edger = clipper2sided(score_exp =  count_norm_edger[,1:r1] , 
#                                    score_back =  count_norm_edger[,-(1:r1)] , 
#                                    FDR = FDR,
#                                    importanceScore_method = 'diff',
#                                    contrastScore_method = 'max',
#                                    # nknockoff = 4,
#                                    FDR_control_method = 'GZ',
#                                    ifpowerful = F)
#   discovery = rownames(count_norm_edger)[re_clipper_edger$results[[1]]$discovery]
#   fdp_clipper_edger = sum(!discovery %in% trueDE)/max(1, length(discovery))
#   pow_clipper_edger = sum(discovery %in% trueDE)/length(trueDE)
#   
#   
#   re_clipper_log_edger = clipper2sided(score_exp = log(base = 2, count_norm_edger[,1:r1] + 1), 
#                                        score_back = log(base = 2, count_norm_edger[,-(1:r1)] + 1), 
#                                        FDR = FDR,
#                                        importanceScore_method = 'diff',
#                                        contrastScore_method = 'max',
#                                        # nknockoff = 4,
#                                        FDR_control_method = 'GZ',
#                                        ifpowerful = F)
#   discovery = rownames(count_norm_edger)[re_clipper_log_edger$results[[1]]$discovery]
#   fdp_clipper_log_edger = sum(!discovery %in% trueDE)/max(1, length(discovery))
#   pow_clipper_log_edger = sum(discovery %in% trueDE)/length(trueDE)
#   
#   
#   fdp = c( deseq2 = fdppow_deseq2[1],
#            clipper_deseq2 = fdp_clipper_deseq2, 
#            clipper_log_deseq2 = fdp_clipper_log_deseq2,
#            edger =  fdppow_edger[1], 
#            clipper_edger = fdp_clipper_edger,
#            clipper_log_edger = fdp_clipper_log_edger,
#            clipper = fdp_clipper,
#            clipper_log = fdp_clipper_log)
#   pow = c( deseq2 = fdppow_deseq2[2],
#            clipper_deseq2 = pow_clipper_deseq2,
#            clipper_log_deseq2 = pow_clipper_log_deseq2,
#            edger =  fdppow_edger[2],
#            clipper_edger = pow_clipper_edger,
#            clipper_log_edger = pow_clipper_log_edger,
#            clipper = pow_clipper,
#            clipper_log = pow_clipper_log)
#   return(list(fdp = fdp, pow = pow))
# },mc.cores = ncores)
# 
# rowMeans(sapply(re_r17, '[[','fdp'))
# rowMeans(sapply(re_r17, '[[','pow'))
# 
# saveRDS(re_r17, file = paste0('DEanalysis_r17_FDR', 100*FDR,'.rds'))
# 
# 
# re_r17 = readRDS(paste0('DEanalysis_r17_FDR', 100*FDR,'.rds'))
# rowMeans(sapply(re_r17, '[[','fdp'))
# rowMeans(sapply(re_r17, '[[','pow'))
# 
# dat = readRDS( paste0('DEanalysis_FDR', 100*FDR,'.rds'))
# fdp = rowMeans(sapply(dat, '[[','fdp'))
# pow = rowMeans(sapply(dat, '[[','pow'))
# fdp
# pow
