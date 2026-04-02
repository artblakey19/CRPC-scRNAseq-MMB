################### lineplots ###################
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
dat = readRDS( 'DEG_spikein111220.rds')

fdp = lapply(dat, '[[','fdp')
fdp = lapply(1:5, function(i){
  rowMeans(sapply(fdp, '[[',i))
})

pow = lapply(dat, '[[','pow')
pow = lapply(1:5, function(i){
  rowMeans(sapply(pow, '[[',i))
})

pval = lapply(dat, '[[','pval')
pval_deseq = sapply(pval,'[[','deseq2')
# nullpval_edger = sapply(nullpval, '[[','edger')
pvaladj_deseq = apply(pval_deseq, 2, function(x){
  p.adjust(x, method = 'BH')
})
FDR = 0.05
false_discovery = apply(pvaladj_deseq, 2, function(x){
  disc = rownames(pval_deseq)[which(x <= FDR)]
  disc[! (disc %in% trueDE)]
})
genes2examine = names(which(table(unlist(false_discovery))>= 9))
# par(mfrow = c(3,3))
textsize = 12
p_ls = vector('list', 25)
for(i in 1:25){
 
  value = pval_deseq[rownames(pval_deseq) %in% genes2examine[i], ]
  value = value[!is.na(value)]
  df = data.frame(value = value)
  p = ggplot(df, aes(x=value)) +
    geom_histogram(aes(y =..density..),bins = 20, colour = 1, fill = NA, size = 0.5) +
    theme(text = element_text(size=textsize, color = "#525252"),
          # panel.grid = element_rect(),
          panel.background = element_rect(fill = 'white', colour = "#525252", size = 0.3),
          axis.text = element_text( colour = "#525252", size = textsize),
          axis.line.y = element_line(size = 0.3, colour = "#525252"),
          axis.line.x = element_line(size = 0.3, colour = "#525252"),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(color = 1),
          # title = element_text(genes2examine[i]),
          legend.key = element_rect(colour = NA, fill = NA)) +
    ggtitle(genes2examine[i])+
    geom_abline(slope = 0, intercept = 1, lty = 'dashed', color ="#525252", size = 0.5) 
  
  if(i > 4){
    p = p + theme(axis.title.y = element_blank())
  }
  
  p = set_panel_size(p, width = unit(2.5, 'in'), height = unit(2, 'in'))
  p_ls[[i]] = p
  
}

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

g1 <- rbind(p_ls[[1]], p_ls[[2]], p_ls[[3]], p_ls[[4]])
g2 <- rbind(p_ls[[5]], p_ls[[6]], p_ls[[7]], p_ls[[8]])
g3 <- rbind(p_ls[[9]], p_ls[[10]], p_ls[[11]], p_ls[[12]])
g4 <- rbind(p_ls[[13]], p_ls[[14]], p_ls[[15]], p_ls[[16]])
g <- cbind(g1, g2, g3, g4)
# g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)



pdf(file =  'hist_pvalue_deseq2.pdf', height = 11, width = 12)
grid.newpage()
grid.draw(g)
dev.off()



pow
transformy = function(x){
  (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
}
dat = cbind.data.frame(fdr = unlist(fdp),
                       pow = unlist(pow),
                       methods = rep(c('DESeq2','DESeq2 (IHW)','edgeR','edgeR (IHW)','Clipper'), each = 10), 
                       q = rep(FDR, 5))
dat$fdr_tr = transformy(dat$fdr)


textsize = 16
pt_shape = c(1:9)


ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels)

col_clipper = rgb(0, 114, 178, max = 255)
col_edger = "#D69840"  # yellow
col_edger_ihw = "#FDB552"
col_clipper_edger =  "#3C8A88"    # light blue
col_deseq2 = "#B03A3D" # red
col_deseq2_ihw = "#E48068"
col_clipper_deseq2 = "#8A91C4"
col_palette = c(
  deseq2 = col_deseq2,
  deseq2_ihw = col_deseq2_ihw,
  clipper_deseq2 = col_clipper_deseq2,
  edger = col_edger,
  edger_ihw = col_edger_ihw,
  clipper_edger = col_clipper_edger,
  clipper = col_clipper
)
names(col_palette) = c('DESeq2','DESeq2 (IHW)','Clipper (DESeq2 normalized)','edgeR','edgeR (IHW)','Clipper (edgeR normalized)','Clipper')
lwd = 0.75
pt_size = 2

methods_sub = c('DESeq2','DESeq2 (IHW)','edgeR','edgeR (IHW)','Clipper')
dat = dat[dat$methods %in% methods_sub,]

pdf(file = paste0('lineplot_DEanalysis_fdr.pdf'), height = 5, width = 5)
p = ggplot(dat, 
           aes(x = q, y = fdr_tr, group = methods)) +
  theme(text = element_text(size=textsize, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text( colour = 1, size = textsize),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        # axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA)) +
  scale_x_continuous(limits = c(0, 0.11),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                     expand = c(0,0),
                     labels = function(x) paste0(x*100)) +
  scale_y_continuous(limits = c(0, transformy(1.05)),
                     breaks = yticks,
                     expand = c(0,0),
                     labels = ylabels*100) + 
  # coord_cartesian( ylim=c(0, 1.1)) +
  geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = fdr_tr, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = fdr_tr, shape = methods, col = methods), stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Actual FDR (%)')
saveRDS(p, file='lineplot_DEanalysis_fdr.rds')

grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

pdf(file = paste0('lineplot_DEanalysis_pow.pdf'), height = 5, width = 5)
p = ggplot(dat, 
           aes(x = q, y = pow, group = methods)) +
  theme(text = element_text(size=textsize, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        # axis.title.x = element_blank(),
        axis.text = element_text( colour = 1, size = textsize),
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y = element_blank()) +
  scale_x_continuous(limits = c(0, 0.11),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                     expand = c(0,0),
                     labels = function(x) paste0(x*100)) +
  scale_y_continuous(limits = c(0, 0.75),
                     breaks = c(0.25, 0.5, 0.75),
                     expand = c(0,0),
                     labels = c(0.25, 0.5, 0.75)*100) + 
  coord_cartesian( ylim=c(0, 0.70)) +
  # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = pow, shape = methods, col = methods),stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)")
saveRDS(p, file='lineplot_DEanalysis_pow.rds')
grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

