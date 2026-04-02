library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
setwd("../data_DEG_analysis/")
protocol_list <- c('drop-seq', '10x')
DE_pair = c('CD4+ T cell','Cytotoxic T cell')
q_tot = seq(0.01, 0.1, 0.01)
textsize = 16
pt_shape = c(1:9)


ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels)

col_clipper = rgb(0, 114, 178, max = 255)
col_edger = "#D69840"  # yellow
col_mast = "#FDB552"
col_wilcoxon =  wes_palette("GrandBudapest1", 4, type = 'discrete')[2]
col_ttest = rgb(125, 22, 37, max = 255) # dark brown
col_monocle3 = "#E48068"

#col_clipper_deseq2 = "#8A91C4"
col_palette = c(
  ttest = col_ttest,
  wilcoxon = col_wilcoxon,
  edger = col_edger,
  clipper = col_clipper,
  mast = col_mast,
  monocle3 = col_monocle3
)
names(col_palette) = c('t-test','Wilcoxon', 'edgeR','Clipper',"MAST", "Monocle3")
lwd = 0.75
pt_size = 2


for (iter1 in 1:2){
  fdppow_ls <- readRDS(paste0('scdesign2_', protocol_list[iter1],'_fdppow_0619.rds'))
  fdr = Reduce(cbind, lapply(fdppow_ls, '[[', 'fdp'))
  pow = Reduce(cbind, lapply(fdppow_ls, '[[', 'pow'))
  
  
  methods = fdppow_ls[[1]]$methods
  q = fdppow_ls[[1]]$q
  re = cbind.data.frame(fdr = rowMeans(fdr), pow = rowMeans(pow), q=q, methods = methods)
  
  transformy = function(x){
    (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
  }
  dat <- re
  dat$fdr_tr = transformy(dat$fdr)
  
  
  methods = c("t-test", "Wilcoxon", "edgeR", "Clipper" ,   "MAST" ,    "Monocle3")
  
  pdf(paste0("Scdesign2_", protocol_list[iter1], "_fdr.pdf"), width = 5, height = 5)
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
          legend.position = 'none',
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
    scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Actual FDR (%)') +
    geom_line(data=dat[dat$methods == "Clipper",], aes(q, y = fdr_tr, col = methods)) + 
    geom_point(data=dat[dat$methods == "Clipper",], aes(q, y = fdr_tr, shape = methods, col = methods), size = pt_size)

  grid.newpage()
  grid.draw(p1)
  dev.off()
  p1 = set_panel_size(p, width  = unit(2, "in"),
                      height = unit(2/0.11*transformy(1), "in"))
 
  
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
    scale_y_continuous(limits = c(0, 1),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) + 
    coord_cartesian( ylim=c(0, 1)) +
    # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    geom_point(aes(q, y = pow, shape = methods, col = methods),stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") +  ylab('Power (%)')
  p = p + theme( axis.title.y = element_text()) + 
    ylab('Power (%)') + geom_line(data=dat[dat$methods == "Clipper",], aes(q, y = pow, col = methods)) + 
    geom_point(data=dat[dat$methods == "Clipper",], aes(q, y = pow, shape = methods, col = methods), size = pt_size)

  p2 = set_panel_size(p, width  = unit(2, "in"),height = unit(2/0.11*transformy(1), "in"))
  pdf(paste0("Scdesign2_", protocol_list[iter1], "_power.pdf"), width = 5, height = 5)
  grid.newpage()
  grid.draw(p2)
  dev.off()
  p_ls = list(p1, p2)               
  pdf(paste0("Scdesign2_", protocol_list[iter1], ".pdf"), width = 10, height = 5)
  grid.arrange(do.call('cbind', p_ls),  ncol = 1,nrow = 1)
  dev.off()
}
