remove(list = ls())
library(ggplot2)
library(egg)
library(grid)
library(Rgb)

format_text = function(x){
  re = sapply(x, function(x_i){
    if(x_i < 10){
      return(as.character(round(x_i, 1)))
    }else{
      return(as.character(round(x_i)))
    }
  })
  return(re)
}

transformy = function(x){
  (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
}

ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels)
textsize = 16

################################################
hic = readRDS( 'Diffinteraction_spikein_nbatcheffect.rds')
FDR = (1:10)/100
iter = 100
fdp = lapply(hic, '[[','fdp')

fdp = lapply(1:4, function(i){
  temp = sapply(fdp, '[[', i)
  fdr = rowMeans(temp)
  fdr_se = apply(temp, 1, sd)/sqrt(iter)
  cbind(fdr, fdr_se)
})

fdp = Reduce('rbind', fdp)

pow = lapply(hic,'[[','pow')
pow = lapply(1:4, function(i){
  temp = sapply(pow, '[[', i)
  pow = rowMeans(temp)
  pow_se = apply(temp, 1, sd)/sqrt(iter)
  cbind(pow, pow_se)
})
pow = Reduce('rbind', pow)


method = c('Clipper','multiHiCcompare','diffHic','FIND')
col_clipper = rgb(0, 114, 178, max = 255) # blue
col_hiccompare = '#E8BA00'
col_diffhic = '#90A93B'
col_find = '#B94345'
col_palette = c(col_clipper, col_hiccompare, col_diffhic, col_find)
names(col_palette) = method
pt_shape = c(1,2, 3,4)
names(pt_shape) = method

dat = cbind.data.frame(methods = factor(rep(method, each = 10), levels = method), 
                       fdp , 
                       pow,
                       q = (rep(FDR, 4)))

dat$fdr_tr = transformy(dat$fdr)



lwd = 0.75
pt_size = 2

pdf(file = paste0('lineplot_DIRanalysis_nbatch_fdr.pdf'), height = 5, width = 5)
p = ggplot(dat, 
           aes(x = q, y = fdr_tr, group = methods)) +
  theme(text = element_text(size=textsize, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text( colour = 1, size = textsize),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        axis.title.y = element_blank(),
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA)) +
  scale_x_continuous(limits = c(0, 0.11),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                     expand = c(0,0),
                     labels = function(x) paste0(x*100)) +
  scale_y_continuous(limits = c(0, transformy(1)),
                     breaks = yticks,
                     expand = c(0,0),
                     labels = ylabels*100) + 
  # coord_cartesian( ylim=c(0, 1.1)) +
  geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = fdr_tr, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = fdr_tr, shape = methods, col = methods), stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)")

saveRDS(p, file='lineplot_DIRanalysis_nbatch_fdr.rds')
grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

pdf(file = paste0('lineplot_DIRanalysis_nbatch_pow.pdf'), height = 5, width = 5)
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
                     breaks = c(0.25, 0.5,0.75,  1),
                     expand = c(0,0),
                     labels = c(0.25, 0.5,0.75,  1)*100) + 
  coord_cartesian( ylim=c(0, 1)) +
  # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = pow, shape = methods, col = methods),stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)")
saveRDS(p, file='lineplot_DIRanalysis_nbatch_pow.rds')
grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

