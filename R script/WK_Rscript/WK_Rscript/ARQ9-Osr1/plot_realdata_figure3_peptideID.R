library(ggplot2)
library(egg)
library(grid)
library(Rgb)
setwd('~/Dropbox/clipper/proteomics')
FDR_ls = (1:10)/100
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

#### power of sequest higher 
fdrpow_qval_ls = readRDS('realdata_proteomics_sequest_qvalue_fdppow_ls.rds')
fdrpow_clipper_ls = readRDS('realdata_proteomics_sequest_clipper_fdppow_ls.rds')
fdrpow_qval_ls
fdrpow_clipper_ls

dat = cbind.data.frame(
  methods = rep(c('SEQUEST','Clipper'), each = length(FDR_ls)),
  rbind(t(fdrpow_qval_ls),t(fdrpow_clipper_ls)),
  q = rep(FDR_ls,2)
  
)
colnames(dat) = c('methods','fdp','pow','q')

dat$fdp_tr = transformy(dat$fdp)




ylim_fdr = 1.2
barwidth = 0.65
label_height = 0.005


ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels)

textsize = 16
lwd = 0.75
pt_size = 2
col_clipper = rgb(0, 114, 178, max = 255) # blue
col_mascot =  "#E8BA00" #yellow
col_palette = c(col_clipper, col_mascot)
names(col_palette) = c('Clipper','Mascot')
pt_shape = c(Clipper = 1, Mascot = 2)


pdf(file = paste0('lineplot_peptideID_fdr.pdf'), height = 5, width = 5)
p = ggplot(dat, 
           aes(x = q, y = fdp_tr, group = methods)) +
  theme(text = element_text(size=textsize, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text( colour = 1, size = textsize),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        # axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA)
        # axis.title.y = element_blank()
        ) +
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
  geom_line(aes(q, y = fdp_tr, col = methods),lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = fdp_tr, shape = methods, col = methods), size = pt_size,stroke = 1) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Actual FDP (%)')
# p
saveRDS(p, file='lineplot_peptideID_fdr.rds')

grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()


pdf(file = paste0('lineplot_peptideID_pow.pdf'), height = 5, width = 5)
p = ggplot(dat, 
           aes(x = q, y = pow, group = methods)) +
  theme(text = element_text(size=textsize, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text( colour = 1, size = textsize),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        # axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y = element_blank()) +
  scale_x_continuous(limits = c(0, 0.11),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                     expand = c(0,0),
                     labels = function(x) paste0(x*100)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0.25, 0.5, 0.75, 1),
                     expand = c(0,0),
                     labels = c(0.25, 0.5, 0.75, 1)*100) + 
  coord_cartesian( ylim=c(0, 0.5)) +
  # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = pow, col = methods),lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = pow, shape = methods, col = methods), stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") 
saveRDS(p, file='lineplot_peptideID_pow.rds')
grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()
