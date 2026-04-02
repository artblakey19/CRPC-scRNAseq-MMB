remove(list = ls())
library(openxlsx)
library(scales)
library(ggplot2)
library(wesanderson)
library(Rgb)
library(egg)
library(grid)
library(gridExtra)
library(gtable)
library(ggsymbol)
setwd("../data_simulation/")
fdrpenal_size_height = unit(3.8, "in")
penal_size_width = unit(3, "in")
penal_size_height = unit(2.95, "in")

method_tot = c('Clipper',
               'BH-pool', 'BH-pair-correct', 'BH-pair-2as1','BH-pair-mis',
               'qvalue-pool','qvalue-pair-correct','qvalue-pair-2as1','qvalue-pair-mis',
               'locfdr-emp','locfdr-swap')
col_clipper = rgb(0, 114, 178, max = 255) # blue

col_pool = rgb(125, 22, 37, max = 255)
col_pair_correct = wes_palette("GrandBudapest1", 4, type = 'discrete')[2] 
col_paired_2as1 = rgb(206, 106, 73, max = 255) # orange ## 2as 1
col_paired_mis = rgb(254, 160, 135, max = 255) # pink ## mis

col_locfdr_e = rgb(230, 159, 0, max = 255) # ginger 
col_locfdr_p = wes_palette("Moonrise1", 4, type = "discrete")[1] # light yellow

col_palette = c(col_clipper, rep(c(col_pool,col_pair_correct,col_paired_2as1,col_paired_mis),2), col_locfdr_e, col_locfdr_p)
names(col_palette) = method_tot

pt_shape = c(17, rep(1, 4), rep(3, 4), 4,4)
pt_size = 4
pt_size2 = 3
pt_stroke = 1.5
names(pt_shape) = method_tot

mydarkgrey = "#252525"
FDR  = 0.05
q_tot = (1:10)/100
textsize = 20
titletextsize = 26
iter = 200
label_height = 0.005

format_method_name = function(methods){
  methods = as.character(methods)
  methods[methods == "clippermax_1ko"] = 'Clipper'
  methods[methods == "clipper"] = 'Clipper'  
  methods[methods == 'pooled'] = 'BH-pool'
  methods[methods == 'pair_c'] = 'BH-pair-correct'
  methods[methods == "pair_2as1"] = 'BH-pair-2as1'
  methods[methods == "pair_m"] = 'BH-pair-mis'
  methods[methods == "pair_mis"] = 'BH-pair-mis'
  methods[methods == 'pooled_q'] = 'qvalue-pool'
  methods[methods == "pair_c_q"] = 'qvalue-pair-correct'
  methods[methods == "pair_2as1_q"] = 'qvalue-pair-2as1'
  methods[methods == "pair_m_q"] = 'qvalue-pair-mis'
  methods[methods == "pair_mis_q"] = 'qvalue-pair-mis'
  
  methods[methods == 'locfdremp'] = 'locfdr-emp'
  methods[methods == 'locfdrper'] = 'locfdr-swap'
  
  methods = factor(methods, levels = method_tot)
  return(methods)
  
}





mydarkgrey = "#252525"
FDR  = 0.05
q_tot = (1:10)/100
textsize = 16
titletextsize = 26
iter = 200
label_height = 0.005
lwd = 1.25
blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()


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

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

transformy = function(x){
  tr =  (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
  tr[x==0] = 0
  return(tr)
}
transformy = function(x, FDR = 0.05){
  tr =  (FDR* log10(x/ FDR) + FDR)*(x > FDR) + x* ( x <= FDR)
  tr[x==0] = 0
  return(tr)
}
xlabels = c(1, 5, 25,  100)/100
xticks = transformy(xlabels)

ylabels = c(0, 0.25, 0.5, 0.75, 1)
method_subset = method_tot

ylabels1 = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels, FDR = 0.1)
yticks1 = transformy(ylabels1, FDR = 0.1)

############### simulations ##################
dist_tot = c(rep('norm', 2), rep('pois', 3), rep('nb', 1))
scenario_tot = c('homo' ,rep('hete', 5))

out_prop = c(0, 0, 0, 0, 0.1, 0)
analysis_tot = c('singlerep', 'singlerep', '1sided_unequal', 'multirep', 'multirep', '2sided')
folder_tot = c('1vs1_enrich', '1vs1_enrich', '2vs1_enrich', '3vs3_enrich', '3vs3_enrich', '3vs3_diff')
N = c(1000, 10000, 0, 0, 0, 0)
grid_ls = data.frame(folder_tot, analysis_tot, dist_tot, scenario_tot, out_prop, N)

p_ls1 = vector(mode = 'list',length = 6)

p_ls2 = vector(mode = 'list',length = 6)
p_ls = vector(mode = 'list',length = 6)

for(i in 1:nrow(grid_ls)){
  grid_i = grid_ls[i,]
  folder = grid_i[1]
  analysis = grid_i[2]
  out = grid_i[5]
  scenario = grid_i[4]
  dist = grid_i[3]
  n = grid_i[6]
  fdppow = readRDS(file = paste0(folder,'/sim_', analysis, '_',dist,'_',scenario,'.rds'))
  if (n!=0){
    dat = fdppow[[as.character(n)]]
  }else{
    dat = fdppow[[as.character(out)]]
  }  
  ####### clipper variants #####
  ylim_fdr = 0.10
  method_subset = c('Clipper',
                    'BH-pool', 
                    'BH-pair-correct',
                    'BH-pair-2as1','BH-pair-mis',
                    'locfdr-emp','locfdr-swap')
  # dat = dat[ !dat$methods == 'clipper',]
  dat$methods = format_method_name(dat$methods)
  
  
  p = ggplot(dat[dat$methods %in% method_subset ,], 
             aes(x = q, y = transformy(fdr, 0.1), group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.title=element_blank(),
          legend.position = 'none',
          panel.border = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.title.y = element_blank()) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, transformy(1,0.1)),
                       breaks = yticks1,
                       expand = c(0,0),
                       labels = ylabels1*100) + 
    geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y =  transformy(fdr,0.1), col = methods), lwd = lwd) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y =  transformy(fdr,0.1), shape = methods, col = methods), size = pt_size2) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  # p
  
  #### if first plot, extract legend and add ylabel and xlabel
  if(i == 6){
    p = p + theme(legend.position = 'right',
                  legend.text = element_text(size = 22),
                  legend.key.size = unit(0.4,'in'))
    legend = get_legend(p)
    p = p + theme(legend.position = 'none')
  }
  if(i%% 3 == 1 ){
    p = p + theme( axis.title.y = element_text( size = 22)) + 
      ylab('Actual FDR (%)')
  }
  ## add xlabel to the last row
  
  p = p + theme( axis.title.x = element_text(size = 22) ) + 
    xlab('Target FDR (%)')
  
  pfdr = set_panel_size(p, width  = penal_size_width, height = fdrpenal_size_height)
  p_ls1[[i]] = pfdr
  
  
  method_subset = method_tot
  dat$methods = format_method_name(dat$methods)
  dat = dat[dat$q %in% FDR,]
  #### fdr
  p = ggplot(dat[dat$methods %in% method_subset ,], 
             aes(x = transformy(fdr), y = pow, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.title=element_blank(),
          legend.position = 'none',
          panel.border = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.title.y = element_blank()) +
    geom_vline(xintercept = transformy(FDR), lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_point(aes(x = transformy(fdr), y = pow,shape = methods, col = methods, stroke = pt_stroke), size = pt_size) +
    scale_x_continuous(limits = c(0, transformy(1)),
                       breaks = xticks,
                       expand = c(0,0),
                       labels = xlabels*100) + 
    scale_y_continuous(limits = c(0, 1.05), expand = c(0,0),labels = ylabels*100) + 
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset]) + 
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset])
  
  
  # p
  
  #### if first plot, extract legend and add ylabel and xlabel

  if(i%% 3 == 1 ){
    p = p + theme( axis.title.y = element_text(size = 22)) + 
      ylab('Power (%)')
  }
  ## add xlabel to the last row
  # if(dist == 'nb'){
  p = p + theme( axis.title.x = element_text(size = 22)) + 
    xlab('Actual FDR (%)')
  # }
  
  ppow = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  p_ls2[[i]] = ppow
  p_ls[[i]] = list(p_ls1[[i]], p_ls2[[i]])
}


############### DE analysis on rna-seq ##################

p_ls = lapply(1:length(p_ls), function(i){
  print(i)
  p = p_ls[[i]]
  do.call('rbind', p)
})
pdf(file =  'figure2.pdf', height = 22, width = 15.3)
p = arrangeGrob(
  ncol = 6,nrow = 15, widths = c( 0.5, 3, 0.5, 3 , 0.5, 3), heights = c(0.4, 0.25, 0.2, 0.2, 0.25, 3.8,3.2, 0.5, 0.4, 0.25, 0.2, 0.2, 0.25, 3.8,3.2),
  textGrob(label = 'a', just = 'left', gp = gpar(fontsize= 36, fontface =2)),blankPlot,
  textGrob(label = 'b', just = 'left', gp = gpar(fontsize= 36, fontface =2)),blankPlot,
  textGrob(label = 'c', just = 'left', gp = gpar(fontsize= 36, fontface =2)),blankPlot,
  blankPlot,textGrob(label = '1vs1 Enrichment',just = 'center', gp = gpar(fontsize= 22, fontface =2)),
  blankPlot,textGrob(label = '1vs1 Enrichment',just = 'center', gp = gpar(fontsize= 22, fontface =2)),
  blankPlot,textGrob(label = '2vs1 Enrichment',just = 'center', gp = gpar(fontsize= 22, fontface =2)),
  blankPlot,textGrob(label = 'Gaussian',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Gaussian',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Poisson',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'd = 1000',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'd = 10000',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'd = 10000',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Homogeneous',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Heterogeneous',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Heterogeneous',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,do.call('cbind', p_ls1[1]), blankPlot,do.call('cbind', p_ls1[2]), blankPlot,do.call('cbind', p_ls1[3]),
  blankPlot,do.call('cbind', p_ls2[1]), blankPlot,do.call('cbind', p_ls2[2]), blankPlot,do.call('cbind', p_ls2[3]),
  blankPlot, blankPlot, blankPlot,blankPlot,blankPlot,blankPlot,
  textGrob(label = 'd', just = 'left', gp = gpar(fontsize= 36, fontface =2)),blankPlot,
  textGrob(label = 'e', just = 'left', gp = gpar(fontsize= 36, fontface =2)),blankPlot,
  textGrob(label = 'f', just = 'left', gp = gpar(fontsize= 36, fontface =2)),blankPlot,
  blankPlot,textGrob(label = '3vs3 Enrichment',just = 'center', gp = gpar(fontsize= 22, fontface =2)),
  blankPlot,textGrob(label = '3vs3 Enrichment',just = 'center', gp = gpar(fontsize= 22, fontface =2)),
  blankPlot,textGrob(label = '3vs3 Differential',just = 'center', gp = gpar(fontsize= 22, fontface =2)),
  blankPlot,textGrob(label = 'Gaussian',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Gaussian',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Negative binomial',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'No outlier',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Outlier',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'No outlier',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Heterogeneous',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Heterogeneous',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,textGrob(label = 'Heterogeneous',hjust = 0.4, gp = gpar(fontsize= 18)),
  blankPlot,do.call('cbind', p_ls1[4]), blankPlot,do.call('cbind', p_ls1[5]), blankPlot,do.call('cbind', p_ls1[6]),
  blankPlot,do.call('cbind', p_ls2[4]), blankPlot,do.call('cbind', p_ls2[5]), blankPlot,do.call('cbind', p_ls2[6])
)
grid.arrange( p, legend, ncol = 2,nrow = 1, widths = c(12, 3.3))
dev.off()
