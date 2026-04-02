setwd('../data_simulation/3vs3_enrich/')
remove(list = ls())
library(parallel)
library(locfdr)
# library(truncnorm)
library(scales)
# library(truncdist)
library(ggplot2)
library(wesanderson)
library(egg)
library(grid)
library(Rgb)
library(scales)
library(ggplot2)
library(gridExtra)
library(gtable)
penal_size_width = unit(2, "in")
fdrpenal_size_height = unit(3.4, "in")
powpenal_size_height = unit(2, "in")

method_tot = c('Clipper','Clipper-diff-BC', "Clipper-max-BC","Clipper-t", 
               'BH-pool', 'BH-pair-correct','BH-pair-Wilcoxon','BH-pair-parametric','BH-pair-permutation',
               'BH-pair-2as1','BH-pair-mis',
               'locfdr-emp','locfdr-swap')
col_clipper = rgb(0, 114, 178, max = 255) # blue
col_clipperdiff_rc =  "#05A6A4" # green
col_clippermax_rc = "#4956AC" # purple
col_clipper_t = "#A7E0E4" # light green 
col_clippermax_bh = "#7F86C2" # light purple
col_clipper_1sided_KO = rgb(0, 20, 138, max = 255) # dark blue
col_clipper_2sided_KO = rgb(20, 60, 200, max = 255) # dark blue
# col_clipper_2s_ko = rgb(0, 80, 183, max = 255)# dark blue2
col_pooled = rgb(125, 22, 37, max = 255) # dark brown
col_pair_correct = wes_palette("GrandBudapest1", 4, type = 'discrete')[2] 
col_paired_2as1 = rgb(206, 106, 73, max = 255) # orange ## 2as 1
col_paired_mis = rgb(254, 160, 135, max = 255) # pink ## mis
col_paired_wilcoxon = rgb(180, 40, 70, max = 255) # orange ## 2as 1
col_paired_permutation = rgb(200, 90, 100, max = 255)# pink ## mis

col_locfdr_e = rgb(230, 159, 0, max = 255) # ginger 
col_locfdr_p = wes_palette("Moonrise1", 4, type = "discrete")[1] # light yellow
# col_palette = c(col_clipper, col_clipperdiff_rc, col_clippermax_rc,col_clipperdiff_bh,col_clippermax_bh,
#                 col_clipper_1sided_KO, col_clipper_2sided_KO,
#                 col_pooled, col_pooled2, col_pair_fdp, col_pair_correct, col_paired_2as1, col_paired_mis,
#                 col_locfdr_e, col_locfdr_p)

col_palette = c(col_clipper, col_clipperdiff_rc, col_clippermax_rc,col_clipper_t,
                col_pooled,  col_pair_correct, col_paired_wilcoxon, col_pair_correct, col_paired_permutation,
                col_paired_2as1, col_paired_mis,
                col_locfdr_e, col_locfdr_p)

names(col_palette) = method_tot
pt_shape = c(1, 2,2,3:4, 5,6,5, 7:11)
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
  methods[methods == 'pooled'] = 'BH-pool'
  methods[methods == 'pooled2'] = 'FDP-pool'
  methods[methods == 'fdp_pair'] = 'FDP-pair'
  methods[methods == "pair_2as1"] = 'BH-pair-2as1'
  methods[methods == "pair_mis"] = 'BH-pair-mis'
  methods[methods == "pair_m"] = 'BH-pair-mis'
  methods[methods == "clipper"] = 'Clipper'
  methods[methods == "clippermax_bc"] = "Clipper-max-BC"
  methods[methods == 'clipperdiff_bc'] = 'Clipper-diff-BC'
  methods[methods == "clipper_t"] = "Clipper-t"
  methods[methods == "clippermax_bh"] = "Clipper-max-BH"
  methods[methods == "clipperdiff_bh"] = "Clipper-diff-BH"
  methods[methods == "clipper_1sided_ko"] = "Clipper-1sided-KO"
  methods[methods == "clipper_2sided_ko"] = "Clipper-2sided-KO"
  methods[methods == "clipper_pow"] = 'Clipper'
  methods[methods == 'pair_wilcoxon'] = 'BH-pair-Wilcoxon'
  methods[methods == "pair_permutation"] = 'BH-pair-permutation'
  methods[methods == "pair_parametric"] = 'BH-pair-parametric'
  methods[methods == 'locfdremp'] = 'locfdr-emp'
  methods[methods == 'pair_c'] = 'BH-pair-correct'
  methods[methods == 'locfdrper'] = 'locfdr-swap'
  methods = factor(methods, levels = method_tot)
  return(methods)
  
}
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
ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels)

############### line plots ##################
dist_tot = c('norm', 'pois', 'nb')
scenario_tot = c('homo', 'hete')
out_prop = c(0, 0.1)

grid_ls = expand.grid(list(out_prop, scenario_tot, dist_tot), stringsAsFactors = F)

p_ls = vector(mode = 'list',length = 12)

for(i in 1:nrow(grid_ls)){
  grid_i = grid_ls[i,]
  out = grid_i[1]
  scenario = grid_i[2]
  dist = grid_i[3]
  fdppow = readRDS(file = paste0('sim_multirep_',dist,'_',scenario,'_t.rds'))
  dat = fdppow[[which(out_prop %in% out)]]
  if (!(dist == 'norm' &scenario == 'hete')){
    fdppow2 = readRDS(file = paste0('sim_multirep_',dist,'_',scenario,'.rds'))
    dat2 = fdppow2[[which(out_prop %in% out)]]
    dat = rbind(dat, dat2)
  }
  ####### clipper variants #####
  ylim_fdr = 0.10
  method_subset = c('Clipper', 'Clipper-t',
                    'BH-pool', 'BH-pair-correct',
                    'BH-pair-2as1','BH-pair-mis',
                    'locfdr-emp','locfdr-swap')
  
  dat$methods = format_method_name(dat$methods)
  #### fdr
  p = ggplot(dat[dat$methods %in% method_subset ,], 
             aes(x = q, y = transformy(fdr), group = methods)) +
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
    scale_y_continuous(limits = c(0, transformy(1)),
                       breaks = yticks,
                       expand = c(0,0),
                       labels = ylabels*100) + 
    geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y =  transformy(fdr), col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y =  transformy(fdr), shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  # p
  
  #### if first plot, extract legend and add ylabel and xlabel
  if(i == 1){
    p = p + theme(legend.position = 'right',
                  legend.text = element_text(size = textsize),
                  legend.key.size = unit(0.4,'in'))
    legend = get_legend(p)
    p = p + theme(legend.position = 'none')
  }
  if(i%% 2 == 1 ){
    p = p + theme( axis.title.y = element_text()) + 
      ylab('Actual FDR (%)')
  }
  
  pfdr = set_panel_size(p,width  = penal_size_width, height = fdrpenal_size_height)
  
  #### power
  p = ggplot(dat[dat$methods %in% method_subset ,], aes(x = q, y = pow, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          legend.position = 'none') +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, 1),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) + 
    coord_cartesian( ylim=c(0, 1)) +
    # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = pow, col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y = pow, shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset]) 
  
  ## if the first plot
  if(i%% 2 == 1 ){
    p = p + theme( axis.title.y = element_text()) + 
      ylab('Power (%)')
  }
  ## add xlabel to the last row
  if(dist == 'nb'){
    p = p + theme( axis.title.x = element_text()) + 
      xlab('Target FDR (%)')
  }
  
  ppow = set_panel_size(p,width  = penal_size_width, height = powpenal_size_height)
  
  
  p_ls[[i]] = list(pfdr, ppow)
}

### combine fdr and power
p_ls = lapply(1:length(p_ls), function(i){
  print(i)
  p = p_ls[[i]]
  do.call('rbind', p)
})

pdf(file =  'lineplot_enrichment_multirep_t.pdf', height = 23.4, width = 15.3)

p = arrangeGrob(
  ncol = 4,nrow = 8, widths = c(0.5, 3, 0.5,3), heights = c(0.3,0.2, 0.3, 3, 0.3, 3, 0.3, 3),
  blankPlot,textGrob(label = 'Homogeneous background',just = 'center', gp = gpar(fontsize= titletextsize, fontface =2)),blankPlot,textGrob(label = 'Heterogeneous background',just = 'center', gp = gpar(fontsize= titletextsize, fontface =2)),
  blankPlot,textGrob(label = 'No outlier                 Outlier',hjust = 0.4, gp = gpar(fontsize= textsize)),blankPlot,textGrob(label = 'No outlier                 Outlier',hjust = 0.5, gp = gpar(fontsize= textsize)),
  textGrob(label = 'a Gaussian',just = 'left', gp = gpar(fontsize= titletextsize, fontface =2)), blankPlot, blankPlot,blankPlot,
  blankPlot, do.call('cbind', p_ls[1:2]), blankPlot, do.call('cbind', p_ls[3:4]), 
  textGrob(label = 'b Poisson',just = 'left', gp = gpar(fontsize= titletextsize, fontface =2)), blankPlot, blankPlot,blankPlot,
  blankPlot, do.call('cbind', p_ls[5:6]), blankPlot, do.call('cbind', p_ls[7:8]), 
  textGrob(label = 'c Negative binomial',just = 'left', gp = gpar(fontsize= titletextsize, fontface =2)), blankPlot, blankPlot, blankPlot,
  blankPlot, do.call('cbind', p_ls[9:10]), blankPlot, do.call('cbind', p_ls[11:12])
)
grid.arrange(p, legend, ncol = 2,nrow = 1, widths = c(12, 3.3))
dev.off()