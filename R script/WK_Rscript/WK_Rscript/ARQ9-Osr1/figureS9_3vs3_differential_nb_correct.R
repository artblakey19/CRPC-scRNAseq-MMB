setwd('../data_simulation/3vs3_diff/')
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
method_tot = c('Clipper-diff(1KO)-GZ','Clipper-diff(3KO)-GZ', 'Clipper-diff(allKO)-GZ',
               'Clipper-max(1KO)-GZ','Clipper-max(3KO)-GZ', 'Clipper-max(allKO)-GZ',  'Clipper-powerful',
               'BH-pool', 'FDR-Pool', 'BH-pair-correct', 'BH-pair-correct(unknown dispersion)', 
               'BH-pair-2as1','BH-pair-mis',
               'locfdr-emp','locfdr-perm')
col_clipper_diff = "#3C8A89" # dark green
col_clipper_diff_3ko = "#05A6A4" # green
col_clipper_diff_1ko = "#A7E0E4" # light green 
col_clipper_max = "#152776" # dark purple 
col_clipper_max_3ko = "#4956AC" # purple
col_clipper_max_1ko = "#7F86C2" # light purple
# col_clipper_diff_1ko = rgb(0, 114, 178, max = 255) # blue
# col_clipper_diff_3ko = rgb(0, 20, 138, max = 255) # dark blue
# col_clipper_diff = rgb(86, 180, 233, max = 255)# sky blue 
# col_clipper_max_1ko = rgb(0, 114, 178, max = 255) # blue
# col_clipper_max_3ko = rgb(0, 20, 138, max = 255) # dark blue
# col_clipper_max = rgb(86, 180, 233, max = 255)# sky blue 
col_clipper_powerful = rgb(0, 114, 178, max = 255)# blue 
col_pooled = rgb(125, 22, 37, max = 255) # dark brown
col_pooled2 = rgb(180, 40, 70, max = 255)
col_pair_correct = wes_palette("GrandBudapest1", 4, type = 'discrete')[2] 
col_pair_correct2 = wes_palette("GrandBudapest1", 4, type = 'discrete')[4] 
col_paired_2as1 = rgb(206, 106, 73, max = 255) # orange ## 2as 1
col_paired_mis = rgb(254, 160, 135, max = 255) # pink ## mis
col_locfdr_e = rgb(230, 159, 0, max = 255) # ginger 
col_locfdr_p = wes_palette("Moonrise1", 4, type = "discrete")[1] # light yellow
col_palette = c(col_clipper_diff_1ko, col_clipper_diff_3ko, col_clipper_diff,
                col_clipper_max_1ko, col_clipper_max_3ko, col_clipper_max, col_clipper_powerful,
                col_pooled, col_pooled2, col_pair_correct, col_pair_correct2, col_paired_2as1, col_paired_mis,
  col_locfdr_e, col_locfdr_p)
names(col_palette) = method_tot
pt_shape = c(rep(2,3),rep(3,3), 4, rep(5,2), rep(6,2), 7:10)
names(pt_shape) = method_tot
mydarkgrey = "#252525"
FDR  = 0.05
q_tot = (1:10)/100
textsize = 20
iter = 200
label_height = 0.005

format_method_name = function(methods){
  methods = as.character(methods)
  methods[methods == 'pooled'] = 'BH-pool'
  methods[methods == 'pooled2'] = 'Pooled-fdr'
  methods[methods == "pair_2as1"] = 'BH-pair-2as1'
  methods[methods == "pair_m"] = 'BH-pair-mis'
  methods[methods == 'pair_c'] = 'BH-pair-correct'
  methods[methods == 'pair_c2'] = 'BH-pair-correct(unknown dispersion)'
  methods[methods == "clippermax_1ko"] = 'Clipper-max(1KO)-GZ'
  methods[methods == "clipperdiff_1ko"] = 'Clipper-diff(1KO)-GZ'
  methods[methods == "clippermax_3ko"] = 'Clipper-max(3KO)-GZ'
  methods[methods == "clipperdiff_3ko"] = 'Clipper-diff(3KO)-GZ'
  methods[methods == "clippermax"] = 'Clipper-max(allKO)-GZ'
  methods[methods == "clipperdiff"] = 'Clipper-diff(allKO)-GZ'
  methods[methods == "clipper_pow"] = 'Clipper-powerful'
  methods[methods == 'locfdremp'] = 'locfdr-emp'
  methods[methods == 'locfdrper'] = 'locfdr-perm'
  methods = factor(methods, levels = method_tot)
  return(methods)
}


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
  tr =  (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
  tr[x==0] = 0
  return(tr)
}
ylabels = c(2,4,6,8,10, 25,50, 75, 100)/100
yticks = transformy(ylabels)

############### line plots ##################
dist_tot = c('norm','pois','nb')
scenario_tot = c('homo','hete')
out_prop = c(0, 0.1)

grid_ls = expand.grid(list(dist_tot, scenario_tot, out_prop), stringsAsFactors = F)
for(i in 1:nrow(grid_ls)){
  grid_i = grid_ls[i,]
  dist = grid_i[1]
  scenario = grid_i[2]
  prop = grid_i[3]
  fdppow = readRDS(file = paste0('sim_2sided_',dist,'_',scenario,'.rds'))
  dat = fdppow[[which(out_prop %in% prop)]]
  
  ####### clipper vs BH based and local fdr ############
  
  method_subset  = c( 'Clipper-powerful',  'BH-pool', 'Pooled-fdr', 'BH-pair-correct','BH-pair-2as1','BH-pair-mis',
                                 'locfdr-emp','locfdr-perm')
  dat$methods = format_method_name(dat$methods)
  pdf(file = paste0('lineplot_fdr_2sided_clipperVSothers_',dist,'_',scenario, '_',prop,'.pdf'), height = 8, width = 6)
  p = ggplot(dat[dat$methods %in% method_subset ,], 
             aes(x = q, y = transformy(fdr), group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text( colour = 1),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text( colour = 1),
          legend.title=element_blank(),
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
    geom_line(aes(q, y = transformy(fdr), col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y = transformy(fdr), shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(2/0.11*transformy(1), "in")))
  dev.off()
  
  
  pdf(file = paste0('lineplot_pow_2sided_clipperVSothers_',dist,'_',scenario, '_',prop,'.pdf'), height = 9, width = 6)
  p = ggplot(dat[dat$methods %in% method_subset ,], aes(x = q, y = pow, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text( colour = 
                                        1),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text( colour = 1),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.title.y = element_blank()) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, 1.05),
                       breaks = c( 0.25, 0.5, 0.75, 1),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) + 
    coord_cartesian( ylim=c(0, 1)) +
    # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = pow, col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y = pow, shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  # p  
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(2, "in")))
  dev.off()
  
  ############ Clipper Variants ############
  method_subset =  c('Clipper-max(1KO)-GZ','Clipper-max(3KO)-GZ', 'Clipper-max(allKO)-GZ', 
                     'Clipper-diff(1KO)-GZ','Clipper-diff(3KO)-GZ', 'Clipper-diff(allKO)-GZ')
  ylim_fdr = 0.10
  dat$methods = format_method_name(dat$methods)
  pdf(file = paste0('lineplot_fdr_2sided_clippervariants_',dist,'_',scenario, '_',prop,'.pdf'), height = 4, width = 6)
  p = ggplot(dat[dat$methods %in% method_subset ,], 
             aes(x = q, y = fdr, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text( colour = 1),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text( colour = 1),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.title.y = element_blank()) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c( 0.05,0.1),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) + 
    coord_cartesian( ylim=c(0, ylim_fdr)) +
    geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = fdr, col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y = fdr, shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(1.8, "in")))
  dev.off()
  
  
  
  pdf(file = paste0('lineplot_pow_2sided_clippervariants_',dist,'_',scenario, '_',prop,'.pdf'), height = 4, width = 6)
  p = ggplot(dat[dat$methods %in% method_subset ,], aes(x = q, y = pow, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text( colour = 
                                        1),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text( colour = 1),
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
    geom_line(aes(q, y = pow, col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y = pow, shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  # p  
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(8, "in")))
  dev.off()
  
}

############### known vs. unknown dispersion #############

dist_tot = c('nb')
scenario_tot = c('homo','hete')
out_prop = c(0, 0.1)

grid_ls = expand.grid(list(dist_tot, scenario_tot, out_prop), stringsAsFactors = F)
for(i in 1:nrow(grid_ls)){
  grid_i = grid_ls[i,]
  dist = grid_i[1]
  scenario = grid_i[2]
  prop = grid_i[3]
  fdppow = readRDS(file = paste0('sim_2sided_',dist,'_',scenario,'.rds'))
  dat = fdppow[[which(out_prop %in% prop)]]
  
  ####### clipper vs BH based and local fdr ############
  
  method_subset  = c('Clipper-powerful', 'BH-pair-correct', 'BH-pair-correct(unknown dispersion)')
  dat$methods = format_method_name(dat$methods)
  nb_correct <- dat[dat$methods %in% method_subset ,]
  pdf(file = paste0('lineplot_fdr_2sided_correct_nb_',scenario, '_',prop,'.pdf'), height = 8, width = 8)
  p = ggplot(nb_correct, 
             aes(x = q, y = transformy(fdr), group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text( colour = 1),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text( colour = 1),
          legend.title=element_blank(),
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
    geom_line(aes(q, y = transformy(fdr), col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset], labels = c('Clipper', 'BH-pair-correct(known dispersion)','BH-pair-correct(unknown dispersion)')) +
    geom_point(aes(q, y = transformy(fdr), shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset], labels = c('Clipper', 'BH-pair-correct(known dispersion)','BH-pair-correct(unknown dispersion)'))
  
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(2/0.11*transformy(1), "in")))
  dev.off()
  
  
  pdf(file = paste0('lineplot_pow_2sided_correct_nb_',scenario, '_',prop,'.pdf'), height = 9, width = 8)
  p = ggplot(nb_correct, aes(x = q, y = pow, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text( colour = 
                                        1),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          axis.title.x = element_blank(),
          axis.text.y = element_text( colour = 1),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.title.y = element_blank()) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, 1.05),
                       breaks = c( 0.25, 0.5, 0.75, 1),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) + 
    coord_cartesian( ylim=c(0, 1)) +
    # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = pow, col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset], labels = c('Clipper', 'BH-pair-correct(known dispersion)','BH-pair-correct(unknown dispersion)')) +
    geom_point(aes(q, y = pow, shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset], labels = c('Clipper', 'BH-pair-correct(known dispersion)','BH-pair-correct(unknown dispersion)'))
  # p  
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(2, "in")))
  dev.off()
  
  
}
