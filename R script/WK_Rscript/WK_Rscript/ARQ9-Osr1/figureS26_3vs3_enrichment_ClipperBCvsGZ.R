setwd('../data_simulation/3vs3_enrich/')
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
method_tot = c('Clipper-BC', 'Clipper-GZ')
col_BC = rgb(0, 114, 178, max = 255) # blue

col_GZ = rgb(0, 20, 138, max = 255) # dark blue
col_palette = c(col_BC, col_GZ)
names(col_palette) = method_tot
pt_shape = c(1, 2)
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
  methods[methods == 'pooled2'] = 'FDP-pool'
  methods[methods == 'fdp_pair'] = 'FDP-pair'
  methods[methods == "pair_2as1"] = 'BH-pair-2as1'
  methods[methods == "pair_mis"] = 'BH-pair-mis'
  methods[methods == "clipper"] = 'Clipper-BC'
  methods[methods == "clippermax_bc"] = "Clipper-max-BC"
  methods[methods == 'clipperdiff_bc'] = 'Clipper-diff-BC'
  methods[methods == "clippermax_bh"] = "Clipper-max-BH"
  methods[methods == "clipperdiff_bh"] = "Clipper-diff-BH"
  methods[methods == "clipper_1sided_ko"] = "Clipper-GZ"
  methods[methods == "clipper_2sided_ko"] = "Clipper-2sided-KO"
  methods[methods == 'locfdremp'] = 'locfdr-emp'
  methods[methods == 'pair_c'] = 'BH-pair-correct'
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
  fdppow = readRDS(file = paste0('sim_multirep_',dist,'_',scenario,'.rds'))
  dat = fdppow[[which(out_prop %in% prop)]]
  
  ####### clipper variants #####
  ylim_fdr = 0.10
  
  method_subset = c('Clipper-BC', "Clipper-GZ")
  
  dat$methods = format_method_name(dat$methods)
  pdf(file = paste0('lineplot_fdr_multirep_GZ_',dist,'_',scenario, '_',prop,'.pdf'), height = 4, width = 6)
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
                       breaks = c(0.05, 0.10, 0.15, 0.2),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) + 
    coord_cartesian( ylim=c(0, ylim_fdr)) +
    geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = fdr, col = methods)) +
    scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
    geom_point(aes(q, y = fdr, shape = methods, col = methods)) +
    scale_shape_manual(values = pt_shape[names(pt_shape) %in% method_subset])
  # p
  grid.newpage()
  grid.draw(set_panel_size(p, width  = unit(2, "in"),
                           height = unit(1.8, "in")))
  dev.off()
  
  
  pdf(file = paste0('lineplot_pow_multirep_GZ_',dist,'_',scenario, '_',prop,'.pdf'), height = 8, width = 6)
  p = ggplot(dat[dat$methods %in% method_subset ,], aes(x = q, y = pow, group = methods)) +
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
                           height = unit(3, "in")))
  dev.off()
  
}


# ############### FDR = 0.05 barplots ##################
# ylim_fdr = 0.1
# barwidth = 0.65
# dist_tot = c('norm','pois','nb')
# scenario_tot = c('homo','hete')
# out_prop = c(0, 0.1, 0.2)
# method_subset = c('Clipper','BH-pool', 'FDP-pool','FDP-pair', 'BH-pair-correct','BH-pair-2as1',
#                   'locfdr-emp', 'locfdr-perm')
# 
# grid_ls = expand.grid(list(dist_tot, scenario_tot), stringsAsFactors = F)
# for(i in 1:nrow(grid_ls)){
#   grid_i = grid_ls[i,]
#   dist = grid_i[1]
#   scenario = grid_i[2]
#   
#   fdppow = readRDS(file = paste0('sim_multirep_',dist,'_',scenario,'.rds'))
#   dat = Reduce(rbind.data.frame, fdppow[which(out_prop %in% out_prop)])
#   dat$methods = format_method_name(dat$methods)
#   dat$prop = factor(rep(out_prop, each = nrow(fdppow[[1]])))
#   dat$fdr_sd  = dat$fdr_sd / sqrt(iter)
#   dat$pow_sd  = dat$pow_sd/sqrt(iter)
#   
  # pdf(file = paste0('barplot_fdr_multirep_',dist,'_', scenario,'.pdf'), height = 2, width = 8.5)
  # p = ggplot(dat[dat$methods %in% method_subset & dat$q == FDR,],
  #            aes(x = prop, y = fdr)) +
  #   theme(text = element_text(size=textsize, color = 'black'),
  #         panel.grid = element_blank(),
  #         panel.background = element_blank(),
  #         axis.text.x = element_text( colour = 1),
  #         axis.line.y = element_line(size = 0.5),
  #         axis.line.x = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.text.y = element_text( colour = 1),
  #         legend.title=element_blank(),
  #         legend.key = element_rect(colour = NA, fill = NA),
  #         axis.title.y = element_blank()) +
  #   coord_cartesian( ylim=c(0, ylim_fdr)) +
  #   geom_bar(aes(x = prop, fill = methods, y = pmin(ylim_fdr,fdr)),
  #            stat="identity",
  #            width = barwidth,
  #            position = position_dodge(barwidth) ) +
  #   geom_errorbar(width = 0.5,
  #                 position = position_dodge(barwidth),
  #                 aes(ymin = fdr, ymax = pmin(fdr + fdr_sd,ylim_fdr + 0.03), color = methods)) +
  #   geom_label(aes(x = prop,
  #                  y = label_height,
  #                  group = methods,
  #                  label = format_text(fdr*100) ),
  #              fill = alpha("#F0F0F0", 0.5),
  #              label.r = unit(0.1, "lines"),
  #              label.padding = unit(0.1, "lines"),
  #              label.size = 0,
  #              size = 3,
  #              color = 1,
  #              vjust = 0,
  #              position=position_dodge(barwidth)) +
  #   scale_y_continuous(limits = c(0, 1),
  #                      breaks = c(0, 0.05, 0.10, 0.15),
  #                      expand = c(0,0),
  #                      labels = function(x) paste0(x*100)) +
  #   scale_x_discrete(expand = c(0.22,0.22),label = c("No missing/outlier", "Outlier", "Missing value")) +
  #   geom_hline(yintercept=0.05,
  #              color ="#525252",
  #              lwd = 0.5,
  #              lty ='dashed') +
  #   scale_color_manual(values = col_palette[ names(col_palette) %in% method_subset]) +
  #   scale_fill_manual(values = col_palette[ names(col_palette) %in% method_subset])
  # grid.newpage()
  # grid.draw(set_panel_size(p, width  = unit(6, "in"),
  #                          height = unit(1.6, "in")))
  # dev.off()
  # 
  # ylim_pow = 1
  # pdf(file = paste0('barplot_pow_multirep_',dist,'_', scenario,'.pdf'), height = 2, width = 8.5)
  # p = ggplot(dat[dat$methods %in% method_subset & dat$q == FDR,],
  #            aes(x = prop, y = pow)) +
  #   theme(text = element_text(size=textsize, color = 'black'),
  #         panel.grid = element_blank(),
  #         panel.background = element_blank(),
  #         axis.text.x = element_text( colour = 1),
  #         axis.line.y = element_line(size = 0.5),
  #         axis.line.x = element_blank(),
  #         axis.title.x = element_blank(),
  #         axis.text.y = element_text( colour = 1),
  #         legend.title=element_blank(),
  #         legend.key = element_rect(colour = NA, fill = NA),
  #         axis.title.y = element_blank()) +
  #   coord_cartesian( ylim=c(0, ylim_pow)) +
  #   geom_bar(aes(x = prop, fill = methods, y = pmin(ylim_pow,pow)),
  #            stat="identity",
  #            width = barwidth,
  #            position = position_dodge(barwidth) )+
  #   geom_errorbar(width = 0.5,
  #                 position = position_dodge(barwidth),
  #                 aes(ymin = pow, ymax = pmin(pow + pow_sd,ylim_pow + 0.03), color = methods)) +
  #   scale_y_continuous(limits = c(0, 1),
  #                      # breaks = c(0, 0.05, 0.10, 0.15),
  #                      expand = c(0,0),
  #                      labels = function(x) paste0(x*100)) +
  #   scale_x_discrete(expand = c(0.22,0.22),label = c("No missing/outlier", "Outlier", "Missing value")) +
  #   geom_label(aes(x = prop, 
  #                  y = label_height*6, 
  #                  group = methods,
  #                  label = format_text(pow*100)),
  #              fill = alpha("#F0F0F0", 0.5),
  #              label.r = unit(0.1, "lines"),
  #              label.padding = unit(0.1, "lines"),
  #              label.size = 0,
  #              size = 3,
  #              color = 1,
  #              vjust = 0,
  #              position=position_dodge(barwidth)) + 
  #   scale_color_manual(values = col_palette[names(col_palette) %in% method_subset]) +
  #   scale_fill_manual(values = col_palette[names(col_palette) %in% method_subset])
  # grid.newpage()
  # grid.draw(set_panel_size(p, width  = unit(6, "in"),
  #                          height = unit(1.6, "in")))
  # dev.off()
  # 
  # 
  
  
# }

