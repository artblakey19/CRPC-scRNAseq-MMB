remove(list = ls())
library(multiHiCcompare)
library(diffHic)
library(lattice)
library(edgeR)
library(FIND)
library(csaw)
library(truncnorm)
library(gridExtra)
source('functions.R')
library(ggplot2)
library(Rgb)
library(scales)

mean_exp = readRDS('mean_exp.rds')
mean_back = readRDS('mean_back.rds')

dat_exp = reshape2::melt(mean_exp, c("x", "y"), value.name = "z")
dat_back = reshape2::melt(mean_back, c("x", "y"), value.name = "z")

library(ggplot2)
library(Rgb)
library(scales)
myblue = rgb(0, 114, 178, max = 255) 
myred = rgb(125, 22, 37, max = 255)
dat = as.data.frame(dat_back)
colnames(dat) = c('x','y','z')
p1 = ggplot(dat, aes(x, y, fill= log10(z))) + 
  theme(text = element_text(size= 15, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 1),
        axis.text = element_text( colour = 1, size = 15),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none',
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA)) +
  scale_x_continuous( expand = c(0,0),position = "top") + 
  scale_y_continuous( expand = c(0,0)) + 
  scale_fill_gradient2(limits = c(0,6), midpoint = 3) +
  geom_tile() 



dat = as.data.frame(dat_exp)
colnames(dat) = c('x','y','z')
p2 = ggplot(dat, aes(x, y, fill= log10(z))) + 
  theme(text = element_text(size= 15, color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 1),
        axis.text = element_text( colour = 1, size = 15),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA)) +
  scale_x_continuous( expand = c(0,0),position = "top") + 
  scale_y_continuous( expand = c(0,0)) + 
  scale_fill_gradient2(limits = c(0,6), midpoint = 3) +
  geom_tile() 


p1 = set_panel_size(p1, width  = unit(2, "in"),
                    height = unit(2, "in"))
p2 = set_panel_size(p2, width  = unit(2, "in"),
                    height = unit(2, "in"))
pdf(file = 'heatmap_hic.pdf', height = 3, width = 7)
grid.arrange(p1, p2,ncol = 2,nrow = 1, widths = c( 3, 3), heights = c(2))
dev.off()

