setwd("../codes_omics_data/DE analysis/intermediate_data/")
remove(list = ls())
library(parallel)
library(DESeq2)
library(edgeR)
library(IHW)
library(idr)
source('clipper042520_w2sided.R')


################### lineplots ###################
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
library(scales)
transformy = function(x){
  (0.1* log10(x/ 0.1) + 0.1)*(x > 0.1) + x* ( x <= 0.1)
}
textsize = 16
titletextsize = 22

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
  clipper_unnorm = col_clipper_deseq2,
  edger = col_edger,
  edger_ihw = col_edger_ihw,
  clipper_edgernorm = col_clipper_edger,
  clipper = col_clipper
)
names(col_palette) = c('DESeq2','DESeq2 (IHW)','Clipper_unnorm','edgeR','edgeR (IHW)', 'Clipper_edgernorm','Clipper')
lwd = 0.75
pt_size = 2
penal_size_width  = unit(2, "in")
penal_size_height = unit(2/0.11*transformy(1), "in")
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
dat2plot1 <- function(dat, use.norm){
  
  fdp = lapply(dat, '[[','fdp')
  fdp = lapply(1:10, function(i){
    rowMeans(sapply(fdp, '[[',i))
  })
  
  pow = lapply(dat, '[[','pow')
  pow = lapply(1:10, function(i){
    rowMeans(sapply(pow, '[[',i))
  })
  
  pow.i = lapply(dat, '[[','pow.i')
  pow.i = lapply(1:10, function(i){
    pow.mat <- sapply(pow.i, '[[',i)
    pow.mat[is.na(pow.mat)] = 0
    rowMeans(pow.mat)
  })
  
  
  dat = cbind.data.frame(fdr = unlist(fdp),
                         pow = unlist(pow),
                         pow.i = unlist(pow.i),
                         methods = rep(c('DESeq2 unnorm','DESeq2 norm','DESeq2 (IHW) unnorm','DESeq2 (IHW) norm',
                                         'edgeR unnorm','edgeR norm', 'edgeR (IHW) unnorm','edgeR (IHW) norm',
                                         'Clipper_unnorm', 'Clipper_edgernorm'), each = 10), 
                         q = rep(FDR, 10))
  dat$fdr_tr = transformy(dat$fdr)
  
  
  if(use.norm == 1){
    methods_sub = c('DESeq2 norm','DESeq2 (IHW) norm','edgeR norm','edgeR (IHW) norm','Clipper_edgernorm')
  }else{
    methods_sub = c('DESeq2 unnorm','DESeq2 (IHW) unnorm','edgeR unnorm','edgeR (IHW) unnorm','Clipper_unnorm')
  }
  dat = dat[dat$methods %in% methods_sub,]
  dat$methods <- rep( c('DESeq2','DESeq2 (IHW)','edgeR','edgeR (IHW)','Clipper'), each = 10)
  p = ggplot(dat, 
             aes(x = q, y = fdr_tr, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          legend.position = 'none',
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
  #saveRDS(p, file='lineplot_DEanalysis_fdr.rds')
  p1 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  
  
  #pdf(file = paste0('lineplot_DEanalysis_pow.pdf'), height = 5, width = 5)
  p = ggplot(dat, 
             aes(x = q, y = pow, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.position = 'none',
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0.25, 0.5, 0.75),
                       expand = c(0,0),
                       labels = c(0.25, 0.5, 0.75)*100) + 
    coord_cartesian( ylim=c(0, 1)) +
    # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    geom_point(aes(q, y = pow, shape = methods, col = methods),stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)")  + ylab("Power (%)")
  p2 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  #saveRDS(p, file='lineplot_DEanalysis_pow.rds')
  # grid.newpage()
  # grid.draw(set_panel_size(p, width  = unit(2, "in"),
  #                          height = unit(2/0.11*transformy(1), "in")))
  # dev.off()
  
  p = ggplot(dat, 
             aes(x = q, y = pow.i, group = methods)) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.position = 'none',
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_x_continuous(limits = c(0, 0.11),
                       breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                       expand = c(0,0),
                       labels = function(x) paste0(x*100)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = c(0.25, 0.5, 0.75),
                       expand = c(0,0),
                       labels = c(0.25, 0.5, 0.75)*100) + 
    coord_cartesian( ylim=c(0, 1)) +
    # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
    geom_line(aes(q, y = pow.i, col = methods), lwd = lwd) +
    scale_color_manual(values = col_palette) +
    geom_point(aes(q, y = pow.i, shape = methods, col = methods),stroke = 1,size = pt_size) +
    scale_shape_manual(values = pt_shape) + xlab("Actual FDR (%)") + ylab("Power (%)")
  
  p = p + theme(legend.position = 'right',
                legend.text = element_text(size = textsize),
                legend.key.size = unit(0.4,'in'))
  legend = get_legend(p)
  p = p + theme(legend.position = 'none')
  p3 = set_panel_size(p,width  = penal_size_width, height = penal_size_height)
  
  p_ls = list(p1, p2, p3)
  return(list(p_ls,legend))
}

pt_shape2 <- c(1,3)
names(pt_shape2) <- c("non-DEG", "true DEG")
dat2plot2 <- function(dat, logfd){
  is.trueDE <-rep("non-DEG", length(logfd))
  is.trueDE[trueDE_idx] = "true DEG"
  dat0 <- cbind(-log(dat[[1]]$pval$edger_norm), abs(logfd))
  dat0[is.na(dat0[,1]),1] = 0
  df <- data.frame(edgeR=dat0[,1], logFC =dat0[,2],
                   edgeRrank=rank(-dat0[,1]),
                   is.trueDE= factor(is.trueDE))
  df <- subset(df, edgeRrank<=100)
  df$logFCrank <- rank(-df$logFC)
  df$logFC[df$logFC == 0] <- rnorm(sum(df$logFC == 0),mean = 0, sd = 1)
  df <- subset(df, logFC!=0)
  df$edgeR <- scale(df$edgeR)
  df$logFC <- scale(df$logFC)
  res <- est.IDR(scale(df[,c(1,2)]), mu=0, sigma=1, rho=0.5, p=0.1)
  df$IDR <-res$IDR
  
  p1 = ggplot(df) + geom_point(aes(x = edgeRrank,y = logFCrank, color =IDR, shape = is.trueDE), size = 2)  +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_shape_manual(values = pt_shape2) + scale_color_gradient(low = "blue", high = "red", breaks = c(0,0.2, 0.4,0.6, 0.8), limits = range(0, 0.8)) +
    ggtitle(paste0("Cor = ", round(cor(df$edgeRrank, df$logFCrank),3))) +
    xlab("edgeR rank") +
    ylab("logFC rank") +
    labs(shape="Gene", color="IDR") + guides(shape = guide_legend(order = 1))
  p1 = set_panel_size(p1,width  = penal_size_width, height = penal_size_height2)
  
  dat0 <- cbind(-log(dat[[1]]$pval$deseq2), abs(logfd))
  dat0[is.na(dat0[,1]),1] = 0
  df <- data.frame(DESeq2=dat0[,1], logFC = dat0[,2],
                   DESeq2rank=rank(-dat0[,1]),
                   is.trueDE= factor(is.trueDE))
  df <- subset(df, DESeq2rank<=100)
  df$logFCrank <- rank(-df$logFC)
  df <- subset(df, logFC!=0)
  df$DESeq2 <- scale(df$DESeq2)
  df$logFC <- scale(df$logFC)
  res <- est.IDR(scale(df[,c(1,2)]), mu=0, sigma=1, rho=0.5, p=0.1)
  
  df$IDR <-res$IDR
  
  p2 = ggplot(df) + geom_point(aes(x = DESeq2rank,y = logFCrank, color =IDR, shape = is.trueDE), size = 2) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA))+
    scale_shape_manual(values = pt_shape2)+ 
    scale_color_gradient(low = "blue", high = "red", breaks = c(0,0.2, 0.4), limits = c(0, 0.5)) + 
    ggtitle(paste0("Cor = ", round(cor(df$DESeq2rank, df$logFCrank),3)))+
    xlab("DESeq2 rank") +
    ylab("logFC rank")  +
    labs(shape="Gene", color="IDR") + guides(shape = guide_legend(order = 1))
  p2 = set_panel_size(p2,width  = penal_size_width, height = penal_size_height2)
  
  
  dat0 <- cbind(-dat[[1]]$pval$clipper_edgernorm_diff_max_9, abs(logfd))
  df <- data.frame(Clipper=dat0[,1], logFC = dat0[,2],
                   Clipperrank=rank(-dat0[,1]),
                   is.trueDE= factor(is.trueDE))
  df <- subset(df, Clipperrank<=100)
  df$logFCrank <- rank(-df$logFC)
  df <- subset(df, logFC!=0)
  df$Clipper <- scale(df$Clipper)
  df$logFC <- scale(df$logFC)
  res <- est.IDR(scale(df[,c(1,2)]), mu=0, sigma=1, rho=0.5, p=0.1)
  
  df$IDR <-res$IDR
  
  p3 = ggplot(df) + geom_point(aes(x = Clipperrank,y = logFCrank, color =IDR, shape = is.trueDE), size = 2)  +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA))+
    scale_shape_manual(values = pt_shape2) + 
    scale_color_gradient(low = "blue", high = "red", breaks = c(0,0.2, 0.4), limits = c(0, 0.5))+
    ggtitle(paste0("Cor = ", round(cor(df$Clipperrank, df$logFCrank),3)))+
    xlab("Clipper rank") +
    ylab("logFC rank")  +
    labs(shape="Gene", color="IDR") + guides(shape = guide_legend(order = 1))
  p3 = set_panel_size(p3,width  = penal_size_width, height = penal_size_height2)
  
  p_ls = list(p1, p2, p3)
  
}

pt_shape2 <- c(1,3)
names(pt_shape2) <- c("non-DEG", "true DEG")
dat2plot3 <- function(dat, logfd){
  dat0 <- cbind(-log(dat[[1]]$pval$edger_norm), abs(logfd))
  dat0[is.na(dat0[,1]),1] = 0
  df <- data.frame(edgeR=dat0[,1], logFC =dat0[,2],
                   edgeRrank=rank(-dat0[,1]))
  df <- subset(df, edgeRrank<=100)
  df$logFCrank <- rank(-df$logFC)
  df <- subset(df, logFC!=0)
  df <- subset(df, !is.nan(df$edgeR))

  p1 = ggplot(df) + geom_point(aes(x = edgeRrank,y = logFCrank), size = 2)  +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA)) +
    scale_shape_manual(values = pt_shape2) + scale_color_gradientn(colors = c("blue", "white", "red"), values=rescale(c(0,0.2,1)),
                                                                   limits=c(0,1), breaks = c(0,0.2, 0.4,0.6, 0.8)) +
    ggtitle(paste0("Pearson Cor = ", round(cor(df$edgeRrank, df$logFCrank, method = "pearson"),3), '\n',
                  'Spearman Cor = ',  round(cor(df$edgeRrank, df$logFCrank, method = "spearman"),3))) +
    xlab("edgeR p-value rank") +
    ylab("logFC rank") +
    labs(shape="Gene", color="IDR") + guides(shape = guide_legend(order = 1))
  p1 = set_panel_size(p1,width  = penal_size_width, height = penal_size_height2)
  
  dat0 <- cbind(-log(dat[[1]]$pval$deseq2_norm), abs(logfd))
  dat0[is.na(dat0[,1]),1] = 0
  df <- data.frame(DESeq2=dat0[,1], logFC = dat0[,2],
                   DESeq2rank=rank(-dat0[,1]))
  df <- subset(df, DESeq2rank<=100)
  df$logFCrank <- rank(-df$logFC)
  df <- subset(df, logFC!=0)
  df <- subset(df, !is.infinite(df$DESeq2))

  p2 = ggplot(df) + geom_point(aes(x = DESeq2rank,y = logFCrank), size = 2) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA))+
    scale_shape_manual(values = pt_shape2)+ 
    scale_color_gradientn(colors = c("blue", "white", "red"), values=rescale(c(0,0.2,1)),
                          limits=c(0,1), breaks = c(0,0.2, 0.4,0.6, 0.8)) +
    ggtitle(paste0("Pearson Cor = ", round(cor(df$DESeq2rank, df$logFCrank, method = "pearson"),3), '\n',
                   'Spearman Cor = ',  round(cor(df$DESeq2rank, df$logFCrank, method = "spearman"),3))) +
    
    xlab("DESeq2 p-value rank") +
    ylab("logFC rank")  +
    labs(shape="Gene", color="IDR") + guides(shape = guide_legend(order = 1))
  p2 = set_panel_size(p2,width  = penal_size_width, height = penal_size_height2)
  
  
  dat0 <- cbind(-dat[[2]]$pval$clipper_edgernorm_diff_max_9, abs(logfd))
  df <- data.frame(Clipper=dat0[,1], logFC = dat0[,2],
                   Clipperrank=rank(-dat0[,1]))
  df <- subset(df, Clipperrank<=100)
  df$logFCrank <- rank(-df$logFC)
  df <- subset(df, logFC!=0)

  p3 = ggplot(df) + geom_point(aes(x = Clipperrank,y = logFCrank), size = 2)  +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA))+
    scale_shape_manual(values = pt_shape2) + 
    scale_color_gradientn(colors = c("blue", "white", "red"), values=rescale(c(0,0.2,1)),
                          limits=c(0,1), breaks = c(0,0.2, 0.4,0.6, 0.8)) +
    ggtitle(paste0("Pearson Cor = ", round(cor(df$Clipperrank, df$logFCrank, method = "pearson"),3), '\n',
                   'Spearman Cor = ',  round(cor(df$Clipperrank, df$logFCrank, method = "spearman"),3))) +
    
    xlab("Clipper contrast score rank") +
    ylab("logFC rank")  +
    labs(shape="Gene", color="IDR") + guides(shape = guide_legend(order = 1))
  p3 = set_panel_size(p3,width  = penal_size_width, height = penal_size_height2)
  
  p_ls = list(p1, p2, p3)
  
}

cs2q <- function(cs){
  fdp = sapply(cs, function(x){
    if (x <= 0){
      q = 1
    }else{
      q = (1 + sum(cs<= -x))/(9*sum(cs>=x))
    }
    return(q)
  })
  x <- sort(unique(cs))
  qvalue <- rep(0, length(x))
  for (i in 1:length(x)){
    qvalue[i] <- min(fdp[cs <= x[i]])
  }
  qvalue <- qvalue[match(cs, x)]
}
dat2plot4 <- function(dat, logfd, trueDE_idx){
  is.trueDE <-rep("non-DEG", length(logfd))
  is.trueDE[trueDE_idx] = "true DEG"
  
  dat0 <- cbind(-log(dat[[1]]$pval$edger_norm), -log(dat[[2]]$pval$edger_norm))
  dat0[is.na(dat0[,1]),1] = 0
  dat0[is.na(dat0[,2]),2] = 0
  df <- as.data.frame(dat0)
  names(df) <- c("Exp1", "Exp2")
  df$trueDE <- is.trueDE
  df$rank1 <- rank(-df$Exp1)
  df$rank2 <- rank(-df$Exp2)
  df <- subset(df, rank1 <= 100|rank2 <= 100)
  df <- subset(df, !is.infinite(Exp1) & !is.infinite(Exp2))
  res <- est.IDR(df[,c(1,2)], mu=4, sigma=1, rho=0.5, p=0.2)
  df$IDR <- res$IDR
  
  llim = min(c(df$Exp1, df$Exp2))
  ulim = max(c(df$Exp1, df$Exp2))
  p1 = ggplot(df) + geom_point(aes(x = Exp1,y = Exp2, color =IDR), size = 2) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA)) + scale_color_gradientn(colors = c("blue", "white", "red"), values=rescale(c(0,0.2,1)),
                                                                   limits=c(0,1), breaks = c(0,0.2, 0.4,0.6, 0.8)) +
    ggtitle(paste0("edgeR\nPearson Cor = ", round(cor(df$Exp1, df$Exp2, method = "pearson"),3), '\n',
                   'Spearman Cor = ',  round(cor(df$Exp1, df$Exp2, method = "spearman"),3))) +
    
    xlab("Dataset 1 -log(p-value)") + xlim(c(llim,ulim)) + 
    ylab("Dataset 2 -log(p-value)") + ylim(c(llim,ulim))
    labs(color="IDR")
  p1 = set_panel_size(p1,width  = penal_size_width, height = penal_size_height2)
  
  dat0 <- cbind(-log(dat[[1]]$pval$deseq2_norm), -log(dat[[2]]$pval$deseq2_norm))
  dat0[is.na(dat0[,1]),1] = 0
  dat0[is.na(dat0[,2]),2] = 0
  df <- as.data.frame(dat0)
  names(df) <- c("Exp1", "Exp2")
  df$trueDE <- is.trueDE
  df$rank1 <- rank(-df$Exp1)
  df$rank2 <- rank(-df$Exp2)
  df <- subset(df, rank1 <= 100|rank2 <= 100)
  df <- subset(df, !is.infinite(Exp1) & !is.infinite(Exp2))
  res <- est.IDR(df[,c(1,2)], mu=4, sigma=1, rho=0.5, p=0.2)
  df$IDR <- res$IDR
  
  llim = min(c(df$Exp1, df$Exp2))
  ulim = max(c(df$Exp1, df$Exp2))
  
  p2 = ggplot(df) + geom_point(aes(x = Exp1,y = Exp2, color =IDR), size = 2) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          legend.position = 'none',
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA)) + scale_color_gradientn(colors = c("blue", "white", "red"), values=rescale(c(0,0.2,1)),
                                                                                     limits=c(0,1), breaks = c(0,0.2, 0.4,0.6, 0.8)) +
    ggtitle(paste0("DESeq2\nPearson Cor = ", round(cor(df$Exp1, df$Exp2, method = "pearson"),3), '\n',
                   'Spearman Cor = ',  round(cor(df$Exp1, df$Exp2, method = "spearman"),3))) +
    xlab("Dataset 1 -log(p-value)") +xlim(c(llim,ulim)) + 
    ylab("Dataset 2 -log(p-value)") + ylim(c(llim,ulim))
    labs(color="IDR")
  p2 = set_panel_size(p2,width  = penal_size_width, height = penal_size_height2)
  
  
  dat0 <- cbind(-dat[[1]]$pval$clipper_edgernorm_diff_max_9, -dat[[2]]$pval$clipper_edgernorm_diff_max_9)
  dat0[is.na(dat0[,1]),1] = 0
  dat0[is.na(dat0[,2]),2] = 0
  
  df <- as.data.frame(dat0)
  names(df) <- c("Exp1", "Exp2")
  # df$q1 <- -log(cs2q(df$Exp1))
  # df$q2 <- -log(cs2q(df$Exp2))
  df$rank1 <- rank(-df$Exp1)
  df$trueDE = is.trueDE
  df$rank2 <- rank(-df$Exp2)
  df <- subset(df, rank1 <= 100|rank2 <= 100)
  df <- subset(df, !is.infinite(Exp1) & !is.infinite(Exp2))
  #df$Exp1[df$Exp1 <0] = 0
  #df$Exp2[df$Exp2 <0] = 0
  res <- est.IDR(df[,c(1,2)], mu=3, sigma=1, rho=0.5, p=0.2)
  df$IDR <- res$IDR
  llim = min(c(df$Exp1, df$Exp2))
  ulim = max(c(df$Exp1, df$Exp2))
  
  
  p3 = ggplot(df) + geom_point(aes(x = Exp1,y = Exp2, color =IDR), size = 2) +
    theme(text = element_text(size=textsize, color = 'black'),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.y = element_line(size = 0.5),
          axis.line.x = element_line(size = 0.5),
          # axis.title.x = element_blank(),
          axis.text = element_text( colour = 1, size = textsize),
          legend.key = element_rect(colour = NA, fill = NA)) + scale_color_gradientn(colors = c("blue", "white", "red"), values=rescale(c(0,0.2,1)),
                                                                                     limits=c(0,1), breaks = c(0,0.2, 0.4,0.6, 0.8)) +
    ggtitle(paste0("Clipper\nPearson Cor = ", round(cor(df$Exp1, df$Exp2, method = "pearson"),3), '\n',
                   'Spearman Cor = ',  round(cor(df$Exp1, df$Exp2, method = "spearman"),3))) +
    xlab("Dataset 1 contrast score") + xlim(c(llim,ulim))
    ylab("Dataset 2 contrast score") + ylim(c(llim,ulim))
    labs(color="IDR")
  legend = get_legend(p3)
  p3 = p3 + theme(legend.position = 'none')
  p3 = set_panel_size(p3,width  = penal_size_width, height = penal_size_height2)
  p_ls = list(p1, p2, p3)
  return(list(p_ls,legend))
  
  
  p_ls = list(p1, p2, p3)
  
}

penal_size_width  = unit(2, "in")
penal_size_height2 = unit(2, "in")
FDR = seq(0.01, 0.1 ,0.01)




blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

dat2plot_combined <- function(dat, trueDE_ls, filename, title){
  p_ls <- dat2plot1(dat, use.norm = 1)
  trueDE_idx <- trueDE_ls[[2]]
  logfd <- log(base = 2, trueDE_ls[[1]])
  true.logdf <- abs(log(base = 2, trueDE_ls[[1]]))[trueDE_idx]
  
  p_ls2 <- dat2plot3(dat, logfd = logfd)
  
  p_ls3 <- dat2plot4(dat, logfd = logfd, trueDE_idx = trueDE_idx)
  pdf(file = paste0(filename, '.pdf'), height = 15, width = 15)
  p = arrangeGrob(ncol = 4,nrow =9, widths = c(4,4,4, 2), heights = c(0.4, 0.2, 5, 0.2, 0.2, 4,0.2, 0.2, 4),
                  blankPlot, textGrob(label = title, gp = gpar(fontsize= titletextsize, fontface =2)), blankPlot,blankPlot,
                  textGrob(label = 'a', hjust = 7, gp = gpar(fontsize= titletextsize, fontface =2)), blankPlot,blankPlot,blankPlot, 
                  p_ls[[1]][[1]], p_ls[[1]][[2]], p_ls[[1]][[3]], p_ls[[2]],
                  textGrob(label = 'b', hjust = 7, gp = gpar(fontsize= titletextsize, fontface =2)),blankPlot,blankPlot,blankPlot,
                  blankPlot,textGrob(label = 'True DEGs among top 100 DEGs identified by each DE method', gp = gpar(fontsize= titletextsize, fontface =2)),blankPlot,blankPlot,
                  p_ls2[[1]], p_ls2[[2]],  p_ls2[[3]],blankPlot,
                  textGrob(label = 'c', hjust = 7, gp = gpar(fontsize= titletextsize, fontface =2)),blankPlot,blankPlot,blankPlot,
                  blankPlot,textGrob(label = 'True DEGs identified among the top 100 DEGs in dataset 1 or dataset 2', gp = gpar(fontsize= titletextsize, fontface =2)),blankPlot,blankPlot,
                  p_ls3[[1]][[1]], p_ls3[[1]][[2]], p_ls3[[1]][[3]],p_ls3[[2]]
  )
  grid.arrange(p,  ncol = 1,nrow = 1)
  dev.off()
  return(0)
}

dat = readRDS(file = 'DEG_spikeinFollow2019_wbatcheffect_17vs17data_070321_fd16.rds')
trueDE_ls <- readRDS("DEG_spikeinFollow2019_wbatcheffect_17vs17data_070321_trueFD&DEidx.rds")
dat2plot_combined(dat, trueDE_ls, filename = "17vs17_2019", title = 'Human monocyte semi-synthetic datasets by simulation strategy 2')
dat = readRDS(file = 'DEG_spikeinFollow2019_wbatcheffect_48vs48data_070321_fd3divby2.rds')
trueDE_ls <- readRDS("DEG_spikeinFollow2019_wbatcheffect_48vs48data_070321_trueFD&DEidx.rds")
dat2plot_combined(dat, trueDE_ls, filename = "48vs48_2019", title = 'Yeast semi-synthetic datasets by simulation strategy 2')
dat = readRDS(file = 'DEG_spikeinNDEGbypermute_17vs17data_070321.rds')
trueDE_ls <- readRDS("DEG_spikeinNDEGbypermute_17vs17data_070321_trueFD&DEidx.rds")
dat2plot_combined(dat, trueDE_ls, filename = "17vs17_permute", title = 'Human monocyte semi-synthetic datasets by simulation strategy 1')
dat = readRDS(file = 'DEG_spikeinNDEGbypermute_48vs48data_070321.rds')
trueDE_ls <- readRDS("DEG_spikeinNDEGbypermute_48vs48data_070321_trueFD&DEidx.rds")
dat2plot_combined(dat, trueDE_ls, filename = "48vs48_permute", title = 'Yeast semi-synthetic datasets by simulation strategy 1')
 
