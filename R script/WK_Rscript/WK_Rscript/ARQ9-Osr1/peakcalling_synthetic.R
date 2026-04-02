source("../../codes_figures/clipper042520_w2sided.R")
combine.peaks <- function(peaks, gap) {
  i <- 1
  while(i < nrow(peaks)) {
    if(peaks[i+1,1] - peaks[i,2]<= gap) {
      peaks[i,2] <- peaks[i+1, 2] 
      peaks[i,3] <- max(peaks[i+1, 3], peaks[i,3])
      peaks <- peaks[-i,]
    }else{
      i = i+1
    }
  }
  return(peaks)
}
merge <- function (x,y) {
  nx <- nrow(x)
  ny <- nrow(y)
  nm <- nx + ny
  m <- matrix(0,nm,4)
  x.p <- 1
  y.p <- 1
  m.p <- 1
  comp <- c(x[1,2], y[1,2])
  start <- 0
  end <- min (comp)
  c <- which.min(comp)
  if (comp[1] == comp[2]) {
    c = 0
  }
  m[m.p,] <- c(start, end,x[x.p,3], y[y.p,3])
  m.p <- m.p +1
  if (c==1) {
    x.p  = x.p + 1
  }
  if ( c==2) {
    y.p = y.p +1
  }
  if (c == 0) {
    x.p = x.p + 1
    y.p = y.p + 1
  }
  while (x.p <= nx & y.p <= ny) {
    comp <- c(x[x.p,2], y[y.p,2])
    start = end
    end =  min (comp)
    c <- which.min(comp)
    if (comp[1] == comp[2]) {
      c = 0
    }
    m[m.p,] <- c(start, end,x[x.p,3], y[y.p,3])
    m.p <- m.p +1
    if (c==1) {
      x.p  = x.p + 1
    }
    if ( c==2) {
      y.p = y.p +1
    }
    if (c == 0) {
      x.p = x.p + 1
      y.p = y.p + 1
    }
  }
  if (x.p > nx & y.p <= ny) {
    start = end
    end = y[y.p,2]
    m[m.p, ] <- c(start, end, 0, y[y.p,3])
    m.p <- m.p +1
    y.p = y.p + 1
    while (y.p <= ny) {
      start = end
      end = y[y.p,2]
      m[m.p, ] <- c(start, end, 0, y[y.p,3])
      m.p <- m.p +1
      y.p = y.p + 1
    }
  }
  if (y.p > ny & x.p <= nx) {
    start = end
    end = x[x.p,2]
    m[m.p,] <- c(start, end, x[x.p,3], 0)
    m.p <- m.p +1
    x.p = x.p + 1
    while (x.p <= nx) {
      start = end
      end = x[x.p,2]
      m[m.p,] <- c(start, end, x[x.p,3], 0)
      m.p <- m.p +1
      x.p = x.p +1
    }
  }
  colnames(m) <- c("start", "end", "x", "y")
  return(m)
}
#this function combines continuous clipper outcomes into one
clip2reg <- function (re, m, i) {
  dis <- re$results[[i]]$discovery
  n <- length(dis)
  start <- rep(0,n)
  end <- rep(0, n)
  score <- rep(0,n)
  k <- 1
  p <- 1
  while (k <= n ) {
    start[p] <- dis[k]
    score[p] <- re$contrastScore[dis[k]]
    while (k < n & dis[k+1] <= dis[k] + 1) {
      k = k+1
      if( re$contrastScore[dis[k]] > score[p]) {
        score[p] <- re$contrastScore[dis[k]]
      }
    }
    end[p] <- dis[k]
    
    p = p + 1
    k = k + 1
  }
  start <- start[start!=0]
  end <- end[end!=0]
  score <- score[score!=0]
  return(data.frame(start = m[start,1], end = m[end,2], contrastscore = score))
}

combine.peaks <- function(peaks, gap) {
  i <- 1
  while(i < nrow(peaks)) {
    if(peaks[i+1,1] - peaks[i,2]<= gap) {
      peaks[i,2] <- peaks[i+1, 2] 
      peaks[i,3] <- max(peaks[i+1, 3], peaks[i,3])
      peaks <- peaks[-i,]
    }else{
      i = i+1
    }
  }
  return(peaks)
}

fdp_regions <-function(peak, truth) {
  over.peak <- sapply(1:nrow(peak), function(i){
    r1 <- as.numeric(peak[i,1:2])
    return(sum(!(r1[1]>truth[,2]|r1[2]<truth[,1])))
  })
  return((nrow(peak)-sum(over.peak>0))/nrow(peak))
}

pow_regions <-function(peak, truth) {
  over.peak <- sapply(1:nrow(truth), function(i){
    r1 <- as.numeric(truth[i,])
    return(sum(!(r1[1]>peak[,2]|r1[2]<peak[,1])))
  })
  return(sum(over.peak>0)/nrow(truth))
}
loverlape <- function(as, ae, bs, be){
  return((!(as > be|bs > ae))&(pmin(ae, be) - pmax(as,bs))>= 1/2 * (as-ae))
}
fdp_regions2 <-function(peak, truth) {
  over.peak <- sapply(1:nrow(peak), function(i){
    r1 <- as.numeric(peak[i,1:2])
    return(sum(loverlape(r1[1], r1[2], truth[,1], truth[,2])))
  })
  return((nrow(peak)-sum(over.peak>0))/nrow(peak))
}

pow_regions2 <-function(peak, truth) {
  over.peak <- sapply(1:nrow(truth), function(i){
    r1 <- as.numeric(truth[i,])
    return(sum(loverlape(peak[,1], peak[,2], r1[1], r1[2])))
  })
  return(sum(over.peak>0)/nrow(truth))
}
fdppow2_regions <- function(peak, truth, score, FDR){
  over.peak <- sapply(1:nrow(peak), function(i){
    r1 <- as.numeric(peak[i,1:2])
    return(sum(loverlape(r1[1], r1[2], truth[,1], truth[,2])))
  })
  n.dis <- cumsum(over.peak[order(score)])
  FDR.i = 1 - n.dis/1:length(score)
  peak0 <- peak[order(score),]
  pow <- sapply(FDR, function(q){
    index.i <- max(which(FDR.i<=q))
    peak1 <- peak0[1:index.i,]
    over.peak2 <- sapply(1:nrow(truth), function(i){
      r1 <- as.numeric(truth[i,])
      return(sum(!(r1[1]>peak1[,2]|r1[2]<peak1[,1])))
    })
    return(sum(over.peak2>0)/nrow(truth))
  })
  return(pow)
}

compute_fdppow2 = function(pvalues, trueidx, q){
  index.p <- which(!is.na(pvalues))
  is.trueDE <- 1:length(pvalues) %in% trueidx
  n.dis <- cumsum(is.trueDE[order(pvalues)])
  FDR.i = 1 - n.dis/1:length(pvalues)
  index.i <- max(which(FDR.i<=q))
  fdp <- FDR.i[index.i]
  pow <- n.dis[index.i]/length(trueidx)
  return(c(fdp, pow))
}
promotor <- read.csv("mart_promotor.txt")
colnames(promotor) <- c("chr", "start", "end", "type")
chrom = "1"
prom.chr <- subset(promotor, chr == chrom)
prom.chr <- prom.chr[,2:3]
prop_prom <- 0.2
prom.chr <- prom.chr[prom.chr$end - prom.chr$start <10000 & prom.chr$end - prom.chr$start>=300, ]

for (j in 1:20){
  sample_prom <- sample(1:nrow(prom.chr), prop_prom*nrow(prom.chr))
  setwd("peakcalling/")
  dat <- data.frame(chr = rep("chr1", length(sample_prom)), start = prom.chr$start[sample_prom],
                    end = prom.chr$end[sample_prom], p_ext = runif(length(sample_prom), 0.6, 1),
                    p_amp	= runif(length(sample_prom), 0.6, 1), energy_A = runif(length(sample_prom), 0, 1))
  write.table(dat, file=paste0("example/rep", j, "/test.tsv"), quote=FALSE, sep='\t', row.names = F)
  write.table(sample_prom, paste0("example/rep", j, "/sample_prom.txt"), quote = F, row.names = F, col.names = F)
}
system(paste0("python3 chipulate-master/chipulate.py --input-file test.tsv --genome-file chr1.fa --chrom-size-file chr1.fa.fai -d 50 --output-dir example/rep", j))
system(paste0("bedops -u example/rep", j, "/test.chip_reads.bed syn.bed > example/syn_with_peak", j, ".bed"))
system(paste0("bedops -u example/rep", j, "/test.control_reads.bed syn.bed > example/syn2_with_peak", j, ".bed"))
system(paste0("macs2 callpeak -t example/syn_with_peak",j, ".bed example/syn2_with_peak",j, ".bed -c 1.bed -f BED -n twosample -B -q 0.5 --outdir example/rep",j))
system(paste0("macs2 callpeak -t example/syn_with_peak",j, ".bed -f BED -n experimental -B -q 1 --outdir example/rep",j))
system(paste0("macs2 callpeak -t example/syn2_with_peak",j, ".bed -f BED -n experimental2 -B -q 1 --outdir example/rep",j))
system(paste0("macs2 callpeak -t 1.bed -f BED -n control -B -q 1 --outdir example/rep",j))


syn.results <- sapply(1:10, function(i) {
  fdppow.i <- sapply(2:20, function(j){
    macs2.syn.peak <- read.table(paste0("example/rep", j, "/twosample_peaks.narrowPeak"))
    sample_prom <-  read.table(paste0("example/rep", j, "/sample_prom.txt"))$V1
    macs2.syn.peak <- macs2.syn.peak[macs2.syn.peak[,9]>=-log(i/100, 10),]
    macs2.syn.peak <- macs2.syn.peak[order(macs2.syn.peak$V2),c(2,3,9)]
    macs2.syn.peak <- combine.peaks(macs2.syn.peak, 50)
    fdp.macs <- fdp_regions(macs2.syn.peak, prom.chr[sample_prom,])
    pow.macs <- pow_regions(macs2.syn.peak, prom.chr[sample_prom,])
    return(c(fdp.macs, pow.macs))
  })
  return(c(mean(fdppow.i[1,]), mean(fdppow.i[2,]), var(fdppow.i[1,]), var(fdppow.i[2,])))
})

saveRDS(syn.results, "macs_fdppow.rds")

clipper_fdppow <- lapply(1:20, function(j){
  control <- read.table(paste0("example/rep", j, "/control_treat_pileup.bdg"))
  experimental <-  read.table(paste0("example/rep", j, "/experimental_treat_pileup.bdg"))
  sample_prom <-  read.table(paste0("example/rep", j, "/sample_prom.txt"))$V1
  m <- merge(experimental[,2:4], control[,2:4])
  m <- m[1:(min(which(m[,1] + m[,2] == 0))-1),]
  re.diff <- clipper(m[,3], m[,4],  FDR = (1:10)/100)
  ls_fdppow <- sapply(1:10, function(i){
    reg.syn <- clip2reg(re.diff, m, i)
    reg.syn.sub <- reg.syn[reg.syn[,2]-reg.syn[,1]>=200,]
    over.clip.peak <- sapply(1:nrow(reg.syn.sub), function(i){
      r1 <- as.numeric(reg.syn.sub[i,])
      return(sum(!(r1[1]>prom.chr[sample_prom,2]|r1[2]<prom.chr[sample_prom,1])))
    })
    fdp.clipper <- (nrow(reg.syn.sub)-sum(over.clip.peak>0))/nrow(reg.syn.sub)
    
    over.peak.clip <- sapply(1:length(sample_prom), function(i){
      r1 <- as.numeric(prom.chr[sample_prom[i],])
      return(sum(!(r1[1]-500>reg.syn.sub[,2]|r1[2]+500<reg.syn.sub[,1])))
    })
    pow.clipper <- sum(over.peak.clip>0)/length(sample_prom)
    return(c(fdp.clipper, pow.clipper))
  })
  return(ls_fdppow)
})

fdp = sapply(2:20, function(j){clipper_fdppow[[j]][1,]})
pow = sapply(2:20, function(j){clipper_fdppow[[j]][2,]})
clipper_fdppow2 <- rbind(apply(fdp,1,mean), apply(pow,1,mean), apply(fdp,1,var), apply(pow,1,var))

saveRDS(clipper_fdppow2, "clipper_fdppow.rds")
thre.i <- read.table("sam/thre_i.txt")
combine.results <- lapply(2:20, function(j){
    macs2.syn.peak <- read.table(paste0("example/rep", j, "/twosample_peaks.narrowPeak"))
    sample_prom <-  read.table(paste0("example/rep", j, "/sample_prom.txt"))$V1
    control <- read.table(paste0("example/rep", j, "/control_treat_pileup.bdg"))
    experimental <-  read.table(paste0("example/rep", j, "/experimental_treat_pileup.bdg"))
    experimental2 <-  read.table(paste0("example/rep", j, "/experimental2_treat_pileup.bdg"))
    s1 <- rep(experimental$V4, experimental$V3- experimental$V2 )
    s2 <- rep(control$V4, control$V3- control$V2)
    s3 <- rep(experimental2$V4, experimental2$V3- experimental2$V2 )
    s1[(length(s1)+1):length(s2)] <- 0
    s3[(length(s3)+1):length(s2)] <- 0
    p1 = (s1+s3)/2 - s2
    p2 = (s1 + s2)/2 - s1
    #re_clipper <- Clipper(score.exp = p1, score.back = p2, analysis = "e", FDR = seq(0.01, 0.1, 0.01))
    # median.peak <- sapply(1:nrow(macs2.syn.peak), function(i){
    #   # start point of the candidate peak
    #   start <- macs2.syn.peak$V2[i]+1
    #   # end point of the candiate peak
    #   end <- macs2.syn.peak$V3[i]
    #   # the different contrast score is the difference between the two condiitons
    #   return(median(1/2*(s1[start:end] + s3[start:end]) -(s2[start:end])))
    # })
    p3 = p1 - p2
    mean.peak <- sapply(1:nrow(macs2.syn.peak), function(i){
      # start point of the candidate peak
      start <- macs2.syn.peak$V2[i]+1
      # end point of the candiate peak
      end <- macs2.syn.peak$V3[i]
      # the different contrast score is the difference between the two condiitons
      seg <- p3[start:end]
      return(ifelse(end-start<300, max(seg), mean(seg[seg>8])))
      #return(mean(p3[start:(end-100)]))
    })
    mean.peak[is.nan(mean.peak)] = 0
    max.peak <- sapply(1:nrow(macs2.syn.peak), function(i){
      # start point of the candidate peak
      start <- macs2.syn.peak$V2[i]+1
      # end point of the candiate peak
      end <- macs2.syn.peak$V3[i]
      # the different contrast score is the difference between the two condiitons
      return(max(p3[start:end]))
    })
    median.peak <- sapply(1:nrow(macs2.syn.peak), function(i){
      # start point of the candidate peak
      start <- macs2.syn.peak$V2[i]+1
      # end point of the candiate peak
      end <- macs2.syn.peak$V3[i]
      # the different contrast score is the difference between the two condiitons
      return(median(p3[start:end]))
    })
    # mean.peak <- sapply(1:nrow(macs2.syn.peak), function(i){
    #   # start point of the candidate peak
    #   start <- macs2.syn.peak$V2[i]+1
    #   # end point of the candiate peak
    #   end <- macs2.syn.peak$V3[i]
    #   # the different contrast score is the difference between the two condiitons
    #   return(mean(1/2*(s1[start:end] + s3[start:end]) -(s2[start:end])))
    # })
    
    ls_fdppow_mean <- sapply(1:10, function(i){
      #macs2.syn.peak2 <- macs2.syn.peak[median.peak >= thre.i[5,i],]
      macs2.syn.peak2 <- macs2.syn.peak[mean.peak >= thre.i[i],]
      #macs2.syn.peak2 <- macs2.syn.peak[mean.peak >= thre.i[5,i]-1,]
      macs2.syn.peak2 <- macs2.syn.peak2[order(macs2.syn.peak2$V2),c(2,3,9)]
      macs2.syn.peak2 <- combine.peaks(macs2.syn.peak2, 50)
      fdp.macs.median <- fdp_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      pow.macs.median <- pow_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      return(c(fdp.macs.median, pow.macs.median))
    })
    ls_fdppow_max <- sapply(1:10, function(i){
      #macs2.syn.peak2 <- macs2.syn.peak[median.peak >= thre.i[5,i],]
      macs2.syn.peak2 <- macs2.syn.peak[max.peak >= thre.i[i],]
      #macs2.syn.peak2 <- macs2.syn.peak[mean.peak >= thre.i[5,i]-1,]
      macs2.syn.peak2 <- macs2.syn.peak2[order(macs2.syn.peak2$V2),c(2,3,9)]
      macs2.syn.peak2 <- combine.peaks(macs2.syn.peak2, 50)
      fdp.macs.median <- fdp_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      pow.macs.median <- pow_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      return(c(fdp.macs.median, pow.macs.median))
    })
    ls_fdppow_median <- sapply(1:10, function(i){
      #macs2.syn.peak2 <- macs2.syn.peak[median.peak >= thre.i[5,i],]
      macs2.syn.peak2 <- macs2.syn.peak[median.peak >= thre.i[i],]
      #macs2.syn.peak2 <- macs2.syn.peak[mean.peak >= thre.i[5,i]-1,]
      macs2.syn.peak2 <- macs2.syn.peak2[order(macs2.syn.peak2$V2),c(2,3,9)]
      macs2.syn.peak2 <- combine.peaks(macs2.syn.peak2, 50)
      fdp.macs.median <- fdp_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      pow.macs.median <- pow_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      return(c(fdp.macs.median, pow.macs.median))
    })
    ls_fdppow_macs2 <- sapply(1:10, function(i){
      macs2.syn.peak2 <- macs2.syn.peak[macs2.syn.peak[,9]>=-log(i/100, 10),]
      macs2.syn.peak2 <- macs2.syn.peak2[order(macs2.syn.peak2$V2),c(2,3,9)]
      macs2.syn.peak2 <- combine.peaks(macs2.syn.peak2, 50)
      fdp.macs <- fdp_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      pow.macs <- pow_regions(macs2.syn.peak2, prom.chr[sample_prom,])
      return(c(fdp.macs, pow.macs))
    })
    
    macs2.syn.peak2 <- macs2.syn.peak[order(macs2.syn.peak$V2),c(2,3,9)]
    ls_pow2_mean <- fdppow2_regions(peak = macs2.syn.peak2, truth = prom.chr[sample_prom,], 
                                  FDR = seq(0.01, 0.1, 0.01), score = -mean.peak)
    ls_pow2_median <- fdppow2_regions(peak = macs2.syn.peak2, truth = prom.chr[sample_prom,], 
                                    FDR = seq(0.01, 0.1, 0.01), score = -median.peak)
    ls_pow2_max <- fdppow2_regions(peak = macs2.syn.peak2, truth = prom.chr[sample_prom,], 
                                    FDR = seq(0.01, 0.1, 0.01), score = -max.peak)
    ls_pow2_macs2 <- fdppow2_regions(peak = macs2.syn.peak2, truth = prom.chr[sample_prom,], 
                                   FDR = seq(0.01, 0.1, 0.01), score = -macs2.syn.peak2[,3])
  return(rbind(ls_fdppow_macs2, ls_pow2_macs2, ls_fdppow_median, ls_pow2_median, ls_fdppow_mean, ls_pow2_mean,
               ls_fdppow_max, ls_pow2_max))
})
combine.results2 <- sapply(1:nrow(combine.results[[1]]), function(i){
  return(rowMeans(sapply(1:19, function(j){combine.results[[j]][i,]})))
})

saveRDS(combine.results2, "clipper_macs_max_fdppow.rds")

for (j in 1:20){
  system(paste0("./homer/bin/makeTagDirectory example/rep", j, "/ example/syn_with_peak", j, ".bed example/syn2_with_peak", j, ".bed"))
}
for (j in 1:20){
  for (i in 1:10){
    system(paste0("./homer/bin/findPeaks example/rep", j,"/  -region -size 150 -fdr ", i/100, " -P 0.1 -LP 0.1 -o example/rep", j,"/syn_rep",j,"_", i, " -i homer/HOMER_1/"))
  }
}
homer.results <- sapply(1:10, function(i) {
  fdppow.i <- sapply(2:20, function(j){
    sample_prom <-  read.table(paste0("example/rep", j, "/sample_prom.txt"))$V1
    homer.peaks <- read.delim(paste0("example/rep", j,"/syn_rep",j,"_", i), skip = 39)
    homer.peaks <- homer.peaks[,c(3,4,6)]
    homer.peaks <- combine.peaks(homer.peaks, 50)
    fdp.homer <- fdp_regions(homer.peaks, prom.chr[sample_prom,])
    pow.homer <- pow_regions(homer.peaks, prom.chr[sample_prom,])
    return(c(fdp.homer, pow.homer))
  })
  return(c(mean(fdppow.i[1,]), mean(fdppow.i[2,]), var(fdppow.i[1,]), var(fdppow.i[2,])))
})
saveRDS(homer.results, "homer_fdppow.rds")

homer.clipper.results <- lapply(2:20, function(j){
  homer.peaks <- read.delim(paste0("example/rep", j,"/syn_rep",j,"_10"), skip = 39)
  homer.peaks <- homer.peaks[,c(3,4,6)]
  sample_prom <-  read.table(paste0("example/rep", j, "/sample_prom.txt"))$V1
  control <- read.table(paste0("example/rep", j, "/control_treat_pileup.bdg"))
  experimental <-  read.table(paste0("example/rep", j, "/experimental_treat_pileup.bdg"))
  experimental2 <-  read.table(paste0("example/rep", j, "/experimental2_treat_pileup.bdg"))
  s1 <- rep(experimental$V4, experimental$V3- experimental$V2 )
  s2 <- rep(control$V4, control$V3- control$V2)
  s3 <- rep(experimental2$V4, experimental2$V3- experimental2$V2 )
  s1[(length(s1)+1):length(s2)] <- 0
  s3[(length(s3)+1):length(s2)] <- 0
  p1 = (s1+s3)/2 - s2
  p2 = (s1 + s2)/2 - s1
  p3 = p1 - p2
  max.peak <- sapply(1:nrow(homer.peaks), function(i){
    # start point of the candidate peak
    start <- homer.peaks[i,1]+1
    # end point of the candiate peak
    end <- homer.peaks[i,2]
    # the different contrast score is the difference between the two condiitons
    return(mean(p3[start:end]))
  })
  ls_fdppow_max <- sapply(1:10, function(i){
    #macs2.syn.peak2 <- macs2.syn.peak[median.peak >= thre.i[5,i],]
    homer.peak2 <- homer.peaks[max.peak >= thre.i[i],]
    #macs2.syn.peak2 <- macs2.syn.peak[mean.peak >= thre.i[5,i]-1,]
    fdp.macs.median <- fdp_regions(homer.peak2, prom.chr[sample_prom,])
    pow.macs.median <- pow_regions(homer.peak2, prom.chr[sample_prom,])
    return(c(fdp.macs.median, pow.macs.median))
  })
  ls_pow2_max <- fdppow2_regions(peak = homer.peaks, truth = prom.chr[sample_prom,], 
                                 FDR = seq(0.01, 0.1, 0.01), score = -max.peak)
  return(rbind(ls_fdppow_max, ls_pow2_max))
})
homer.clipper.results2 <- sapply(1:nrow(homer.clipper.results[[1]]), function(i){
  return(rowMeans(sapply(1:19, function(j){homer.clipper.results[[j]][i,]})))
})
saveRDS(homer.clipper.results2, "homer_clipper_fdppow.rds")

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
method = c('Homer', 'Homer + Clipper', 'MACS2', 'Clipper', 'MACS + Clipper', 'MACS + Clipper(mean)', 'MACS + Clipper(median)')
col_clipper = rgb(0, 200, 200, max = 255) # blue
col_clipper_max =  rgb(0, 114, 178, max = 255)   # light blue
col_clipper_mean = rgb(20, 200, 150, max = 255)
col_clipper_median = '#90A93B'
col_macs2 = '#E8BA00'
col_homer = "#D69840"  # yellow
col_clipper_homer = "#FDB552"

col_palette = c(col_homer, col_clipper_homer, col_macs2, col_clipper, col_clipper_max, 
                col_clipper_mean, col_clipper_median)
names(col_palette) = method
pt_shape = c(1:7)
names(pt_shape) = method

lwd = 0.75
pt_size = 2
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
############### line plots ##################
ylim_fdr = 0.30

dat1 = readRDS('homer_fdppow.rds')
dat2 = readRDS('homer_clipper_fdppow.rds')
dat3 = readRDS('clipper_fdppow.rds')
dat4 = readRDS('clipper_macs_max_fdppow.rds')
dat4 <- rbind(dat4[,1:3], dat4[,4:6], dat4[,7:9], dat4[,10:12])

dat <- rbind(t(dat1)[,1:2], dat2[,1:2], t(dat3)[,1:2], dat4[,1:2])
dat <- as.data.frame(dat)
colnames(dat) <- c('fdr', 'pow')
dat$methods <- rep(c('Homer', 'Homer + Clipper',  'Clipper', 'MACS2', 'MACS + Clipper(median)', 
                     'MACS + Clipper(mean)', 'MACS + Clipper'), each = 10)
dat$q <- rep(seq(0.01, 0.1, 0.01), 7)

dat$fdr_tr = transformy(dat$fdr)
dat$fdr_tr[is.nan(dat$fdr_tr)] = 0
####### clipper variants #####
dat <- dat[dat$methods %in% c('Homer', 'Homer + Clipper', 'MACS2', 'MACS + Clipper'),]

pdf(file = paste0('~/Dropbox/clipper/peak calling/lineplot_peakcalling_fdr_0624.pdf'), height = 5, width = 7)
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
  geom_line(aes(q, y = fdr_tr, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = fdr_tr, shape = methods, col = methods), stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab('Actual FDR (%)')
#saveRDS(p, file='lineplot_peakcalling_fdr.rds')

grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

pdf(file = paste0('~/Dropbox/clipper/peak calling/lineplot_peakcalling_pow_0624.pdf'), height = 5, width = 7)
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
        legend.key = element_rect(colour = NA, fill = NA)) +
  scale_x_continuous(limits = c(0, 0.11),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                     expand = c(0,0),
                     labels = function(x) paste0(x*100)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0.25, 0.5, 0.75, 1),
                     expand = c(0,0),
                     labels = c(0.25, 0.5, 0.75, 1)*100) + 
  coord_cartesian( ylim=c(0, 1)) +
  # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = pow, shape = methods, col = methods), stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Target FDR (%)") + ylab("Power (%)")

grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

dat <- dat4[,3]
dat <- as.data.frame(dat)
colnames(dat) <- c('pow')
dat$methods <- rep(c('MACS2', 'MACS + Clipper(median)', 
                     'MACS + Clipper(mean)', 'MACS + Clipper(max)'), each = 10)
dat$q <- rep(seq(0.01, 0.1, 0.01), 4)
pdf(file = paste0('~/Dropbox/clipper/peak calling/lineplot_peakcalling_pow2_0624.pdf'), height = 5, width = 7)
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
        legend.key = element_rect(colour = NA, fill = NA)) +
  scale_x_continuous(limits = c(0, 0.11),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10),
                     expand = c(0,0),
                     labels = function(x) paste0(x*100)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0.25, 0.5, 0.75, 1),
                     expand = c(0,0),
                     labels = c(0.25, 0.5, 0.75, 1)*100) + 
  coord_cartesian( ylim=c(0, 1)) +
  # geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5) + 
  geom_line(aes(q, y = pow, col = methods), lwd = lwd) +
  scale_color_manual(values = col_palette) +
  geom_point(aes(q, y = pow, shape = methods, col = methods), stroke = 1,size = pt_size) +
  scale_shape_manual(values = pt_shape) + xlab("Actual FDR (%)") + ylab("Power (%)")
grid.newpage()
grid.draw(set_panel_size(p, width  = unit(2, "in"),
                         height = unit(2/0.11*transformy(1), "in")))
dev.off()

