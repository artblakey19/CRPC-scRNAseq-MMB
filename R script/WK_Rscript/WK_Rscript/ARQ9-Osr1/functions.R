library(HiTC)
library(multiHiCcompare)
library(HiCcompare)
library(diffHic)
library(edgeR)
library(FIND)
library(csaw)
library(mixtools)
convertfull2bin = function(x){
  ij = which(x != 0, arr.ind = T)
  values = x[x!=0]
  values = values[ij[,1] <= ij[,2]]
  ij = ij[ij[,1] <= ij[,2], ]
  return(cbind(ij, values))

}
# temp = convertfull2bin(convert2matrix(rep1_binid))
# all(temp == rep1_binid)
# dat_binid = rep1_binid
# my_icenorm = function(dat_binid, resolution = 10^6){
#   m = convert2matrix(dat_binid)
#   n = nrow(m)
#   colnames(m ) = as.character(1:nrow(m))
#   rownames(m) = as.character(1:nrow(m))
#   rsm = rowSums(m)
#   m = m[rsm > 0 ,rsm > 0 ]
#   allregions = GRanges('chr1',
#                        ranges = IRanges(((0:(n-1))*resolution + 1)[rsm > 0 ], ((1:n)*resolution)[rsm > 0 ],
#                        names = as.character(1:n)[rsm > 0 ]))
#   # names(ranges(allregions)) =  as.character(1:nrow(m))
#   
#   m_HTCexp = new('HTCexp',as(m,'dgeMatrix'), allregions, allregions)
#   m_norm = normICE(m_HTCexp)
#   m_norm = m_norm@intdata
#   m_norm = as.matrix(m_norm)
#   return(m_norm)
# }
# 
# temp = my_icenorm(rep1_binid)

averagescores = function(dat_ls){
  feature_col = dat_ls[[1]][,1:2]
  scores = rowMeans(sapply(dat_ls, function(x){x[,3]}))
  return(cbind(feature_col, scores))
}

my_multihiccompare = function(score_exp_ls, score_back_ls, true_binidx_char, FDR, resolution = 10^6, ifnormalize = T){
  n_exp = length(score_exp_ls)
  n_back = length(score_back_ls)
  score_exp_ls = lapply(score_exp_ls, function(x){
    x = convert2bp(x)
    x = cbind(1, x)
    colnames(x) = c('chr', 'region1', 'region2', 'IF')
    return(as.data.frame(x))
  })
  
  score_back_ls = lapply(score_back_ls, function(x){
    x = convert2bp(x)
    x = cbind(1, x)
    colnames(x) = c('chr', 'region1', 'region2', 'IF')
    return(as.data.frame(x))
  })
  score_ls = c(score_exp_ls, score_back_ls)
  hicexp1 = make_hicexp(score_ls[[1]], score_ls[[2]], score_ls[[3]], score_ls[[4]] , 
                        groups = c(rep(0, n_exp), rep(1, n_back)),
                        filter = F)
  # data("hicexp2")
  if(ifnormalize){
    hicexp1 <- cyclic_loess(hicexp1, verbose = FALSE,
                            parallel = FALSE, span = 0.2)
  }
  
  norm_count = hic_table(hicexp1)[,-(1:4)]
  
  # MD_hicexp(hicexp1)
  hicexp1 <- hic_exactTest(hicexp1, p.method = 'fdr',
                           parallel = FALSE)
  region1 = results(hicexp1)$region1/resolution + 1
  region2 = results(hicexp1)$region2/resolution + 1
  features = convert2char(region1, region2)
  discovery_ls = sapply(FDR, function(FDR_i){
    
    features[results(hicexp1)$p.adj <= FDR_i]
  })
  
  fdppow = compute_fdppow(discovery_ls, trueidx = true_binidx_char)
  re = list(fdppow = fdppow, norm_count = norm_count)
  
  temp = cbind(region1, region2, -log10(results(hicexp1)$p.adj))
  print(levelplot(convert2matrix(temp)))
  # temp2 = cbind(mean_back[,1:2], abs(log(mean_exp[,3]/mean_back[,3])))
  # levelplot(convert2matrix(temp2))
  
  return(re)
}

# 
# score_exp_ls = dat_ls_binned$exp
# score_back_ls = dat_ls_binned$back
# design = c(rep(1,n_exp), rep(2, n_back))
# features = totrow_idx
# FDR = FDR
# iffilter = F
# truerow_idx = truerow_idx
my_diffhic = function( score_exp_ls, score_back_ls,  design, FDR, truerow_idx, iffilter = T, ifnormalize = T){
  # features = 1:nrow(score_exp_ls[[1]])
  cm_exp_ls = lapply(score_exp_ls,function(x){
    my_ContactMatrix(x)
  })
  
  cm_back_ls = lapply(score_back_ls,function(x){
    my_ContactMatrix(x)
  })
  cm_ls = c(cm_exp_ls, cm_back_ls)
  # idx_extract = lapply(cm_ls, function(x){
  #   as.matrix(x )!=0
  # })
  # idx_extract = Reduce( '|', idx_extract)
  # # levelplot(as.matrix(cm_ls[[4]]))
  # iset_ls = lapply(cm_ls, function(x){
  #   deflate(x, extract = idx_extract)
  # })
  data <- do.call(mergeCMs, cm_ls)
  anchor1 = data@interactions@anchor1
  anchor2 = data@interactions@anchor2
  features = convert2char(anchor2, anchor1) ### anchor2 <= anchor1
  
 #  interactions(data) <- as(interactions(data), "ReverseStrictGInteractions")
 # temp = data@assays@data@listData[[1]]
 #  temp = cbind(score_exp_ls[[1]][,1:2], temp[,1])
 #  levelplot(convert2matrix(temp))
 #  
 #  anchor1 == rep1_binid[,1]
 #  temp = cbind(anchor1, anchor2, temp [,1])
  # iset1 <- deflate(cm1, extract=to.keep)
  # iset2 <- deflate(cm2, extract=to.keep)
  # data <- cbind(iset1, iset2)
  # 
  # data  = do.call( mergeCMs, cm_ls)
  

  
  if(iffilter){
    keep <- aveLogCPM(asDGEList(data)) > 0
    data <- data[keep,]
    features <- features[keep]
  }
  count_norm = data@assays@data@listData$counts
  if(ifnormalize){
    data <- normOffsets(data, method ="loess", se.out = T)
    logoffset = data@assays@data@listData$offset
    count_norm = log2(count_norm + 0.5)- assay(data, "offset")/log(2)
    count_norm = 2^count_norm
    # counts = data@assays@data@listData$counts
    # count_norm = counts/exp(logoffset)
    # data@assays@data@listData$counts = count_norm
    # data@assays@data@listData$offset = matrix(0, nrow = nrow(count_norm), ncol = ncol(count_norm))
  }
  
  y <- asDGEList(data)
  if(!is.factor(design)){
    design = factor(design)
  }
  # colnames(y@.Data[[1]]) = paste0('V',1:4)
  design <- model.matrix(~design)
  y <- estimateDisp(y, design)
  # plotBCV(y)
  # fit = glmQLFit(y, design)
  # qlf = glmLRT(fit)
  fit <- glmQLFit(y, design, robust=TRUE)
  # plotQLDisp(fit)
  qlf <- glmQLFTest(fit)
  
  # all(as.numeric(rownames(qlf@.Data[[1]])) %in% 1:length(rownames(qlf@.Data[[1]])))
  

  fdppow_ls = sapply(1:length(FDR), function(i_FDR){
    
    qlf_i <- topTags(qlf, n = 10^10, p.value = FDR[i_FDR], sort.by = 'none')
   
    discovery = features[as.numeric(rownames(qlf_i))]
    # discovery_ls[[i_FDR]] = discovery
    fdp = sum(!discovery %in% true_binidx_char)/max(1, length(discovery))
    pow = sum(discovery %in% true_binidx_char)/length(true_binidx_char)
    c(fdp = fdp, pow =pow)
  })
  fdppow_ls
  
  # qlf_i <- topTags(qlf, n = 10^10, p.value = 1, sort.by = 'none')
  # temp = cbind(anchor1, anchor2, qlf_i$table$FDR)
  # levelplot((convert2matrix(temp)))
  return( list(fdppow = fdppow_ls,features = features, count_norm = count_norm))
}


# adjustforFDR_diffhic = function(qlf, qval, FDR, features,truerow_idx){
#   score= sort(unique(as.vector(qval)))
#   score = score - min(score)
#   score = score/max(score)
#   score_sorted = sort(unique(score),decreasing = T)
#   n = length(score_sorted)
#   i = 1
#   fdp = 2
#   fdp_ls = rep(NA, n )
#   pow_ls = fdp_ls
#   while( i < n & fdp > min(FDR) ){
#     thre = score_sorted[i]
#     qlf_i <- topTags(qlf, n = 10^10, p.value = FDR_i)
#     # View(qlf_i$table)
#     discovery = features[as.numeric(rownames(qlf_i))]
#     fdp = sum(!as.numeric(discovery) %in% truerow_idx)/max(1, length(discovery))
#     pow = sum(discovery %in% truerow_idx)/length(truerow_idx)
#     fdp_ls[i] = fdp
#     pow_ls[i] = pow
#     i = i + 1
#   }
# }
# temp = cbind(rep1_binid[features,1:2], qlf_i@.Data[[1]][,'FDR'])
# temp2 = convert2matrix(temp)
# levelplot(-log10(temp2))
# score_exp_ls = dat_ls_binned$exp
# score_back_ls = dat_ls_binned$back
# FDR = FDR
# true_binidx_char = true_binidx_char
################ read FIND paper ################ 
################ read hiccompare ################ 
################ read hicdiff ################ 
################ add normalization 
my_find = function(score_exp_ls, score_back_ls, FDR, true_binidx_char){
  
  hic1 = lapply(score_exp_ls, function(x){convert2matrix(x)})
  hic2 = lapply(score_back_ls, function(x){convert2matrix(x)})
  hic1 = lapply(hic1, as, 'dgCMatrix')
  hic2 = lapply(hic2, as, 'dgCMatrix')
  dci = getDCIs_fromMat(Hic1 = hic1, Hic2 = hic2,  windowSize = 3, 
                        alpha = 0.7, method = "hardCutof", qvalue = 1, 
                        isrOP_qval = FALSE, resolution = 10^6)
  qval = as(dci@qvals,'matrix') 
  discovery_ls = lapply(FDR, function(FDR_i){
    ind = which(qval <= FDR_i, arr.ind = T)
    ind = ind[ ind[,1] <= ind[,2],]  ## extract upper triangular indices
    ind = convert2char(ind[,1],ind[,2])
    ind
  }) ################ 
  
  re = compute_fdppow(discovery_ls, trueidx = true_binidx_char)
  return(re)
}

# dci = readRDS(file = 'dci.rds')
# qval = as(dci@qvals,'matrix') 
# levelplot(-log10(qval))
# dat_binid = rep1_binid
my_ContactMatrix = function(dat_binid, resolution = 10^6){
  # dat_binid = convert2binned(dat)
  m = convert2matrix(dat_binid)
  allregions = GRanges('chr1',
                       IRanges((0:(nrow(m)-1))*resolution, (1:nrow(m))*resolution))
  cm = ContactMatrix( m, anchor1 = allregions, anchor2 = allregions)
  # temp = as.matrix(cm)
  # levelplot(temp)
  return(cm)
}

convert2matrix = function(dat){
  n = max(dat[,1:2])
  m = matrix(0, nrow = n, ncol = n )
  for(i in 1:nrow(dat)){
    m[dat[i,1],dat[i,2]] = dat[i,3]
    m[dat[i,2],dat[i,1]] = dat[i,3]
  }
  m[is.na(m)] = 0 
  return(m)
}

convert2char = function(start1, start2){
  n = length(start1)
  char_idx = sapply(1:n, function(i){
    paste(as.character(c(start1[i], start2[i])), collapse = ':')
  })
  return(char_idx)
}

convert2binned = function(dat,resolution = 10^6){
  as.matrix(cbind(dat[,1:2]/resolution + 1,dat[,3]))
} 

convert2bp = function(dat,resolution = 10^6){
  as.matrix(cbind((dat[,1:2]-1)*resolution,dat[,3]))
}

generate_data_frommean = function(mean, totrow_idx){
  mean = as.matrix(mean)
  simudat = sapply(mean[,3], function(x){
    # max(rpois(n = 1, lambda = x),1)
    max(rnbinom(n = 1, mu = x, size = 1000),1)
  })
  re = cbind(mean[,1:2],simudat)
  rownames(re) = totrow_idx
  return(re)
}

compute_spikein_mean = function(dat, true_idx, lgfc){
  dat = as.matrix(dat)
  n = nrow(true_idx)
  true_idx_char = sapply(1:n, function(i){
    paste(as.character(unlist(true_idx[i,])), collapse = ':')
  })
  dat_idx_char = sapply(1:nrow(dat), function(i){
    paste(as.character(unlist(dat[i,1:2])), collapse = ':')
  })
  mean_trueidx = dat[match(true_idx_char, dat_idx_char),3]*exp(lgfc)
  dat_new = dat
  dat_new[match(true_idx_char, dat_idx_char),3] = mean_trueidx
  dat_new = round(dat_new)
  return(dat_new)
}
# mean_exp = compute_spikein_mean(rep1, true_idx, lgfc)
# m_exp = convert2matrix(mean_exp)
# levelplot(log(m_exp))



# mean_exp = compute_spikein_mean(rep1, true_idx, lgfc)
# mean_back = rep2
# n_exp = 1
# n_back = 2
generatedata = function(totrow_idx, mean_exp, mean_back,n_exp, n_back){
  score_exp = lapply(1:n_exp, function(i){
    generate_data_frommean(mean = mean_exp, totrow_idx = totrow_idx)
  })
  score_back = lapply(1:n_back, function(i){
    generate_data_frommean(mean = mean_back, totrow_idx = totrow_idx)
  })  
  return(list(exp = score_exp, back = score_back))
}


# dat_ls = dat_ls_binned
format2clipperinput = function(dat_ls, n_exp, n_back){
  dat_ls = unlist(dat_ls, recursive = F)
  n_dat = length(dat_ls)
  idx_char_ls = lapply(dat_ls, function(x){
    convert2char(start1 = x[,1],start2 = x[,2])
  })
  idx_char_tot = Reduce('union',idx_char_ls)
  count_ls = lapply(1:n_dat, function(i_dat){
    count_new = dat_ls[[i_dat]][match(idx_char_tot, idx_char_ls[[i_dat]]),3]
    count_new[is.na(count_new)] = 0
    count_new
  })
  re =  list(feature = idx_char_tot,
             exp = Reduce('cbind',count_ls[1:n_exp]), back = Reduce('cbind',count_ls[-(1:n_exp)]))
  return(re)
}


compute_fdppow = function(discovery_ls, trueidx){
  sapply(discovery_ls, function(discovery){
    fdp = sum(!discovery %in% trueidx )/max(length(discovery),1)
    pow = sum(discovery %in% trueidx)/length(trueidx)
    return(c(fdp, pow))
  })
}


### can only deal with two sparse matrix
my_hiccompare = function(score_exp, score_back, FDR, true_binidx_char){
  chr1.table <- create.hic.table(score_exp, score_back, chr = 'chr1')
  # head(chr1.table)
  hic.table <- hic_loess(chr1.table, Plot = F, Plot.smooth = FALSE) ## contains the normalized counts
  hic.table <- hic_compare(hic.table, A.min = 15, adjust.dist = TRUE, p.method = 'fdr', Plot = F)
  start1 = hic.table$start1/resolution + 1
  start2 = hic.table$start2/resolution + 1
  feature = convert2char(start1, start2)
  discovery_ls = sapply(FDR, function(FDR_i){
    feature[hic.table$p.adj <= FDR_i]
  })
  fdppow = compute_fdppow(discovery_ls, trueidx = true_binidx_char)
  norm_count = list(exp = hic.table$adj.IF1, back = hic.table$adj.IF2)
  re = list(fdppow = fdppow, norm_count = norm_count)
  return(re)
}

normalize = function(score_exp, score_back, pos_thre = 0.8){
  logscore_exp = log10(rowMeans(score_exp))
  logscore_back = log10(rowMeans(score_back))
  
  logscore_tot = cbind(exp = logscore_exp, back = logscore_back)
  emfit_ls = apply(logscore_tot,2, function(x){
    hist(x, breaks = 100)
    emfit = normalmixEM((x))
    pos = emfit$posterior[, which.min(emfit$mu)]
    idx = pos > pos_thre
    which(idx)
  })
  idx_inter = Reduce(intersect, emfit_ls)
  
  dat = as.data.frame(logscore_tot[idx_inter,])
  # names(dat) = c('y','x')
  # plot(dat$x, dat$y, cex = 0.3, main = 'sample', xlim = c(1, 4), ylim =c(1,4))
  # points(x = log10(mean_back[-truerow_idx,3]), y = log10(mean_exp[-truerow_idx,3]), cex = 0.3, col = 2)
  # abline(a = 0, b = 1, col = 3)
  
  fit1 = lm(exp~back, data = dat)
  # plot(fit1)
  # summary(fit1)
  # plot(fit1$fitted.values,residuals(fit1))
  # abline(h = 0.3)
  # abline(h = -0.3)
  res = residuals(fit1)
  res_se = sqrt(sum(res^2)/length(fit1$fitted.values))
  idx2keep = (abs(residuals(fit1)/res_se) <= 2 )
  # plot(y = res/res_se, x =fit1$fitted.values)
  # abline(h = -2)
  # plot(dat[idx2keep,])
  fit2 = lm(exp~back, data = dat[idx2keep,])
  coefs = coefficients(fit2)
  
  # logscore_exp_norm = (logscore_exp - coefs[1])/(coefs[2])
  # plot(back = logscore_back[idx_inter][idx2keep], 
  #      exp = logscore_exp_norm[idx_inter][idx2keep],
  #      cex = 0.3)
  # abline( a = 0, b = 1, col = 2 )
  # lm()
  score_exp_norm = 10^(log10(score_exp) - coefs[1])/(coefs[2])
  re = list(logcoefs = coefs, score_exp_norm = score_exp_norm, score_back_norm = score_back)
  return(re)
}


