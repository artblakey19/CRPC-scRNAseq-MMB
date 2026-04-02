library(ggplot2)
library(parallel)
library(RColorBrewer)
library(vioplot)
library(arules)
compute_fdppow = function(discovery_ls, trueidx){
  sapply(discovery_ls, function(discovery){
    fdp = sum(!discovery %in% trueidx )/max(length(discovery),1)
    pow = sum(discovery %in% trueidx)/length(trueidx)
    return(c(fdp, pow))
  })
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

sample_fd = function(n, fd){
  fd_uniq = sort(unique(fd))
  n_fd_uniq = length(fd_uniq)
  int_len = fd_uniq[-1] - fd_uniq[-n_fd_uniq]
  lfend_idx = sample(size = n, x = 1:(n_fd_uniq-1), replace = T)
  unif = runif(n)
  return(fd_uniq[lfend_idx] + int_len[lfend_idx]*unif)
  # hist(fd_uniq, breaks = 100)
  # hist( fd_uniq[lfend_idx] + int_len[lfend_idx]*unif, breaks = 50)
  # plot(fd_uniq[lfend_idx],  fd_uniq[lfend_idx] + int_len[lfend_idx]*unif)
}

######## functions ######
# count1 = round(dat_i$subcount1)
# count2 = round(dat_i$subcount2)
mydeseq2 = function(count1, count2, FDR , trueDE, ifnorm = F){
  
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  ihw_cov = rowMeans(dat)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  dds <- DESeqDataSetFromMatrix(dat, DataFrame(cond_idx), ~ cond_idx)
  if(!ifnorm){
    normalizationFactors(dds) <- matrix(1, nrow = nrow(dat), ncol = ncol(dat) )
  }
  dds <- DESeq(dds)
  count_norm = counts(dds,normalized = T, replaced = F)
  
  res1 <- results(dds, alpha = 0.99999)
  pval = res1$pvalue
  names(pval) = rownames(res1)
  pvaladj = res1$padj
  names(pvaladj) = rownames(res1)
  
  fdppow_ls = sapply(FDR, function(FDR_i){
    res <- results(dds, alpha = FDR_i)
    discovery = rownames(res)[which(res$padj <= FDR_i)]
    fdp = sum(!discovery %in% trueDE)/max(1, length(discovery))
    pow = sum(discovery %in% trueDE)/length(trueDE)
    
    ihw_re = ihw(res$pvalue, covariates = ihw_cov, alpha = FDR_i)
    discovery = rownames(res)[which(adj_pvalues(ihw_re) <= FDR_i)]
    fdp_ihw = sum(!discovery %in% trueDE)/max(1, length(discovery))
    pow_ihw = sum(discovery %in% trueDE)/length(trueDE)
    c(fdp = fdp, pow = pow, fdp_ihw = fdp_ihw, pow_ihw = pow_ihw)
  })
  
  ihw_re = ihw(res1$pvalue, covariates = ihw_cov, alpha = 0.05)
  pvalues.adj = adj_pvalues(ihw_re)
  
  re = list(fdppow_ls = fdppow_ls, 
            pval = pval,
            pvaladj = pvaladj,
            pval.ihw = pvalues.adj,
            count_norm = count_norm)
  return(re)
}
# count1 = dat_i$subcount1
# count2 = dat_i$subcount2 
# FDR = FDR
# trueDE = trueDE
myedger = function(count1, count2, FDR , trueDE, ifnorm = F){
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  ### extract normalized matrix
  y <- DGEList(counts=dat,group=cond_idx)
  
  if(!ifnorm){
    y <- calcNormFactors(y, method = 'none')
  }else{
    y <- calcNormFactors(y, method = 'TMM')
  }
  
  count_norm = cpm(y)  ### when the method = 'none' in calcNormFactors, cpm converts counts to counts per million
  
  ############ recommended pipeline
  
  y <- DGEList(counts=dat,group=cond_idx)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  if(!ifnorm){
    y <- calcNormFactors(y, method = 'none')
  }else{
    y <- calcNormFactors(y, method = 'TMM')
  }
  count_norm_filtered = cpm(y)
  design <- model.matrix(~cond_idx)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  ihw_cov = rowMeans(count_norm_filtered)
  qlf1 = topTags(qlf, n = nrow(dat), p.value = 1)@.Data[[1]]
  pvalues = topTags(qlf, n = nrow(dat), p.value = 1)@.Data[[1]][,"PValue"]
  genes_left = rownames(topTags(qlf, n = nrow(dat), p.value = 1)@.Data[[1]])
  names(pvalues) = genes_left
  
  
  fdppow_ls = sapply(FDR, function(FDR_i){
    qlf_i <- topTags(qlf, n = nrow(dat), p.value = FDR_i)
    discovery = rownames(qlf_i)
    fdp = sum(!discovery %in% trueDE)/max(1, length(discovery))
    pow = sum(discovery %in% trueDE)/length(trueDE)
    
    ihw_re = ihw(pvalues, covariates = ihw_cov, alpha = FDR_i)
    discovery = genes_left[which(adj_pvalues(ihw_re) <= FDR_i)]
    fdp_ihw = sum(!discovery %in% trueDE)/max(1, length(discovery))
    pow_ihw = sum(discovery %in% trueDE)/length(trueDE)
    
    c(fdp = fdp, pow = pow, fdp_ihw = fdp_ihw, pow_ihw = pow_ihw)
  })
  ihw_re = ihw(pvalues, covariates = ihw_cov, alpha = 0.05)
  pvalues = qlf1[match(rownames(count1), rownames(qlf1)),"PValue"]
  pvalues.adj = adj_pvalues(ihw_re)[match(rownames(count1), rownames(qlf1))]
  # re = list(fdppow_ls = fdppow_ls, pval = pvalues,count_norm = count_norm, pval.ihw = pvalues.adj)
  
  re = list(fdppow_ls = fdppow_ls, 
            pval = pvalues,
            pval.ihw = pvalues.adj, ## different length 
            count_norm = count_norm, 
            count_norm_filtered = count_norm_filtered)
  return(re)
}

plotclipper = function(re_clipper, trueDE_idx){
  constr = re_clipper$contrastScore$tau * (as.numeric(re_clipper$contrastScore$kappa)*2 - 1)
  
  pdf(file = 'clipper contrast score.pdf', width = 5, height = 3)
  hist(constr[trueDE_idx], 
       breaks = 100, freq = F, col = alpha(4, alpha = 1))
  hist(constr[-trueDE_idx],
       breaks = 50, freq = F, add = T, col = alpha(2, alpha = 0.3))
  dev.off()
  return(NULL)
}


generate_subcount_spikeinDEG = function(count, r1, r2, trueDE_idx, fd_trueDE){
  
  lgfd = log(fd_trueDE)
  n_degene = length(trueDE_idx)
  idx_cond1 = as.numeric(runif(n_degene) <= 1/2)
  
  scalar1 = exp(idx_cond1*lgfd)
  scalar2 =  fd_trueDE/scalar1
  
  scaling_m1 = matrix(scalar1, ncol = 1) %*% matrix(rep(1, r1), nrow =1, byrow = T)
  scaling_m2 = matrix(scalar2, ncol = 1) %*% matrix(rep(1, r2), nrow =1, byrow = T)
  
  r1_tot = ncol(count)
  subcount = count[, sample(r1_tot, r1 + r2)]
  subcount1 = subcount[, 1:r1]
  subcount2 = subcount[, -(1:r1)]
  
  subcount1[trueDE_idx, ] = subcount1[trueDE_idx, ]*scaling_m1
  subcount2[trueDE_idx, ] = subcount2[trueDE_idx, ]*scaling_m2
  
  
  re = list(subcount1 = as.matrix(round(subcount1)), subcount2 = as.matrix(round(subcount2)))
  
  return(re)
  
}

generate_subcount_spikeNDEGbypermute = function(count1, count2, r1, r2){
  
  
  # ## generate subcounts from condition 1 by random sampling
  # subcount1 = t(apply(count1, 1, sample, size = r1, replace = T))
  # subcount2 = t(apply(count1, 1, sample, size = r2, replace = T))
  # 
  # ## add spike-in true DEG
  # subcount1[trueDE_idx,] = t(apply(count1_trueDE,1, sample, size = r1, replace = T))
  # subcount2[trueDE_idx,] = t(apply(count2_trueDE,1, sample, size = r2, replace = T))
  # re = list(subcount1 = subcount1, subcount2 = subcount2)
  
  subcount1 = count1[,sample(ncol(count1), r1, replace = F)]
  subcount2 = count2[,sample(ncol(count2), r2, replace = F)]
  
  subcount_tot_nDE = cbind(subcount1[-trueDE_idx,], subcount2[-trueDE_idx,])
  subcount_tot_nDE = t(apply(subcount_tot_nDE, 1, sample, size = r1 + r2, replace = F))
  subcount1[-trueDE_idx,] = subcount_tot_nDE[,1:r1]
  subcount2[-trueDE_idx,] = subcount_tot_nDE[,-(1:r1)]

  
  
  # lgm1 = log(base =2, rowMeans(subcount1) + 1 )
  # lgm2 = log(base =2, rowMeans(subcount2) + 1 )
  # plot(lgm1, lgm2, cex = 0.3, col = 1)
  # points(lgm1[trueDE_idx], lgm2[trueDE_idx], col = 2, cex = 0.3)
  # points(lgm1[-trueDE_idx], lgm2[-trueDE_idx], col = 3, cex = 0.3)
  re = list(subcount1 = as.matrix(round(subcount1)), subcount2 = as.matrix(round(subcount2)))
  return(re)
}

Compare = function(spikeInTarget){
  mclapply(1:iter, function(it){
    set.seed(it)
    if(spikeInTarget == 'DEG'){
      dat_i = generate_subcount_spikeinDEG(count = count2,r1, r2, trueDE_idx = trueDE_idx, fd_trueDE = fd_trueDE)
    }
    if(spikeInTarget == 'NDEG'){
      dat_i = generate_subcount_spikeNDEGbypermute(count1 = count1, count2 = count2, r1, r2)
      # saveRDS(dat_i, file  = 'new.rds')
      
    }
    # subcount1 = dat_i$subcount1
    # subcount2 = dat_i$subcount2
    # plot(rowMeans(subcount1), rowMeans(subcount2), cex = 0.3)
    # points(rowMeans(subcount1)[trueDE_idx],rowMeans(subcount2)[trueDE_idx], cex = 0.3, col =2)
    # abline(a = 0, b = 1)
    ##### deseq2 unnorm
    re_deseq2 = suppressMessages(mydeseq2(count1 = round(dat_i$subcount1), count2 = round(dat_i$subcount2), FDR = FDR, trueDE = trueDE, ifnorm = F))
    fdppow_deseq2_unnorm = re_deseq2$fdppow_ls
    fdppow_deseq2_unnorm_ihw = fdppow_deseq2_unnorm[ 3:4, ]
    fdppow_deseq2_unnorm = fdppow_deseq2_unnorm[ 1:2, ]
    pval_deseq2_unnorm = re_deseq2$pval
    fdppow.i_deseq2_unnorm <- sapply(1:10, function(j){
      compute_fdppow2(pval_deseq2_unnorm, trueDE_idx, FDR[j])
    })
    pval_deseq2_unnorm_ihw = re_deseq2$pval.ihw
    fdppow.i_deseq2_unnorm_ihw <- sapply(1:10, function(j){
      compute_fdppow2(pval_deseq2_unnorm_ihw, trueDE_idx, FDR[j])
    })
    ##### deseq2 norm
    re_deseq2 = suppressMessages(mydeseq2(count1 = round(dat_i$subcount1), count2 = round(dat_i$subcount2), FDR = FDR, trueDE = trueDE, ifnorm = T))
    fdppow_deseq2_norm = re_deseq2$fdppow_ls
    fdppow_deseq2_norm_ihw = fdppow_deseq2_norm[ 3:4, ]
    fdppow_deseq2_norm = fdppow_deseq2_norm[ 1:2, ]
    pval_deseq2_norm = re_deseq2$pval
    # count_norm_deseq2 = re_deseq2$count_norm
    fdppow.i_deseq2_norm <- sapply(1:10, function(j){
      compute_fdppow2(pval_deseq2_norm, trueDE_idx, FDR[j])
    })
    pval_deseq2_norm_ihw = re_deseq2$pval.ihw
    fdppow.i_deseq2_norm_ihw <- sapply(1:10, function(j){
      compute_fdppow2(pval_deseq2_norm_ihw, trueDE_idx, FDR[j])
    })
    ##### edgeR unnorm
    re_edger = myedger(count1 = dat_i$subcount1, count2 = dat_i$subcount2 , FDR = FDR, trueDE = trueDE, ifnorm = F)
    fdppow_edger_unnorm = re_edger$fdppow_ls
    fdppow_edger_unnorm_ihw = fdppow_edger_unnorm[3:4, ]
    fdppow_edger_unnorm = fdppow_edger_unnorm[1:2,]
    pval_edger_unnorm = re_edger$pval
    fdppow.i_edger_unnorm <- sapply(1:10, function(j){
      compute_fdppow2(pval_edger_unnorm, trueDE_idx, FDR[j])
    })
    pval_edger_unnorm_ihw = re_edger$pval.ihw
    fdppow.i_edger_unnorm_ihw <- sapply(1:10, function(j){
      compute_fdppow2(pval_edger_unnorm_ihw, trueDE_idx, FDR[j])
    })
    
    ##### edgeR norm
    re_edger = myedger(count1 = dat_i$subcount1, count2 = dat_i$subcount2 , FDR = FDR, trueDE = trueDE, ifnorm = T)
    fdppow_edger_norm = re_edger$fdppow_ls
    fdppow_edger_norm_ihw = fdppow_edger_norm[3:4, ]
    fdppow_edger_norm = fdppow_edger_norm[1:2,]
    pval_edger_norm = re_edger$pval
    count_norm_edger = re_edger$count_norm
    fdppow.i_edger_norm <- sapply(1:10, function(j){
      compute_fdppow2(pval_edger_norm, trueDE_idx, FDR[j])
    })
    pval_edger_norm_ihw = re_edger$pval.ihw
    fdppow.i_edger_norm_ihw <- sapply(1:10, function(j){
      compute_fdppow2(pval_edger_norm_ihw, trueDE_idx, FDR[j])
    })
    # ##### clipper unnorm diff max 1
    # re_clipper = clipper2sided(score_exp = log(base = 2, dat_i$subcount1 + 1), 
    #                            score_back = log(base = 2, dat_i$subcount2 + 1), 
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'max',
    #                            nknockoff = 1,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames(dat_i$subcount1)[x$discovery]
    # })
    # fdppow_clipper_unnorm_diff_max_1 = compute_fdppow(discovery, trueidx = trueDE)
    
    ##### clipper unnorm diff max 9
    re_clipper = clipper2sided(score_exp = log(base = 2, dat_i$subcount1 + 1), 
                               score_back = log(base = 2, dat_i$subcount2 + 1), 
                               FDR = FDR,
                               importanceScore_method = 'diff',
                               contrastScore_method = 'max',
                               nknockoff = 9,
                               FDR_control_method = 'GZ',
                               ifpowerful = F)
    discovery = lapply(re_clipper$results,function(x){
      rownames(dat_i$subcount1)[x$discovery]
    })
    fdppow_clipper_unnorm_diff_max_9 = compute_fdppow(discovery, trueidx = trueDE)
    contrast_clipper_unnorm_diff_max_9 <- -(2*re_clipper$contrastScore$kappa - 1) * re_clipper$contrastScore$tau
    fdppow.i_clipper_unnorm_diff_max_9  <- sapply(1:10, function(j){
      compute_fdppow2(contrast_clipper_unnorm_diff_max_9, trueDE_idx, FDR[j])
    })
    
    # ##### clipper unnorm diff diff 1
    # re_clipper = clipper2sided(score_exp = log(base = 2, dat_i$subcount1 + 1), 
    #                            score_back = log(base = 2, dat_i$subcount2 + 1), 
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'diff',
    #                            nknockoff = 1,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames(dat_i$subcount1)[x$discovery]
    # })
    # fdppow_clipper_unnorm_diff_diff_1 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    # ##### clipper unnorm diff diff 9
    # re_clipper = clipper2sided(score_exp = log(base = 2, dat_i$subcount1 + 1), 
    #                            score_back = log(base = 2, dat_i$subcount2 + 1), 
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'diff',
    #                            nknockoff = 9,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames(dat_i$subcount1)[x$discovery]
    # })
    # fdppow_clipper_unnorm_diff_diff_9 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    # ##### clipper deseq2norm diff max 1
    # re_clipper = clipper2sided(score_exp = log(base = 2, count_norm_deseq2[, 1:3] + 1), 
    #                            score_back = log(base = 2, count_norm_deseq2[, 4:6] + 1), 
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'max',
    #                            nknockoff = 1,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames( count_norm_deseq2)[x$discovery]
    # })
    # fdppow_clipper_deseq2norm_diff_max_1 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    # ##### clipper deseq2norm diff max 9
    # re_clipper = clipper2sided(score_exp = log(base = 2, count_norm_deseq2[, 1:3] + 1),
    #                            score_back = log(base = 2, count_norm_deseq2[, 4:6] + 1),
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'max',
    #                            nknockoff = 9,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames( count_norm_deseq2)[x$discovery]
    # })
    # fdppow_clipper_deseq2norm_diff_max_9 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    
    # ##### clipper edgernorm diff max 1
    # re_clipper = clipper2sided(score_exp = log(base = 2, count_norm_edger[, 1:3] + 1), 
    #                            score_back = log(base = 2, count_norm_edger[, 4:6] + 1), 
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'max',
    #                            nknockoff = 1,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames( count_norm_edger)[x$discovery]
    # })
    # fdppow_clipper_edgernorm_diff_max_1 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    
    ##### clipper edgernorm diff max 9
    re_clipper = clipper2sided(score_exp = log(base = 2, count_norm_edger[, 1:3] + 1),
                               score_back = log(base = 2, count_norm_edger[, 4:6] + 1),
                               FDR = FDR,
                               importanceScore_method = 'diff',
                               contrastScore_method = 'max',
                               nknockoff = 9,
                               FDR_control_method = 'GZ',
                               ifpowerful = F)
    discovery = lapply(re_clipper$results,function(x){
      rownames( count_norm_edger)[x$discovery]
    })
    fdppow_clipper_edgernorm_diff_max_9 = compute_fdppow(discovery, trueidx = trueDE)
    contrast_clipper_edgernorm_diff_max_9 <- -(2*re_clipper$contrastScore$kappa - 1) * re_clipper$contrastScore$tau
    fdppow.i_clipper_edgernorm_diff_max_9  <- sapply(1:10, function(j){
      compute_fdppow2(contrast_clipper_edgernorm_diff_max_9, trueDE_idx, FDR[j])
    })
    
    # ##### clipper edgernorm diff diff 1
    # re_clipper = clipper2sided(score_exp = log(base = 2, count_norm_edger[, 1:3] + 1), 
    #                            score_back = log(base = 2, count_norm_edger[, 4:6] + 1), 
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'diff',
    #                            nknockoff = 1,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames( count_norm_edger)[x$discovery]
    # })
    # fdppow_clipper_edgernorm_diff_diff_1 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    # 
    # ##### clipper edgernorm diff diff 9
    # re_clipper = clipper2sided(score_exp = log(base = 2, count_norm_edger[, 1:3] + 1),
    #                            score_back = log(base = 2, count_norm_edger[, 4:6] + 1),
    #                            FDR = FDR,
    #                            importanceScore_method = 'diff',
    #                            contrastScore_method = 'diff',
    #                            nknockoff = 9,
    #                            FDR_control_method = 'GZ',
    #                            ifpowerful = F)
    # discovery = lapply(re_clipper$results,function(x){
    #   rownames( count_norm_edger)[x$discovery]
    # })
    # fdppow_clipper_edgernorm_diff_diff_9 = compute_fdppow(discovery, trueidx = trueDE)
    # 
    # 
    

    
    fdp = list( deseq2_unnorm = fdppow_deseq2_unnorm[1,],
                deseq2_norm = fdppow_deseq2_norm[1,],
                deseq2_unnorm_ihw = fdppow_deseq2_unnorm_ihw[1, ],
                deseq2_norm_ihw = fdppow_deseq2_norm_ihw[1, ],
                
                edger_unnorm =  fdppow_edger_unnorm[1,], 
                edger_norm =  fdppow_edger_norm[1,], 
                edger_unnorm_ihw = fdppow_edger_unnorm_ihw[1, ],
                edger_norm_ihw = fdppow_edger_norm_ihw[1, ],
                
                clipper_unnorm_diff_max_9 = fdppow_clipper_unnorm_diff_max_9[1, ],
                # clipper_unnorm_diff_max_1 = fdppow_clipper_unnorm_diff_max_1[1, ],
                # clipper_unnorm_diff_diff_9 = fdppow_clipper_unnorm_diff_diff_9[1, ],
                # clipper_unnorm_diff_diff_1 = fdppow_clipper_unnorm_diff_diff_1[1, ],
                # clipper_deseq2norm_diff_max_9 = fdppow_clipper_deseq2norm_diff_max_9[1, ],
                # clipper_deseq2norm_diff_max_1 = fdppow_clipper_deseq2norm_diff_max_1[1, ],
                clipper_edgernorm_diff_max_9 = fdppow_clipper_edgernorm_diff_max_9[1, ]
                # clipper_edgernorm_diff_max_1 = fdppow_clipper_edgernorm_diff_max_1[1, ],
                # clipper_edgernorm_diff_diff_9 = fdppow_clipper_edgernorm_diff_diff_9[1,],
                # clipper_edgernorm_diff_diff_1 = fdppow_clipper_edgernorm_diff_diff_1[1,]
                
    )
    pow = list( deseq2_unnorm = fdppow_deseq2_unnorm[2,],
                deseq2_norm = fdppow_deseq2_norm[2,],
                deseq2_unnorm_ihw = fdppow_deseq2_unnorm_ihw[2, ],
                deseq2_norm_ihw = fdppow_deseq2_norm_ihw[2, ],
                
                edger_unnorm =  fdppow_edger_unnorm[2,], 
                edger_norm =  fdppow_edger_norm[2,], 
                edger_unnorm_ihw = fdppow_edger_unnorm_ihw[2, ],
                edger_norm_ihw = fdppow_edger_norm_ihw[2, ],
                clipper_unnorm_diff_max_9 = fdppow_clipper_unnorm_diff_max_9[2, ],
                # clipper_unnorm_diff_max_1 = fdppow_clipper_unnorm_diff_max_1[2, ],
                # clipper_unnorm_diff_diff_9 = fdppow_clipper_unnorm_diff_diff_9[2, ],
                # clipper_unnorm_diff_diff_1 = fdppow_clipper_unnorm_diff_diff_1[2, ],
                # clipper_deseq2norm_diff_max_9 = fdppow_clipper_deseq2norm_diff_max_9[2, ],
                # clipper_deseq2norm_diff_max_1 = fdppow_clipper_deseq2norm_diff_max_1[2, ],
                clipper_edgernorm_diff_max_9 = fdppow_clipper_edgernorm_diff_max_9[2, ]
                # clipper_edgernorm_diff_max_1 = fdppow_clipper_edgernorm_diff_max_1[2, ],
                # clipper_edgernorm_diff_diff_9 = fdppow_clipper_edgernorm_diff_diff_9[2,],
                # clipper_edgernorm_diff_diff_1 = fdppow_clipper_edgernorm_diff_diff_1[2,]
    )
    pow.i = list( deseq2_unnorm = fdppow.i_deseq2_unnorm[2,],
                  deseq2_norm = fdppow.i_deseq2_norm[2,],
                  deseq2_unnorm_ihw =  fdppow.i_deseq2_unnorm_ihw[2, ],
                  deseq2_norm_ihw =  fdppow.i_deseq2_norm_ihw[2, ],
                  # clipper_deseq2 = fdppow_clipper_deseq2[2,],
                  edger_unnorm =  fdppow.i_edger_unnorm[2,], 
                  edger_norm =  fdppow.i_edger_norm[2,], 
                  edger_unnorm_ihw = fdppow.i_edger_unnorm_ihw[2, ],
                  edger_norm_ihw = fdppow.i_edger_norm_ihw[2, ],
                  
                  # clipper_edger = fdppow_clipper_edger[2,],
                  clipper_unnorm_diff_max_9  = fdppow.i_clipper_unnorm_diff_max_9[2,],
                  clipper_edgernorm_diff_max_9  = fdppow.i_clipper_edgernorm_diff_max_9[2,])
    pval = list(deseq2_unnorm = pval_deseq2_unnorm,
                deseq2_norm = pval_deseq2_norm,
                edger_unnorm = pval_edger_unnorm,
                edger_norm = pval_edger_norm,
                clipper_unnorm_diff_max_9 = contrast_clipper_unnorm_diff_max_9,
                clipper_edgernorm_diff_max_9 = contrast_clipper_edgernorm_diff_max_9)
    # pval = list(deseq2 = pval_deseq2,
    #             edger = pval_edger)
    return(list(fdp = fdp, pow = pow, pow.i = pow.i, pval = pval))
  },mc.cores = ncores)
  
}


coarseplot = function(re){
  
  # dat = re
  # print(length(re[[1]]$fdp))
  fdp = lapply(re, '[[','fdp')
  fdp = lapply(1:length(re[[1]]$fdp), function(i){
    rowMeans(sapply(fdp, '[[',i))
  })
  
  pow = lapply(re, '[[','pow')
  pow = lapply(1:length(re[[1]]$fdp), function(i){
    rowMeans(sapply(pow, '[[',i))
  })
  
  fdp
  pow
  method_names = names(re[[1]]$fdp)
  dat = cbind.data.frame(fdr = unlist(fdp),
                         pow = unlist(pow),
                         methods = rep(method_names, each = 10),
                         q = rep(FDR, length(method_names)))
  p1 = ggplot(dat,aes(x = q, y = fdr, group = methods)) +
    geom_line(aes(q, y = fdr, col = methods), lwd = 0.5) + 
    geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5)  +
    scale_color_manual(values = brewer.pal(length(method_names),'Paired'))
  
  p1 = set_panel_size(p1, width  = unit(2, "in"),
                      height = unit(4, "in"))
  p2 = ggplot(dat,aes(x = q, y = pow, group = methods)) +
    geom_line(aes(q, y = pow, col = methods), lwd = 0.5) +
    scale_color_manual(values = brewer.pal(length(method_names),'Paired'))
  p2 = set_panel_size(p2, width  = unit(2, "in"),
                      height = unit(4, "in"))
  return(list(p1 = p1, p2 = p2))
}

test_uniformpval = function(re, trueDE_idx){
  method_ls = c('deseq2_norm','edger_norm')
  pval_re = lapply(method_ls, function(m){
    pval = lapply(re, '[[','pval')
    pval = lapply(pval, '[[', m)
    pval = as.matrix(as.data.frame(pval))
    pval = pval[-trueDE_idx, ]
    # i = 1
    
    pval_ifuniform = apply(pval, 1, function(x){
      # print(i)
      # x = unique(x)
      # if(length(unique(x)) < length(x)){
      #   return(NA)
      # }else{
      #   ks.test(x,"punif",0,1)$p.value
      #   
      # }
      na_prop = mean(is.na(x))
      x = x[!is.na(x)]
      if(length(x) ==0 ){
        return(c(1, NA))
      }else{
        x = discretize(c(x), method = 'fixed', breaks = (0:10)/10)
        re_test = chisq.test(table(x), p = rep(0.1, 10))
        return(c(na_prop, re_test$p.value))
      }
      # i = i+1
    
    })
  })
  
  names(pval_re) = method_ls
  return(pval_re)
  
  
}
