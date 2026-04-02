remove(list = ls())
#setwd('~/Dropbox/clipper/source042520/2sided')
#source('~/Dropbox/clipper/source042520/clipper042520_w2sided.R', echo=F)
# library(trunc)
#source('clipper042520_w2sided.R', echo=F)
source('clipper102520.R', echo=F)
library(parallel)
library(locfdr)
library(truncdist)
library(ks)
library(MASS)
library(qvalue)

ncores = detectCores() - 2
iter = 200

n <- 10000
q_tot = (1:10)/100
out_prop <- c(0, 0.1)


####### functions ######

compute_pooled_pval_2s = function(back, exp){
  back_median <- median(back, na.rm = T)
  re <- sapply(exp, function(x){
    if (is.na(x)){return(NA)}
    if (x >=back_median) {
      return(min(1,2*mean(back >= x, na.rm = T)))
    }else{
      return(min(1,2*mean(back <= x, na.rm = T)))
    }})
  return(re)
}

compute_pooled_fdr_2s = function(back, exp, pit = 1){
  back_median <- median(back, na.rm = T)
  re <- sapply(exp, function(x){
    if (is.na(x)){return(NA)}
    if (x >=back_median) {
      return(sum(back >= x, na.rm = T)*pit/max(1, sum(exp >= x, na.rm = T)))
    }else{
      return(sum(back <= x, na.rm = T)*pit/max(1, sum(exp <= x, na.rm = T)))
    }})
  return(re)
}

find_discovery_wpval = function(pval, q){
  pval.adj = p.adjust(pval, method = 'BH')
  discovery_ls = lapply(q, function(q_i){
    discovery = which(pval.adj <= q_i)
  })
  return(discovery_ls)
}

find_discovery_wqval = function(pval, q){
  qval = try(qvalue(pval), silent = T)
  if(class(qval) == "try-error"){
    qval = qvalue(pval, pi0 = 1)
  }else{
    discovery_ls = lapply(q, function(q_i){
      discovery = which(qval$qvalues <= q_i)
    })
  }
  return(discovery_ls)
}

compute_fdppow = function(discovery_ls, trueidx){
  sapply(discovery_ls, function(discovery){
    fdp = sum(!discovery %in% trueidx )/max(length(discovery),1)
    pow = sum(discovery %in% trueidx)/length(trueidx)
    return(c(fdp, pow))
  })
}

compute_paired_pval_norm = function(exp, back, sd =1){
  pnorm(exp-back, mean = 0, sd = sd, lower.tail = F)
}

find_discovery_wfdr = function(fdr, q){
  discovery_ls = lapply(q, function(q_i){
    discovery = which(fdr<= q_i)
  })
  return(discovery_ls)
}
find_discovery_wlocfdr = function(loc.fdr, q ){
  if(class(loc.fdr)== "try-error"){
    discovery_ls = rep(list(which(1<0)), length(q))
  }else{
    discovery_ls = lapply(q, function(q_i){
      discovery = which(loc.fdr <= q_i)
    })
  }
  
  return(discovery_ls)
}


compute_paired_pval_nb_c = function(x,y, nb.p = 1/2){
  k = x-y
  size_hat = mean(c(x,y))
  y_ls = 0:100
  tailprob =  sapply(y_ls, function(y_i){
    # print(y)
    prod1 = pnbinom(y_i + k - 1, size = size_hat, prob = nb.p, lower.tail = F)
    prod2 = dnbinom(y_i, size = size_hat, prob = nb.p)
    if(identical(prod1, 0) | identical(prod2, 0)){
      return(0)
    }else{
      return(prod1*prod2)
    }
  })
  sum(tailprob)
}



compute_paired_pval_pois_mis = function(exp, back, sd = 1){
  logexp = log(exp + 0.01)
  logback =  log(back + 0.01)
  pnorm( logexp - logback, mean = 0, sd = sd, lower.tail = F)
}

localfdr_multrep_p = function(exp, back){
  
  Z = 0.5* (exp[,1] + exp[,2] - back[,1] - back[,2])
  z0 = 0.5* (exp[,1] + back[,1]) - 0.5 * (exp[,2] + back[,2])
  model1 <- kde(Z)
  model2 <- kde(z0)
  Z_p <- pmin(1, (1-pkde(Z, model2))/(1-pkde(Z, model1)))
  return(Z_p)
}

localfdr_multrep_p_miss = function(exp, back){
  zandz0 <- sapply(1:nrow(exp), function(j) {
    exp0 <- na.omit(exp[j,])
    back0 <- na.omit(back[j,])
    if(length(exp0 > 0 & length(back0 >0))) {
      Z <- 0.5* (exp0[1] + exp0[2] - back0[1] - back0[2])
    } else {Z <- NA}
    if(length(exp0 > 1 & length(back0 >1))) {
      z0 <- 0.5* (exp0[1] + back0[1] - exp0[2] - back0[2])
    } else {z0 <- NA}
    return(c(Z, z0))
  })
  Z <- zandz0[1,]
  z0 <- zandz0[2,]
  model1 <- kde(Z[!is.na(Z)])
  model2 <- kde(z0[!is.na(z0)])
  Z[is.na(Z)] <- 0
  z0[is.na(z0)] <- 0
  Z_p <- pmin(1, dkde(Z, model2)/dkde(Z, model1))
  return(Z_p)
}

########## distribution of difference between two nb ###########
twosample_nb = function(x,s1,s2,p1,p2)  {
  p = 0
  for (y in 0:(s2*3)) {
    p = p + pnbinom(y + x -1, s1, p1, lower.tail = F) * dnbinom(y, s2, p2 )
  }
  
  return(p)
}

twosample_nb_2s = function(x,s1,s2,p1,p2)  {
  x = abs(x)
  p = 0
  for (y in 0:(s2*3)) {
    p = p + pnbinom(y + x -1, s1, p1, lower.tail = F) * dnbinom(y, s2, p2 )
  }
  for (y in 0:(s2*3)) {
    p = p + pnbinom(y - x , s1, p1, lower.tail = T) * dnbinom(y, s2, p2 )
  }
  
  return(p)
}

trueidx = 1:(0.2*n)
positiveid = 1:(0.1 * n)
negativeid = (0.1 * n + 1):(0.2*n)

########################################################################################
######################   Result 1  #########################
######################  Fig1 multiple replicates, generated from norm   #############
########################################################################################

##############outlier################
outlier_norm = function(n, mean.norm){
  a = qnorm(0.99, mean.norm )
  rtrunc(n, spec = 'norm', mean = mean.norm, a = a)
}

set.seed(1)

#simulate mean of each point from Gaussian
mean_back = rep(0, n*3)
mean_exp = rep(0, n*3)
mean_exp[1:(length(positiveid)*3)] = rep(rnorm(length(positiveid), mean = 5, sd = 1), each = 3)
mean_exp[(length(negativeid)*3+1):(length(trueidx)*3)] = rep(rnorm(length(negativeid), mean = - 5, sd = 1), each = 3)


sim_2sided_norm_homo <- lapply(1:length(out_prop), function(i_out){
  prop =  out_prop[i_out]
  print(prop)
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rnorm(n*3, mean = 0, sd = 1) + mean_back
    score_exp <- rnorm(n*3, mean = 0, sd = 1) + mean_exp
    if (prop == 0.1) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- outlier_norm((3*n)*prop, mean.norm = 0)
      idout <- sample(1:(3*n), (3*n)*prop)
      score_exp[idout] <- outlier_norm((3*n)*prop, mean.norm = mean_exp[idout])
    }
    if (prop == 0.2) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- NA
      score_exp[sample(1:(3*n), (3*n)*prop)] <- NA
    }
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)

    ##### 2-sided clipper max with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### 2-sided clipper max with all knockoff and powerful
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper_pow = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### use BH pooled
    pval_pooled = compute_pooled_pval_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    pval_pooled[is.na(pval_pooled)] <- 1
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use qvalue pooled
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH pooled2
    fdr_pooled = compute_pooled_fdr_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wfdr(fdr = fdr_pooled, q = q_tot)
    fdppow_pooled2 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "two.sided", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired m2 misspecification
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], mu = mean(score_back[j,!(is.na(score_back[j,]))]), alternative = "two.sided")$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "two.sided", var.equal = F)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})

    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### permutation locfdr
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    fdp = list(pooled = fdppow_pooled[1,],
               pooled2 = fdppow_pooled2[1,],
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_c = fdppow_pair_c[1,],
               pair_m = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,],
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               pair_m_q = fdppow_pair_mis_q[1,],
               clippermax_1ko = fdppow_clippermax_1ko[1,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[1,],
               clippermax = fdppow_clippermax[1,],
               clipperdiff = fdppow_clipperdiff[1,],
               clippermax_3ko = fdppow_clippermax_3ko[1,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[1,],
               clipper_pow = fdppow_clipper_pow[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])

    pow = list(pooled = fdppow_pooled[2,],
               pooled2 = fdppow_pooled2[2,],
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_c = fdppow_pair_c[2,],
               pair_m = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,],
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               pair_m_q = fdppow_pair_mis_q[2,],
               clippermax_1ko = fdppow_clippermax_1ko[2,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[2,],
               clippermax = fdppow_clippermax[2,],
               clipperdiff = fdppow_clipperdiff[2,],
               clippermax_3ko = fdppow_clippermax_3ko[2,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[2,],
               clipper_pow = fdppow_clipper_pow[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])


    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)

  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')

  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr,
                             q = rep(q_tot, length(methods)),
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_2sided_norm_homo) = out_prop
saveRDS(sim_2sided_norm_homo, 'sim_2sided_norm_homo.rds')




set.seed(1)

###simulate mean of each point from Gaussian
mean_back = rep(rnorm(n, mean = 0, sd = 2), each = 3)
mean_exp = mean_back
mean_exp[1:(length(positiveid)*3)] = rep(rnorm(length(positiveid), mean = 5, sd = 1), each = 3)
mean_exp[(length(negativeid)*3+1):(length(trueidx)*3)] = rep(rnorm(length(negativeid), mean = - 5, sd = 1), each = 3)

sim_2sided_norm_hete <- lapply(1:length(out_prop), function(i_out){

  prop =  out_prop[i_out]
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rnorm(n*3, mean = 0, sd = 1) + mean_back
    score_exp <- rnorm(n*3, mean = 0, sd = 1) + mean_exp
    if (prop == 0.1) {
      idout <- sample(1:(3*n), (3*n)*prop)
      score_back[idout] <- outlier_norm((3*n)*prop, mean.norm = mean_back[idout])
      idout <- sample(1:(3*n), (3*n)*prop)
      score_exp[idout] <- outlier_norm((3*n)*prop, mean.norm = mean_exp[idout])
    }
    if (prop == 0.2) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- NA
      score_exp[sample(1:(3*n), (3*n)*prop)] <- NA
    }
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)

    ##### 2-sided clipper max with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### 2-sided clipper max with all knockoff and powerful
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper_pow = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH pooled
    pval_pooled = compute_pooled_pval_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH pooled2
    fdr_pooled = compute_pooled_fdr_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wfdr(fdr = fdr_pooled, q = q_tot)
    fdppow_pooled2 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "two.sided", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired m2 misspecification
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], mu = mean(score_back[j,!(is.na(score_back[j,]))]), alternative = "two.sided")$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "two.sided", var.equal = F)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})

    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### permutation locfdr
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    fdp = list(pooled = fdppow_pooled[1,],
               pooled2 = fdppow_pooled2[1,],
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_c = fdppow_pair_c[1,],
               pair_m = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,],
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               pair_m_q = fdppow_pair_mis_q[1,],
               clippermax_1ko = fdppow_clippermax_1ko[1,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[1,],
               clippermax = fdppow_clippermax[1,],
               clipperdiff = fdppow_clipperdiff[1,],
               clippermax_3ko = fdppow_clippermax_3ko[1,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[1,],
               clipper_pow = fdppow_clipper_pow[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])

    pow = list(pooled = fdppow_pooled[2,],
               pooled2 = fdppow_pooled2[2,],
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_c = fdppow_pair_c[2,],
               pair_m = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,],
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               pair_m_q = fdppow_pair_mis_q[2,],
               clippermax_1ko = fdppow_clippermax_1ko[2,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[2,],
               clippermax = fdppow_clippermax[2,],
               clipperdiff = fdppow_clipperdiff[2,],
               clippermax_3ko = fdppow_clippermax_3ko[2,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[2,],
               clipper_pow = fdppow_clipper_pow[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])


    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)

  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')

  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr,
                             q = rep(q_tot, length(methods)),
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})
names(sim_2sided_norm_hete) = out_prop
saveRDS(sim_2sided_norm_hete, 'sim_2sided_norm_hete.rds')


######################  Fig1 multiple replicates, generated from Pois   #############
########################################################################################

set.seed(1)

outlier_pois = function(n, mean.pois){
  a = qpois(0.99, mean.pois )
  rtrunc(n, spec = 'pois', lambda = mean.pois, a = a)
}


#simulate mean of each point from Poisson
mean_back = rep(20, n*3)
mean_exp = rep(20, n*3)
mean_exp[1:(length(positiveid)*3)] = rep(rpois(length(positiveid), lambda = 40), each = 3)
mean_exp[(length(negativeid)*3+1):(length(trueidx)*3)] = rep(rpois(length(negativeid), lambda = 5), each = 3)


sim_2sided_pois_homo <- lapply(1:length(out_prop), function(i_out){

  prop =  out_prop[i_out]
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rpois(n*3, lambda = mean_back)
    score_exp <- rpois(n*3, lambda = mean_exp)
    if (prop == 0.1) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- outlier_pois((3*n)*prop, mean.pois =  20)
      idout <- sample(1:(3*n), (3*n)*prop)
      score_exp[idout] <- outlier_pois((3*n)*prop, mean.pois = mean_exp[idout])
    }
    if (prop == 0.2) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- NA
      score_exp[sample(1:(3*n), (3*n)*prop)] <- NA
    }
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)

    ##### 2-sided clipper max with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### 2-sided clipper max with all knockoff and powerful
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper_pow = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### use BH pooled
    pval_pooled = compute_pooled_pval_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH pooled2
    fdr_pooled = compute_pooled_fdr_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wfdr(fdr = fdr_pooled, q = q_tot)
    fdppow_pooled2 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "two.sided")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return((poisson.test(sum(score_exp[j,], na.rm = T), sum(score_back[j,], na.rm = T)* sum(!is.na(score_exp[j,]))/sum(!is.na(score_back[j,])), alternative = "two.sided")$p.value))
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (length(unique(na.omit(score_exp[j,])))<2|length(unique(na.omit(score_back[j,])))<2) {
        return(1)} else {
          return(t.test(-log(score_exp[j,!(is.na(score_exp[j,]))]+0.01), -log(score_back[j,!(is.na(score_back[j,]))]+0.01), alternative = "two.sided", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})

    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### permutation locfdr
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    fdp = list(pooled = fdppow_pooled[1,],
               pooled2 = fdppow_pooled2[1,],
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_c = fdppow_pair_c[1,],
               pair_m = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,],
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               pair_m_q = fdppow_pair_mis_q[1,],
               clippermax_1ko = fdppow_clippermax_1ko[1,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[1,],
               clippermax = fdppow_clippermax[1,],
               clipperdiff = fdppow_clipperdiff[1,],
               clippermax_3ko = fdppow_clippermax_3ko[1,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[1,],
               clipper_pow = fdppow_clipper_pow[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,],
               pooled2 = fdppow_pooled2[2,],
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_c = fdppow_pair_c[2,],
               pair_m = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,],
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               pair_m_q = fdppow_pair_mis_q[2,],
               clippermax_1ko = fdppow_clippermax_1ko[2,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[2,],
               clippermax = fdppow_clippermax[2,],
               clipperdiff = fdppow_clipperdiff[2,],
               clippermax_3ko = fdppow_clippermax_3ko[2,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[2,],
               clipper_pow = fdppow_clipper_pow[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)

  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')

  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr,
                             q = rep(q_tot, length(methods)),
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_2sided_pois_homo) = out_prop
saveRDS(sim_2sided_pois_homo, 'sim_2sided_pois_homo.rds')



set.seed(1)
mean_back = rep(rnbinom(n, size = 20, prob = 0.5), each = 3)
mean_exp = mean_back
mean_exp[1:(length(positiveid)*3)] = rep(rpois(length(positiveid), lambda = 40), each = 3)
mean_exp[(length(negativeid)*3+1):(length(trueidx)*3)] = rep(rpois(length(negativeid), lambda = 5), each = 3)

sim_2sided_pois_hete <- lapply(1:length(out_prop), function(i_out){

  prop =  out_prop[i_out]
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rpois(n*3, lambda = mean_back)
    score_exp <- rpois(n*3, lambda = mean_exp)
    if (prop == 0.1) {
      idout <- sample(1:(3*n), (3*n)*prop)
      score_back[idout] <- outlier_pois((3*n)*prop, mean.pois = mean_back[idout])
      idout <- sample(1:(3*n), (3*n)*prop)
      score_exp[idout] <- outlier_pois((3*n)*prop, mean.pois = mean_exp[idout])
    }
    if (prop == 0.2) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- NA
      score_exp[sample(1:(3*n), (3*n)*prop)] <- NA
    }
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)

    ##### 2-sided clipper max with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 1,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper max with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### 2-sided clipper diff with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'diff',
                                nknockoff = 3,
                                ifpowerful = F,
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)


    ##### 2-sided clipper max with all knockoff and powerful
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot,
                                contrastScore_method = 'max',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper_pow = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH pooled
    pval_pooled = compute_pooled_pval_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH pooled2
    fdr_pooled = compute_pooled_fdr_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wfdr(fdr = fdr_pooled, q = q_tot)
    fdppow_pooled2 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired correct specification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "two.sided")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return((poisson.test(sum(score_exp[j,], na.rm = T), sum(score_back[j,], na.rm = T)* sum(!is.na(score_exp[j,]))/sum(!is.na(score_back[j,])), alternative = "two.sided")$p.value))
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (length(unique(na.omit(score_exp[j,])))<2|length(unique(na.omit(score_back[j,])))<2) {
        return(1)} else {
          return(t.test(-log(score_exp[j,!(is.na(score_exp[j,]))]+0.01), -log(score_back[j,!(is.na(score_back[j,]))]+0.01), alternative = "two.sided", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})

    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### permutation locfdr
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    fdp = list(pooled = fdppow_pooled[1,],
               pooled2 = fdppow_pooled2[1,],
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_c = fdppow_pair_c[1,],
               pair_m = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,],
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               pair_m_q = fdppow_pair_mis_q[1,],
               clippermax_1ko = fdppow_clippermax_1ko[1,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[1,],
               clippermax = fdppow_clippermax[1,],
               clipperdiff = fdppow_clipperdiff[1,],
               clippermax_3ko = fdppow_clippermax_3ko[1,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[1,],
               clipper_pow = fdppow_clipper_pow[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,],
               pooled2 = fdppow_pooled2[2,],
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_c = fdppow_pair_c[2,],
               pair_m = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,],
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               pair_m_q = fdppow_pair_mis_q[2,],
               clippermax_1ko = fdppow_clippermax_1ko[2,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[2,],
               clippermax = fdppow_clippermax[2,],
               clipperdiff = fdppow_clipperdiff[2,],
               clippermax_3ko = fdppow_clippermax_3ko[2,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[2,],
               clipper_pow = fdppow_clipper_pow[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)

  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')

  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr,
                             q = rep(q_tot, length(methods)),
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_2sided_pois_hete) = out_prop
saveRDS(sim_2sided_pois_hete, 'sim_2sided_pois_hete.rds')

######################  Fig1 single replicates, generated from NB   #############
########################################################################################
nb.p = 1/2
##############outlier################
outlier_nb = function(n, size.nb){
  a = qnbinom(0.99, size.nb, nb.p )
  rtrunc(n, spec = 'nbinom', size = size.nb,prob = nb.p, a = a)
}

set.seed(1)

mean_back = rep(30, n*3)
mean_exp = rep(30, n*3)
mean_exp[1:(length(positiveid)*3)] = rep(rnbinom(length(positiveid), size = 70, prob = nb.p), each = 3)
mean_exp[(length(negativeid)*3+1):(length(trueidx)*3)] = rep(pmax(1,rnbinom(length(negativeid), size = 7, prob = nb.p)), each = 3)

sim_2sided_nb_homo <- lapply(1:length(out_prop), function(i_out){
  prop = out_prop[i_out]

  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rnbinom(n*3, size = mean_back, prob = nb.p)
    score_exp <- rnbinom(n*3, size = mean_exp, prob = nb.p)
    if (prop == 0.1) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- outlier_nb((3*n)*prop, size.nb =  20)
      idout <- sample(1:(3*n), (3*n)*prop)
      score_exp[idout] <- outlier_nb((3*n)*prop, size.nb = mean_exp[idout])
    }
    if (prop == 0.2) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- NA
      score_exp[sample(1:(3*n), (3*n)*prop)] <- NA
    }
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)

    ##### 2-sided clipper max with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                nknockoff = 1,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper diff with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'diff', 
                                nknockoff = 1,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper max with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper diff with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'diff', 
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper max with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                nknockoff = 3,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper diff with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'diff', 
                                nknockoff = 3,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### 2-sided clipper max with all knockoff and powerful
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper_pow = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH pooled2
    fdr_pooled = compute_pooled_fdr_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wfdr(fdr = fdr_pooled, q = q_tot)
    fdppow_pooled2 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      score_exp0 <- na.omit(score_exp[j,])
      score_back0 <- na.omit(score_back[j,])
      if (length(score_exp0) == 0|length(score_back0) == 0) {return(1)} else {
        y <- c(score_exp0, score_back0)
        indi <- factor(c(rep(1, length(score_exp0)), rep(0,length(score_exp0))))
        re.glm <- try(summary(glm.nb(y~indi)), silent = T)
        if(class(re.glm)[1]=="try-error") {
          return(1)
        }else{
        return(re.glm$coefficients[2,4])
        }
      }
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j) {
      score_exp0 <- na.omit(score_exp[j,])
      score_back0 <- na.omit(score_back[j,])
      if (length(score_exp0) == 0|length(score_back0) == 0) {return(1)} else {
        re.glm <- try(glm.nb(score_back0~1), silent = T)
        if(class(re.glm)[1]=="try-error") {
          return(1)
        }else{
          mean_hat <- exp(re.glm$coefficients)
          size_hat <- re.glm$theta
          return(2*min(pnbinom(sum(score_exp0), mu = length(score_exp0)*mean_hat, size = size_hat*length(score_exp0), lower.tail =T),
                       pnbinom(sum(score_exp0), mu = length(score_exp0)*mean_hat, size = size_hat*length(score_exp0), lower.tail =F)))
        }
      }
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use model misspecification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "two.sided")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})

    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### permutation locfdr
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    fdp = list(pooled = fdppow_pooled[1,],
               pooled2 = fdppow_pooled2[1,],
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_c = fdppow_pair_c[1,],
               pair_m = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,],
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               pair_m_q = fdppow_pair_mis_q[1,],
               clippermax_1ko = fdppow_clippermax_1ko[1,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[1,],
               clippermax = fdppow_clippermax[1,],
               clipperdiff = fdppow_clipperdiff[1,],
               clippermax_3ko = fdppow_clippermax_3ko[1,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[1,],
               clipper_pow = fdppow_clipper_pow[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,],
               pooled2 = fdppow_pooled2[2,],
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_c = fdppow_pair_c[2,],
               pair_m = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,],
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               pair_m_q = fdppow_pair_mis_q[2,],
               clippermax_1ko = fdppow_clippermax_1ko[2,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[2,],
               clippermax = fdppow_clippermax[2,],
               clipperdiff = fdppow_clipperdiff[2,],
               clippermax_3ko = fdppow_clippermax_3ko[2,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[2,],
               clipper_pow = fdppow_clipper_pow[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)

  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')

  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr,
                             q = rep(q_tot, length(methods)),
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_2sided_nb_homo) = out_prop
saveRDS(sim_2sided_nb_homo, 'sim_2sided_nb_homo.rds')


set.seed(1)
mean_back = rep(rnbinom(n, size = 10, prob = 0.25), each = 3)
mean_exp = mean_back
mean_exp[1:(length(positiveid)*3)] = rep(rnbinom(length(positiveid), size = 70, prob = nb.p), each = 3)
mean_exp[(length(negativeid)*3+1):(length(trueidx)*3)] = rep(pmax(1,rnbinom(length(negativeid), size = 7, prob = nb.p)), each = 3)

sim_2sided_nb_hete <- lapply(1:length(out_prop), function(i_out){
  prop = out_prop[i_out]

  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rnbinom(n*3, size = mean_back, prob = nb.p)
    score_exp <- rnbinom(n*3, size = mean_exp, prob = nb.p)
    if (prop == 0.1) {
      idout <- sample(1:(3*n), (3*n)*prop)
      score_back[idout] <- outlier_nb((3*n)*prop, size.nb = mean_back[idout])
      idout <- sample(1:(3*n), (3*n)*prop)
      score_exp[idout] <- outlier_nb((3*n)*prop, size.nb = mean_exp[idout])
    }
    if (prop == 0.2) {
      score_back[sample(1:(3*n) , (3*n)*prop)] <- NA
      score_exp[sample(1:(3*n), (3*n)*prop)] <- NA
    }
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)

    ##### 2-sided clipper max with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                nknockoff = 1,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper diff with 1 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'diff', 
                                nknockoff = 1,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_1ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper max with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper diff with all knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'diff', 
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper max with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                nknockoff = 3,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### 2-sided clipper diff with 3 knockoff
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'diff', 
                                nknockoff = 3,
                                ifpowerful = F, 
                                seed = it)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_3ko = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### 2-sided clipper max with all knockoff and powerful
    re_clipper <- clipper2sided(score_exp = score_exp,
                                score_back = score_back,
                                FDR = q_tot, 
                                contrastScore_method = 'max', 
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper_pow = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH pooled2
    fdr_pooled = compute_pooled_fdr_2s(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wfdr(fdr = fdr_pooled, q = q_tot)
    fdppow_pooled2 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      score_exp0 <- na.omit(score_exp[j,])
      score_back0 <- na.omit(score_back[j,])
      if (length(score_exp0) == 0|length(score_back0) == 0) {return(1)} else {
        y <- c(score_exp0, score_back0)
        indi <- factor(c(rep(1, length(score_exp0)), rep(0,length(score_exp0))))
        re.glm <- try(summary(glm.nb(y~indi)), silent = T)
        if(class(re.glm)[1]=="try-error") {
          return(1)
        }else{
          return(re.glm$coefficients[2,4])
        }
      }
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j) {
      score_exp0 <- na.omit(score_exp[j,])
      score_back0 <- na.omit(score_back[j,])
      if (length(score_exp0) == 0|length(score_back0) == 0) {return(1)} else {
        re.glm <- try(glm.nb(score_back0~1), silent = T)
        if(class(re.glm)[1]=="try-error") {
          return(1)
        }else{
          mean_hat <- exp(re.glm$coefficients)
          size_hat <- re.glm$theta
          return(2*min(pnbinom(sum(score_exp0), mu = length(score_exp0)*mean_hat, size = size_hat*length(score_exp0), lower.tail =T),
                       pnbinom(sum(score_exp0), mu = length(score_exp0)*mean_hat, size = size_hat*length(score_exp0), lower.tail =F)))
        }
       }
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use model misspecification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "two.sided")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})

    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### permutation locfdr
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    fdp = list(pooled = fdppow_pooled[1,],
               pooled2 = fdppow_pooled2[1,],
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_c = fdppow_pair_c[1,],
               pair_m = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,],
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               pair_m_q = fdppow_pair_mis_q[1,],
               clippermax_1ko = fdppow_clippermax_1ko[1,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[1,],
               clippermax = fdppow_clippermax[1,],
               clipperdiff = fdppow_clipperdiff[1,],
               clippermax_3ko = fdppow_clippermax_3ko[1,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[1,],
               clipper_pow = fdppow_clipper_pow[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,],
               pooled2 = fdppow_pooled2[2,],
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_c = fdppow_pair_c[2,],
               pair_m = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,],
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               pair_m_q = fdppow_pair_mis_q[2,],
               clippermax_1ko = fdppow_clippermax_1ko[2,],
               clipperdiff_1ko = fdppow_clipperdiff_1ko[2,],
               clippermax = fdppow_clippermax[2,],
               clipperdiff = fdppow_clipperdiff[2,],
               clippermax_3ko = fdppow_clippermax_3ko[2,],
               clipperdiff_3ko = fdppow_clipperdiff_3ko[2,],
               clipper_pow = fdppow_clipper_pow[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)

  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')

  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr,
                             q = rep(q_tot, length(methods)),
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_2sided_nb_hete) = out_prop
saveRDS(sim_2sided_nb_hete, 'sim_2sided_nb_hete.rds')

########################################################################################





