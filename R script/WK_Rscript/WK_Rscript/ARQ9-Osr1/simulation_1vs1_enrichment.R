remove(list = ls())
library(parallel)
library(locfdr)
library(qvalue)
source('clipper102520.R', echo=F)
ncores = detectCores() - 2
iter = 200

N <- c(1000, 10000)
q_tot = (1:10)/100
####### functions ######

compute_pooled_pval = function(back, exp){
  re = sapply(exp, function(x){
    mean(back >= x)
  })
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


find_discovery_wlocfdr = function(loc.fdr, q ){
  if(class(loc.fdr)== "try-error"){
    discovery_ls = rep(list(which(1<0)), length(q))
  }else{
    discovery_ls = lapply(q, function(q_i){
      discovery = which(loc.fdr$fdr <= q_i)
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



########################################################################################
######################   Result 1  #########################
######################  Fig1 single replicates, generated from norm   #############
########################################################################################


set.seed(1)
trueidx_tot = 1:(0.1*max(N))

#simulate mean of each point from Gaussian
mean_back_tot = rep(0, max(N))
mean_exp_tot = rep(0, max(N))
mean_exp_tot[trueidx_tot] = rnorm(length(trueidx_tot), mean = 5, sd = 1)


sim_singlerep_norm_homo <- lapply(1:length(N), function(i_N){
  
  n =  N[i_N]
  trueidx = 1:(0.1*n)
  idx = c(trueidx, (0.1*max(N)+ 1):(0.1*max(N) + 0.9*n) )
  
  mean_back = mean_back_tot[idx]
  mean_exp = mean_exp_tot[idx]
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnorm(n, mean = 0, sd = 1) + mean_back
    score_exp <- rnorm(n, mean = 0, sd = 1) + mean_exp
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = score_back, exp = score_exp)
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use qvalue pooled
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 invalid pvalue
    p <-  mapply(function(exp, back){
      1 - pnorm(exp, mean = back, sd = sd(score_back))
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use qvalue pooled
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- mapply(function(exp, back){
      compute_paired_pval_norm(exp, back, sd = sqrt(var(score_exp) + var(score_back)))
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    ##### use qvalue pooled
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                          score_back = score_back, 
                          FDR = q_tot, 
                          importanceScore_method = 'diff', 
                          FDR_control_method = 'BC',
                          ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper diff BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper max BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### emp null
    loc.fdr <- try(locfdr((score_exp-score_back)/sqrt(var(score_exp-score_back)),plot = 0 ), silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    fdp = list(pooled = fdppow_pooled[1,], 
            pair_2as1 = fdppow_pair_2as1[1,],
            pair_mis = fdppow_pair_mis[1,],
            pooled_q = fdppow_pooled_q[1,], 
            pair_2as1_q = fdppow_pair_2as1_q[1,],
            pair_mis_q = fdppow_pair_mis_q[1,],
            clipper = fdppow_clipper[1,],
            clippermax_bc = fdppow_clippermax_bc[1,],
            clipperdiff_bc = fdppow_clipperdiff_bc[1,],
            clippermax_bh = fdppow_clippermax_bh[1,],
            clipperdiff_bh = fdppow_clipperdiff_bh[1,],
            locfdremp = fdppow_locfdremp[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
            pair_2as1 = fdppow_pair_2as1[2,],
            pair_mis = fdppow_pair_mis[2,],
            pooled_q = fdppow_pooled_q[2,], 
            pair_2as1_q = fdppow_pair_2as1_q[2,],
            pair_mis_q = fdppow_pair_mis_q[2,],
            clipper = fdppow_clipper[2,],
            clippermax_bc = fdppow_clippermax_bc[2,],
            clipperdiff_bc = fdppow_clipperdiff_bc[2,],
            clippermax_bh = fdppow_clippermax_bh[2,],
            clipperdiff_bh = fdppow_clipperdiff_bh[2,],
            locfdremp = fdppow_locfdremp[2,])

    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply(fdppow_ls, '[[', 'fdp')
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

names(sim_singlerep_norm_homo) = N
saveRDS(sim_singlerep_norm_homo, 'sim_singlerep_norm_homo.rds')
sim_singlerep_norm_homo = readRDS('sim_singlerep_norm_homo.rds')


set.seed(1)
trueidx_tot = 1:(0.1*max(N))

#simulate mean of each point from Gaussian
mean_back_tot = rnorm(max(N), mean = 0, sd = 2)
mean_exp_tot = mean_back_tot
mean_exp_tot[trueidx_tot] = rnorm(length(trueidx_tot), mean = 5, sd = 1)

sim_singlerep_norm_hete <- lapply(1:length(N), function(i_N){
  
  n =  N[i_N]
  trueidx = 1:(0.1*n)
  idx = c(trueidx, (0.1*max(N)+ 1):(0.1*max(N) + 0.9*n) )
  
  mean_back = mean_back_tot[idx]
  mean_exp = mean_exp_tot[idx]
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnorm(n, mean = 0, sd = 1) + mean_back
    score_exp <- rnorm(n, mean = 0, sd = 1) + mean_exp
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = score_back, exp = score_exp)
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use qvalue pooled
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 invalid pvalue
    p <-  mapply(function(exp, back){
      1 - pnorm(exp, mean = back, sd = sd(score_back))
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- mapply(function(exp, back){
      compute_paired_pval_norm(exp, back, sd = sqrt(var(score_exp) + var(score_back)))
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper diff BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper max BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### emp null
    loc.fdr <- try(locfdr((score_exp-score_back)/sqrt(var(score_exp-score_back)),plot = 0 ), silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               clippermax_bh = fdppow_clippermax_bh[1,],
               clipperdiff_bh = fdppow_clipperdiff_bh[1,],
               locfdremp = fdppow_locfdremp[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               clippermax_bh = fdppow_clippermax_bh[2,],
               clipperdiff_bh = fdppow_clipperdiff_bh[2,],
               locfdremp = fdppow_locfdremp[2,])
    
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
names(sim_singlerep_norm_hete) = N
saveRDS(sim_singlerep_norm_hete, 'sim_singlerep_norm_hete.rds')


######################  Fig1 single replicates, generated from Pois   #############
########################################################################################

set.seed(1)
trueidx_tot = 1:(0.1*max(N))

# simulate mean of each point from Gaussian
mean_back_tot = rep(20, max(N))
mean_exp_tot = rep(20, max(N))
mean_exp_tot[trueidx_tot] = rpois(length(trueidx_tot), lambda = 40)


sim_singlerep_pois_homo <- lapply(1:length(N), function(i_N){

  n =  N[i_N]
  trueidx = 1:(0.1*n)
  idx = c(trueidx, (0.1*max(N)+ 1):(0.1*max(N) + 0.9*n) )

  mean_back = mean_back_tot[idx]
  mean_exp = mean_exp_tot[idx]

  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){

    set.seed(it)
    score_back <- rpois(n, lambda = mean_back)
    score_exp <- rpois(n, lambda = mean_exp)

    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = score_back, exp = score_exp)
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired 2as1
    p <-  mapply(function(exp, back){
      1 - ppois(exp-1, lambda = back)
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    

    ##### use BH paired m2 misspecification as variance 1
    p <- mapply(function(exp, back){
      compute_paired_pval_pois_mis(exp, back, sd = sqrt( var(log(score_exp + 0.01)) + var(log(score_back + 0.01)) ))
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper diff BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper max BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### emp null
    loc.fdr <- try(locfdr((score_exp-score_back)/sqrt(var(score_exp-score_back)),plot = 0 ), silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               clippermax_bh = fdppow_clippermax_bh[1,],
               clipperdiff_bh = fdppow_clipperdiff_bh[1,],
               locfdremp = fdppow_locfdremp[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               clippermax_bh = fdppow_clippermax_bh[2,],
               clipperdiff_bh = fdppow_clipperdiff_bh[2,],
               locfdremp = fdppow_locfdremp[2,])
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

names(sim_singlerep_pois_homo) = N
saveRDS(sim_singlerep_pois_homo, 'sim_singlerep_pois_homo.rds')
sim_singlerep_pois_homo = readRDS('sim_singlerep_pois_homo.rds')


set.seed(1)
trueidx_tot = 1:(0.1*max(N))
mean_back_tot = rnbinom(max(N), size = 20, prob = 0.5)
mean_exp_tot = mean_back_tot
mean_exp_tot[trueidx_tot] = rpois(length(trueidx_tot), lambda = 40)


sim_singlerep_pois_hete <- lapply(1:length(N), function(i_N){
  
  n =  N[i_N]
  trueidx = 1:(0.1*n)
  idx = c(trueidx, (0.1*max(N)+ 1):(0.1*max(N) + 0.9*n) )
  
  mean_back = mean_back_tot[idx]
  mean_exp = mean_exp_tot[idx]
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rpois(n, lambda = mean_back)
    score_exp <- rpois(n, lambda = mean_exp)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = score_back, exp = score_exp)
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired 2as1
    p <-  mapply(function(exp, back){
      1 - ppois(exp-1, lambda = back)
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    
    ##### use BH paired m2 misspecification as variance 1
    p <- mapply(function(exp, back){
      compute_paired_pval_pois_mis(exp, back, sd = sqrt( var(log(score_exp + 0.01)) + var(log(score_back + 0.01)) ))
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper diff BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper max BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### emp null
    loc.fdr <- try(locfdr((score_exp-score_back)/sqrt(var(score_exp-score_back)),plot = 0 ), silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               clippermax_bh = fdppow_clippermax_bh[1,],
               clipperdiff_bh = fdppow_clipperdiff_bh[1,],
               locfdremp = fdppow_locfdremp[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               clippermax_bh = fdppow_clippermax_bh[2,],
               clipperdiff_bh = fdppow_clipperdiff_bh[2,],
               locfdremp = fdppow_locfdremp[2,])
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

names(sim_singlerep_pois_hete) = N
saveRDS(sim_singlerep_pois_hete, 'sim_singlerep_pois_hete.rds')
sim_singlerep_pois_hete = readRDS('sim_singlerep_pois_hete.rds')

######################  Fig1 single replicates, generated from NB   #############
########################################################################################
set.seed(1)
trueidx_tot = 1:(0.1*max(N))
nb.p = 1/2


mean_back_tot = rep(20, max(N))
mean_exp_tot = rep(20, max(N))
mean_exp_tot[trueidx_tot] = rnbinom(length(trueidx_tot), size = 45, prob = nb.p)

sim_singlerep_nb_homo <- lapply(1:length(N), function(i_N){
  
  n =  N[i_N]
  trueidx = 1:(0.1*n)
  idx = c(trueidx, (0.1*max(N)+ 1):(0.1*max(N) + 0.9*n) )
  
  mean_back = mean_back_tot[idx]
  mean_exp = mean_exp_tot[idx]
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnbinom(n, size = mean_back, prob = nb.p)
    score_exp <- rnbinom(n, size = mean_exp, prob = nb.p)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = score_back, exp = score_exp)
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 invalid pvalue 
    p <-  mapply(function(exp, back){
      1 - pnbinom(exp-1, size = back, nb.p)
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### use BH paired m2 by misspecifying NB as poisson
    p <- sapply(1:n, function(j){
      poisson.test(c(score_exp[j],score_back[j]),c(1,1), alternative = "greater")$p.value
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper diff BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper max BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    
    
    ##### emp null
    loc.fdr <- try(locfdr((score_exp-score_back)/sqrt(var(score_exp-score_back)),plot = 0 ), silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               clippermax_bh = fdppow_clippermax_bh[1,],
               clipperdiff_bh = fdppow_clipperdiff_bh[1,],
               locfdremp = fdppow_locfdremp[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               clippermax_bh = fdppow_clippermax_bh[2,],
               clipperdiff_bh = fdppow_clipperdiff_bh[2,],
               locfdremp = fdppow_locfdremp[2,])
    
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

names(sim_singlerep_nb_homo) = N
saveRDS(sim_singlerep_nb_homo, 'sim_singlerep_nb_homo.rds')
sim_singlerep_nb_homo = readRDS('sim_singlerep_nb_homo.rds')


set.seed(1)
trueidx_tot = 1:(0.1*max(N))

# rg = 1:50
# plot(rg, dnbinom(rg, size = 20, prob = 0.5))
# points(rg, dnbinom(rg, size = 50, prob = 0.5))
mean_back_tot = rnbinom(max(N), size = 20, prob = nb.p)
mean_exp_tot = mean_back_tot
mean_exp_tot[trueidx_tot] = rnbinom(length(trueidx_tot), size = 45, prob = nb.p)

sim_singlerep_nb_hete <- lapply(1:length(N), function(i_N){
  
  n =  N[i_N]
  trueidx = 1:(0.1*n)
  idx = c(trueidx, (0.1*max(N)+ 1):(0.1*max(N) + 0.9*n) )
  
  mean_back = mean_back_tot[idx]
  mean_exp = mean_exp_tot[idx]
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnbinom(n, size = mean_back, prob = nb.p)
    score_exp <- rnbinom(n, size = mean_exp, prob = nb.p)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = score_back, exp = score_exp)
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 invalid pvalue 
    p <-  mapply(function(exp, back){
      1 - pnbinom(exp-1, size = back, nb.p)
    }, exp = score_exp, back = score_back, SIMPLIFY = T)
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### use BH paired m2 by misspecifying NB as poisson
    p <- sapply(1:n, function(j){
      poisson.test(c(score_exp[j],score_back[j]),c(1,1), alternative = "greater")$p.value
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### clipper diff BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'diff', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper max BH
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method = 'max', 
                                FDR_control_method = 'BH',
                                ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bh = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    ##### emp null
    loc.fdr <- try(locfdr((score_exp-score_back)/sqrt(var(score_exp-score_back)),plot = 0 ), silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               clippermax_bh = fdppow_clippermax_bh[1,],
               clipperdiff_bh = fdppow_clipperdiff_bh[1,],
               locfdremp = fdppow_locfdremp[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               clippermax_bh = fdppow_clippermax_bh[2,],
               clipperdiff_bh = fdppow_clipperdiff_bh[2,],
               locfdremp = fdppow_locfdremp[2,])
    
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

names(sim_singlerep_nb_hete) = N
saveRDS(sim_singlerep_nb_hete, 'sim_singlerep_nb_hete.rds')
sim_singlerep_nb_hete = readRDS('sim_singlerep_nb_hete.rds')

########################################################################################





