rm(list = ls())
library(openxlsx)
library(parallel)
# library(rgl)
ncores = detectCores() - 1
sequest = read.xlsx(xlsxFile ='~/Dropbox/database search in proteomics/data/Archaea database Search Results - updated/Archaea Sequest Targets 100_FDR.xlsx')
sequest_decoy = read.xlsx(xlsxFile ='~/Dropbox/database search in proteomics/data/Archaea database Search Results - updated/Archaea Sequest DECOYS 100_FDR.xlsx')

FDR_ls = (1:10)/100




run_clipper_proteomics = function(searchalg){
  
  dat = eval(parse(text = searchalg))
  truematch = which(dat$Label %in% c('Contaminant','Standard'))

  dat_decoy = eval(parse(text = paste0(searchalg,'_decoy')))
  
  scoreColTitle = c(
    'Percolator.PEP')
  # scoreColTitle = c(
  #   'XCorr')
  scannumColTitle = rep('ScanNum',1)
  sequenceColTitle = c(
    'Sequence')
  
  fileColTitle = c("Spectrum.File")
  fdppow_ls_qval = sapply(FDR_ls,function(FDR_i){
    discovery = which(dat[, 'Percolator.q-Value'] <= FDR_i)
    c(sum(!discovery %in% truematch)/max(1, length(discovery)), mean(truematch %in% discovery))
  }) 
  saveRDS( fdppow_ls_qval, file = paste0('realdata_proteomics_',searchalg,'_qvalue_fdppow_ls.rds'))
  
  
  
  
  
  query = paste( dat[,fileColTitle], dat[,scannumColTitle],sep = ':')
  match = paste( dat[,fileColTitle], dat[,scannumColTitle], dat[,sequenceColTitle],sep = ':')
  score = -log(dat[, scoreColTitle]+ 0.01)
  
  dat = cbind.data.frame(query = query, match = match, score = score, stringsAsFactors = F)
  
  query_decoy = paste(dat_decoy[, fileColTitle], dat_decoy[,scannumColTitle],sep = ':')
  match_decoy = paste(dat_decoy[, fileColTitle], dat_decoy[,scannumColTitle], dat_decoy[,sequenceColTitle],sep = ':')
  score_decoy = -log(dat_decoy[, scoreColTitle]+ 0.01)
  
  ### get rid of the false decoy matches and queries
  query_decoy = query_decoy[!match_decoy %in% match]
  match_decoy = match_decoy[!match_decoy %in% match]
  score_decoy = score_decoy[!match_decoy %in% match]
  
  ### re-order decoy queries and scores to match the target output
  score_decoy = score_decoy[match(dat$query, query_decoy)]
  match_decoy = match_decoy[match(dat$query, query_decoy)]
  query_decoy = query_decoy[match(dat$query, query_decoy)]
  mean(query_decoy %in% query)
  length(score_decoy)
  
 
  
  
  
  ### run clipper with diff contrast score
  re = clipper1sided(score_exp = score,
               score_back = score_decoy,
               FDR = FDR_ls,
               ifpowerful = F)
  re = sapply(re$results, function(x){
    c(mean(! x$discovery %in% truematch),mean(truematch  %in%  x$discovery))
  })
  saveRDS(re, file = paste0('realdata_proteomics_',searchalg,'_clipper_fdppow_ls.rds'))
  ### use Jackknife estimate for estimating FDR
  n = length(score)
  idx_tot = 1:n
  fdppow_ls = mclapply(1:n, function(i){
    idx_i = idx_tot[-i]
    re_i = clipper1sided(score_exp = score[-i], 
                   score_back = score_decoy[-i], 
                   FDR = FDR_ls,
                   ifpowerful = F)
    fdppow_ls_i = sapply(re_i$results, function(x){
      # c(sum(!x$discovery %in% truematch)/max(1, length(x$discovery)), mean(truematch %in% x$discovery))
      
      c(sum(!idx_i[x$discovery] %in% truematch)/max(1, length(x$discovery)), mean(truematch %in% idx_i[x$discovery]))
    })
    rownames(fdppow_ls_i) = c('fdp','pow')
    colnames(fdppow_ls_i) = FDR_ls
    return(fdppow_ls_i)
  },mc.cores = 15)
  
  saveRDS(fdppow_ls, file = paste0('realdata_proteomics_',searchalg,'_clipper_jackknife_fdppow_ls.rds'))
  return(NULL)
}

run_clipper_proteomics('sequest')
### baseline fdp and power
# discovery = which(dat[,scoreColTitle] <= 0.01)
# c(mean(!discovery %in% truematch), mean(truematch %in% discovery))



#### power of sequest higher 
fdrpow_qval_ls = readRDS('~/Dropbox/clipper/proteomics/realdata_proteomics_sequest_qvalue_fdppow_ls.rds')
fdrpow_clipper_ls = readRDS('~/Dropbox/clipper/proteomics/realdata_proteomics_sequest_clipper_fdppow_ls.rds')
fdrpow_qval_ls
fdrpow_clipper_ls
