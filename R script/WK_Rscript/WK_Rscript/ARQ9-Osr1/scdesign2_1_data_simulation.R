library(scDesign2)
library(parallel)
data_filename_list <- c('pbmc2_hca_drop-seq_sub.rds', 'pbmc2_hca_10x_sub.rds')
data_filename_list <- paste0("../dataset_sub/", data_filename_list)
# cell_type_sel <- c('B cell', 'CD14+ monocyte',
#                    'Natural killer cell', 'CD4+ T cell',
#                    'Cytotoxic T cell')
cell_type_sel_list <- matrix(rep(c('CD4+ T cell',
                                   'Cytotoxic T cell'), length(data_filename_list)),
                             byrow = FALSE, ncol = length(data_filename_list))
cell_type_sel_list
protocol_list <- c('drop-seq', '10x')
nde = c(200, 1000, 1500)
#protocol_list <- c('SMART-Seq2')
ncores <- detectCores()-4

for(iter1 in 1:length(data_filename_list)){
  data_filename <- data_filename_list[iter1]
  data_mat_original <- readRDS(data_filename)
  cell_type_name <- colnames(data_mat_original)
  
  # simulate 
  data_param <- lapply(1:nrow(cell_type_sel_list), function(iter2){
    cell_type_sel <- cell_type_sel_list[iter2, iter1]
    output_suffix <- paste0(protocol_list[iter1], '_', cell_type_sel)
    print(output_suffix)
    
    data_mat <- data_mat_original[, cell_type_name == cell_type_sel]
    cat('cell number =', ncol(data_mat), '\n')
    
    copula_model <- readRDS(file = paste0('copula_model_', output_suffix, '.rds'))
    
    list(copula_model = copula_model[[1]], n_cell = ncol(data_mat))
  })
  
  n_gene <- length(data_param[[1]]$copula$gene_sel1) + length(data_param[[1]]$copula$gene_sel2) +
    length(data_param[[1]]$copula$gene_sel3)
  
  param_mat <- matrix(0.0, nrow = n_gene, ncol = 6)
  param_mat[, c(2,5)] <- Inf
  param_mat[data_param[[1]]$copula_model$gene_sel1, 1:3] <-
    data_param[[1]]$copula_model$marginal_param1
  param_mat[data_param[[1]]$copula_model$gene_sel2, 1:3] <-
    data_param[[1]]$copula_model$marginal_param2
  param_mat[data_param[[2]]$copula_model$gene_sel1, 4:6] <-
    data_param[[2]]$copula_model$marginal_param1
  param_mat[data_param[[2]]$copula_model$gene_sel2, 4:6] <-
    data_param[[2]]$copula_model$marginal_param2
  
  abs_mean_diff <- abs(log(param_mat[, 3]+1) -log(param_mat[, 6]+1))
  abs_mean_diff_order <- order(abs_mean_diff, decreasing = TRUE)
  for (iter3 in 1:3){
    n_de <- nde[iter3]
    gene_idx_change <- abs_mean_diff_order[-(1:n_de)]
    gene_idx_remain <- abs_mean_diff_order[1:n_de]
    mat_temp <- t(apply(param_mat[gene_idx_change, ], 1, function(x){
      # print(x)
      x_new <- x
      m <- mean(c(x_new[3], x_new[6]))
      x_new[3] <- m
      x_new[6] <- m
      # print(x_new)
      x_new
    }))
    param_mat[gene_idx_change, ] <- mat_temp
    # sum(abs(param_mat[,3]-param_mat[,6]) >1e-5)
    # sum(abs(param_mat[gene_idx_change,3]-param_mat[gene_idx_change,6]) >1e-5)
    # sum(abs(param_mat[gene_idx_remain,3]-param_mat[gene_idx_remain,6]) >1e-5)
    
    data_param_new <- data_param
    data_param_new[[1]]$copula_model$gene_sel3 <- which(param_mat[, 3] < 1e-5)
    data_param_new[[1]]$copula_model$gene_sel2 <-
      (1:n_gene)[-c(data_param_new[[1]]$copula_model$gene_sel1,
                    data_param_new[[1]]$copula_model$gene_sel3)]
    data_param_new[[1]]$copula_model$marginal_param1 <-
      param_mat[data_param_new[[1]]$copula_model$gene_sel1, 1:3]
    data_param_new[[1]]$copula_model$marginal_param2 <-
      param_mat[data_param_new[[1]]$copula_model$gene_sel2, 1:3]
    
    data_param_new[[2]]$copula_model$gene_sel3 <- which(param_mat[, 6] < 1e-5)
    data_param_new[[2]]$copula_model$gene_sel2 <-
      (1:n_gene)[-c(data_param_new[[2]]$copula_model$gene_sel1,
                    data_param_new[[2]]$copula_model$gene_sel3)]
    data_param_new[[2]]$copula_model$marginal_param1 <-
      param_mat[data_param_new[[2]]$copula_model$gene_sel1, 4:6]
    data_param_new[[2]]$copula_model$marginal_param2 <-
      param_mat[data_param_new[[2]]$copula_model$gene_sel2, 4:6]
    
    save(gene_idx_change, gene_idx_remain, data_param_new,
         file = paste0('data_param_new_', protocol_list[iter1], '_', nde[iter3], '_de.rda'))
    for(iter2 in 1:nrow(cell_type_sel_list)){
      cell_type_sel <- cell_type_sel_list[iter2, iter1]
      output_suffix <- paste0(protocol_list[iter1], '_', cell_type_sel)
      print(output_suffix)
      
      set.seed(1)
      model_temp <- list(data_param_new[[iter2]]$copula_model)
      names(model_temp) <- cell_type_sel
      data_sim_copula <- mclapply(1:15, function(i){
        simulate_count_scDesign2(model_temp, data_param_new[[iter2]]$n_cell)
      }, mc.cores = ncores)
      saveRDS(data_sim_copula, file = paste0('data_sim_copula_', output_suffix, '_', nde[iter3], '_de.rds'))
      
      modified_model <- data_param_new[[iter2]]$copula_model
      modified_model$marginal_param1[, 1] <- 0.0
      modified_model$marginal_param2[, 1] <- 0.0
      set.seed(2)
      model_temp <- list(modified_model)
      names(model_temp) <- cell_type_sel
      data_sim_nozi <- mclapply(1:15, function(i){
        simulate_count_scDesign2(model_temp, data_param_new[[iter2]]$n_cell)
      })
      saveRDS(data_sim_nozi, file = paste0('data_sim_nozi_', output_suffix, '_', nde[iter3], '_de.rds'))
    }
    
  }
}






