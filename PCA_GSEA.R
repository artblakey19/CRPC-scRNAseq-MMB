# PCA로 뽑은 PC를 50위까지 GSEA하는 함수
# Gene Collection은 BP 사용

library(clusterProfiler)
library(org.Hs.eg.db)

run_pc_gsea <- function(seurat_obj, num_pcs = 50, ont = "BP", pvalueCutoff = 0.05) {
  pc_loadings <- Loadings(seurat_obj, reduction = "pca")
  results <- list()
  
  for (i in 1:num_pcs) {
    pc_name <- paste0("PC_", i)
    # PCA에서 얻은 Weights로 Ranking
    ranked_genes <- sort(pc_loadings[, pc_name], decreasing = TRUE)
    
    tryCatch({
      gsea_res <- gseGO(
        geneList = ranked_genes,
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = ont,
        pvalueCutoff = pvalueCutoff,
        verbose = FALSE
      )
      
      if (nrow(gsea_res@result) > 0) {
        top <- head(gsea_res@result, 3)
        results[[pc_name]] <- data.frame(
          PC = pc_name,
          Pathway = top$Description,
          NES = round(top$NES, 2),
          pvalue = signif(top$p.adjust, 3)
        )
        cat(pc_name, ": ", nrow(gsea_res@result), " significant pathways\n")
      } else {
        cat(pc_name, ": no significant pathways\n")
      }
    }, error = function(e) {
      cat(pc_name, ": error -", e$message, "\n")
    })
  }
  
  summary_df <- do.call(rbind, results)
  rownames(summary_df) <- NULL
  return(summary_df)
}

# 실행
# pc_gsea_summary <- run_pc_gsea(데이터 이름, num_pcs = 50)
# View(pc_gsea_summary)