setwd("/Users/chenyiling/Dropbox/clipper/hi-C/rcode")
rep1 = read.table('/Users/chenyiling/Dropbox/clipper/hi-C/data/rawcount_primary.txt', header = F)
rep2 = read.table('/Users/chenyiling/Dropbox/clipper/hi-C/data/rawcount_replicate.txt', header = F)
rep1_features = convert2char(rep1[,1],rep1[,2])
rep2_features = convert2char(rep2[,1],rep2[,2])
rep2_sub = rep2[rep2_features %in% rep1_features, ]
rep1_sub = rep1[rep1_features %in% rep2_features, ]
all(rep2_sub[,1] == rep1_sub[,1])
all(rep2_sub[,2] == rep1_sub[,2])

rownames(rep1_sub ) = rep1_features[rep1_features %in% rep2_features]
rownames(rep2_sub) = rep2_features[rep2_features %in% rep1_features]
saveRDS(rep1_sub, file = 'primary_rep_processed.rds')
saveRDS(rep2_sub, file = 'replicate_rep_processed.rds')
