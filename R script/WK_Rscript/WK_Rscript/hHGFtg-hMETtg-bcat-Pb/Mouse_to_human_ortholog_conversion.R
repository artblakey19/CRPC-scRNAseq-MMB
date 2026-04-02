#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named Symbol
your_rank_file <- rename(your_rank_file, "symbol" = "gene")
genelist = your_rank_file$symbol
ortho_your_rank_file <- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

your_rank_file <- inner_join(your_rank_file,ortho_your_rank_file,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(your_rank_file,"converted_gsea_rank.txt")