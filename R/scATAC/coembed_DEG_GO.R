# coembed对象
#### 细胞类型差异表达基因和GO富集分析
Idents(coembed) <- coembed$celltype
all_markers <- FindAllMarkers(coembed, only.pos = F, logfc.threshold = 0.1, min.pct = 0.1)
all_markers_sig <- subset(all_markers, p_val_adj < 0.01)
all_markers_sig_fc <- subset(all_markers_sig, avg_log2FC > 1)

library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
topn.markers = 150
cluster_markers <- all_markers_sig_fc %>% group_by(cluster) %>% top_n(topn.markers, wt = avg_log2FC) #tibble
clusters <- unique(cluster_markers$cluster)
enrichment_results <- list()
for (cluster in clusters) {
  cluster_genes_symbol <- cluster_markers[cluster_markers$cluster == cluster, "gene"]
  print(head(cluster_genes_symbol))
  cluster_genes_entrez <- mapIds(org.Mm.eg.db, as.character(cluster_genes_symbol$gene), column = "ENTREZID", keytype = "SYMBOL")
  cluster_genes_entrez <- cluster_genes_entrez[!is.na(cluster_genes_entrez)]
  if (length(cluster_genes_entrez) == 0) {
    next
  }
  res1 <- enrichGO(gene = cluster_genes_entrez,
                   OrgDb = org.Mm.eg.db, 
                   readable = TRUE,
                   ont = "BP",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01)
  enrichment_results[[as.character(cluster)]] <- res1
}
result_df <- data.frame(cluster = character(),
                        description = character(),
                        pvalue = numeric(),
                        qvalue = numeric(),
                        stringsAsFactors = FALSE)
for (i in 1:length(enrichment_results)) {
  cluster <- clusters[i]
  enr <- enrichment_results[[i]]
  
  cluster_result <- enr@result
  cluster_result$cluster <- cluster
  
  result_df <- rbind(result_df, cluster_result)
}
show_terms <- c("GO:0044282", "GO:0006094",
                "GO:0001935", "GO:0048010",
                "GO:0006909", "GO:0071674",
                "GO:0046631", "GO:0007159", 
                "GO:0042113", "GO:0002377",
                "GO:0030593", "GO:0071347",
                "GO:0042267", "GO:0002223",
                "GO:0072175", "GO:0006814",
                "GO:0030198", "GO:0030199"
)
show_terms_vector <- c()
for (i in show_terms) {
  id <- grep(i, result_df$ID)
  show_terms_vector <- c(show_terms_vector, id)
}
show_terms_df <- result_df[show_terms_vector, ]
show_terms_df$cluster <- factor(show_terms_df$cluster, levels = levels(coembed$celltype))
show_terms_df$Description <- factor(show_terms_df$Description, levels = unique(result_df[show_terms_vector, ]$Description))

gp <- ggplot(show_terms_df, aes(Count, Description)) + 
  geom_point(aes(size = Count, color = -1 * log10(p.adjust))) +
  viridis::scale_color_viridis(option = "C")+
  labs(color = expression(-log[10](p.adjust)), size = "Count",
       x = "", y = "", title = "") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

# END！