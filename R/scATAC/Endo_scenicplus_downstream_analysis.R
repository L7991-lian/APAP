# 内皮细胞调控网络下游分析
# eRegulon_metadata_result3：内皮细胞scenic+调控网络结构

### TF表达Dotplot
Endo_TFs <- unique(eRegulon_metadata_result3$TF)
Endo_scRNA1$Two_group <- as.character(Endo_scRNA1$Two_group)
Endo_scRNA1$Two_group[which(Endo_scRNA1$Two_group == "Cdkn1a+ Lsec")] <- "Apold1+ Lsec(Stress)"
Endo_scRNA1$Two_group[which(Endo_scRNA1$Two_group == "Other")] <- "Other"
Endo_scRNA1$Two_group <- factor(Endo_scRNA1$Two_group, levels = c("Apold1+ Lsec(Stress)", "Other"))

library(scRNAtoolVis)
# scRNA
id_obj = 'Two_group'
p_tfscrna <- jjDotPlot(object = Endo_scRNA1,
                       gene = rev(Endo_TFs), 
                       col.min = 0, 
                       base_size = 14, textSize = 14, 
                       bar.width = 2,segWidth = 0.3, 
                       x.text.angle = 45,  
                       x.text.hjust = 1,x.text.vjust = 1,
                       xtree = F, ytree = F, 
                       legend.position = "right",
                       id = id_obj) + coord_flip()
# 排序TF基因
a=p_tfscrna$data
a$id <- as.character(a$id)
gene_df <- data.frame("gene" = unique(a$gene))
seq_c <- seq(1, nrow(a),2) #2是两次重复，两个组
for (i in seq_along(seq_c)) {
  value_avg_exp <- a$avg.exp[seq_c[i]]-a$avg.exp[c(seq_c[i]+1)]
  gene_df$cluster[i] <- ifelse(value_avg_exp > 0, a$id[seq_c[i]], a$id[c(seq_c[i]+1)])
}
# cluster select
id_v <- unique(gene_df$cluster)
id_1_df <- subset(a, id %in% c(id_v[1]))
id_2_df <- subset(a, id %in% c(id_v[2]))
# gene select
id_1_df <- subset(id_1_df, gene %in% subset(gene_df, cluster == id_v[1])$gene)
id_2_df <- subset(id_2_df, gene %in% subset(gene_df, cluster == id_v[2])$gene)
# order
id_1_df <- id_1_df[order(id_1_df$pct.exp,decreasing = T),]
id_2_df <- id_2_df[order(id_2_df$pct.exp,decreasing = T),]
Lsec_Hbegf_specific_tfs_df <- id_1_df
# 排序后绘图
id_obj = 'Two_group'
p_tfscrna <- jjDotPlot(object = Endo_scRNA1,
                       gene = c(id_1_df$gene,id_2_df$gene), 
                       col.min = 0, 
                       base_size = 14, textSize = 14, 
                       segWidth = 0.3, bar.width = 2.5,
                       plot.margin = c(0,0.5,0,0),
                       x.text.angle = 45,  
                       x.text.hjust = 1,x.text.vjust = 1,
                       xtree = F, ytree = F, 
                       legend.position = "top",
                       id = id_obj)

### Apold1+Lsec的核心调控网络
Apold1_clusters_tfs <- id_2_df$gene
id_vector2 <- c()
for (i in Apold1_clusters_tfs) {
  id <- grep(i, eRegulon_metadata_result3$TF)
  id_vector2 <- c(id_vector2, id)
}
Apold1_clusters_tfs_eregulon_metadata <- eRegulon_metadata_result3[id_vector2,]
# 保留target为差异表达基因
Endo_scRNA_clusters_markers  <- read.table("Endo_scRNA_Clusters_FindAllMarkers.txt")
Apold1_clusters_markers <- subset(Endo_scRNA_clusters_markers, cluster == "Lsec_2") #Apold1+ Lsec的差异基因
Apold1_clusters_markers$pct_d <- Apold1_clusters_markers$pct.1 - Apold1_clusters_markers$pct.2
# Apold1_clusters_markers$pct_d
# 0%    25%    50%    75%   100% 
# -0.050  0.076  0.098  0.124  0.291 
Apold1_clusters_markers <- Apold1_clusters_markers[Apold1_clusters_markers$pct_d > 0.124,]
Apold1_clusters_markers <- Apold1_clusters_markers[Apold1_clusters_markers$avg_log2FC > 1, ]
# Apold1+ Lsec亚群在时间上的差异表达基因
Endo_scRNA_Apold1_obj <- subset(Endo_scRNA1, Clusters == "Apold1+ Lsec(Stress)")
# 96h只有4个细胞
Endo_scRNA_Apold1_obj <- subset(Endo_scRNA_Apold1_obj, Group != "96h")
Endo_scRNA_Apold1_obj$Group <- droplevels(Endo_scRNA_Apold1_obj$Group)
Idents(Endo_scRNA_Apold1_obj) <- Endo_scRNA_Apold1_obj$Group
Endo_scRNA_Apold1_allmarkers <- FindAllMarkers(Endo_scRNA_Apold1_obj, only.pos = T, min.pct = 0.75,logfc.threshold = 1)
Endo_scRNA_Apold1_allmarkers_sig <- subset(Endo_scRNA_Apold1_allmarkers, p_val_adj < 0.01)
Endo_scRNA_Apold1_allmarkers_sig_delete0h <- subset(Endo_scRNA_Apold1_allmarkers_sig, cluster != "0h")
a=intersect(unique(Endo_scRNA_Apold1_allmarkers_sig_delete0h$gene),intersect(unique(Apold1_clusters_markers$gene), unique(Apold1_clusters_tfs_eregulon_metadata$Gene)))
id_vector <- c()
for (i in a) {
  id <- grep(i, Apold1_clusters_tfs_eregulon_metadata$Gene)
  id_vector <- c(id_vector, id)
}
Apold1_clusters_tfs_eregulon_metadata_DEG_target <- Apold1_clusters_tfs_eregulon_metadata[id_vector,]
Apold1_clusters_tfs_eregulon_metadata_DEG_target <- Apold1_clusters_tfs_eregulon_metadata_DEG_target[-which(Apold1_clusters_tfs_eregulon_metadata_DEG_target$TF == Apold1_clusters_tfs_eregulon_metadata_DEG_target$Gene),]
write.table(Apold1_clusters_tfs_eregulon_metadata_DEG_target, "Apold1_clusters_tfs_eregulon_metadata_DEG_target.txt")
# TF-GENE network
# net
net_tf_GENE <- Apold1_clusters_tfs_eregulon_metadata_DEG_target[, c("TF", "Gene", "TF2G_importance", "TF2G_regulation")]
colnames(net_tf_GENE) <- c("Source", "Target", "TF2G_importance", "Regulation")
net_tf_GENE <- net_tf_GENE[!duplicated(net_tf_GENE),]
net_tf_GENE$Regulation <- ifelse(net_tf_GENE$Regulation > 0, "positive", "negative")
#node
net_tf_GENE_info <- data.frame(name = c(unique(net_tf_GENE$Target), unique(net_tf_GENE$Source)))
net_tf_GENE_info$node_type <- "Gene"
net_tf_GENE_info <- net_tf_GENE_info[!duplicated(net_tf_GENE_info),]
net_tf_GENE_info$node_type[match(Apold1_clusters_tfs,net_tf_GENE_info$name)] <- "TF"
# save
write.csv(net_tf_GENE, "Apold1_clusters_tfs_eregulon_metadata_DEG_target_network.csv", quote = F)
write.csv(net_tf_GENE_info, "Apold1_clusters_tfs_eregulon_metadata_DEG_target_node.csv", quote = F)
write.csv(net_tf_GENE[,c("Source", "Target", "Regulation", "TF2G_importance")], "Apold1_clusters_tfs_eregulon_metadata_DEG_target_edge.csv", quote = F)

###










