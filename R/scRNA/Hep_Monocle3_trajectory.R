# 肝细胞使用monocle3进行轨迹分析
# Hep_scRNA为scRNA-seq的肝细胞亚群聚类和注释后的结果

library(Seurat)
library(ggpubr)
library(Hmisc)
library(Seurat)
library(harmony)
library(ggpubr)
library(monocle3)
library(ggplot2)
seurat.obj=Hep_scRNA
##创建CDS对象并预处理数据
data <- GetAssayData(seurat.obj, assay = 'SCT', slot = 'counts')
cell_metadata <- seurat.obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 30)
cds <- align_cds(cds, alignment_group = "Group",useNames = TRUE )
cds <- reduce_dimension(cds, preprocess_method = "PCA")

#从seurat导入umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat.obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
#从seurat导入Clusters信息
cds@clusters$UMAP$clusters <- Hep_scRNA@meta.data[rownames(colData(cds)),"anno1"]
cds@clusters$UMAP$partitions <- factor(x = rep(1, length(rownames(colData(cds)))), levels = 1)
names(cds@clusters$UMAP$partitions) <- rownames(colData(cds))
#Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(cds)

#主图
tmp_euclidean_distance_ratio=5
tmp_minimal_branch_len=400
tmp_geodesic_distance_ratio=0.5
tmp_ncenter=100
cds <- learn_graph(cds,close_loop = F, use_partition = F,learn_graph_control=list("euclidean_distance_ratio"=tmp_euclidean_distance_ratio,
                                                                                  "minimal_branch_len"=tmp_minimal_branch_len,
                                                                                  "geodesic_distance_ratio"=tmp_geodesic_distance_ratio,
                                                                                  "ncenter"=tmp_ncenter))
#我们设定stress那群细胞中Graph的中间分支点那个细胞作为root
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                  label_leaves = FALSE, trajectory_graph_color = "purple", cell_size = 0.8,
                  label_branch_points = FALSE, 
                  label_roots = F) + 
  theme(legend.position="right", 
        # legend.key.size =  unit(0.2, "inches"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16), 
        text = element_text(size = 16), axis.text = element_text(size = 16)) + 
  # theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))  + 
  tidydr::theme_dr()+
  scale_color_gradientn(colours = c("#0088C3","#6FBCDB","#BCE5F1","#F2EFDA","#FFC482","#FF7B4C","#FF5831"))

#选取root起点到再生或者正常的那群细胞中一个节点的细胞（候选细胞），用作后续查找拟时变化基因和动态调控网络分析
cds_tmp <- choose_graph_segments(cds, clear_cds = F)
saveRDS(cds_tmp, "cds_choose_graph_segments.Rds")
cds_sub33 <- cds_tmp
id_s2 <- as.vector(sapply(rownames(colData(cds_sub33)), function(x) grep(x,rownames(colData(cds)))))
# 取出轨迹中的候选细胞
cds@colData$seurat_clusters <- Hep_scRNA$anno1
Trajectory2 <- as.character(cds@colData$seurat_clusters)
Trajectory2[setdiff(1:16017, id_s2)] <- "not_select"
cds@colData$Trajectory2 <- Trajectory2
cds@colData$Trajectory2 <- factor(cds@colData$Trajectory2, levels = c(levels(pseudotime_df1$annotation), "not_select"))
# 对应的候选细胞的seurat对象
sce_sub33 <- Hep_scRNA[, rownames(colData(cds_sub33))]
sce_sub33$pseudotime <- cds_sub33@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
saveRDS(sce_sub33, "Hep_choose_graph_segments.Rds")

pseudotime_df1 <- as.data.frame(sce_sub33$pseudotime)
colnames(pseudotime_df1) <- "pseudotime"
pseudotime_df1$annotation <- sce_sub33$anno1
pseudotime_df1$Group <- sce_sub33$Group
# 样本在拟时间的分布
ppsed3 <- ggplot(pseudotime_df1, aes(x = pseudotime, y = annotation_order, fill = annotation_order)) +#数据
  geom_density_ridges(stat = "density_ridges", from = 0) +#绘制山脊图
  theme(legend.position = "right", #图例去除
        axis.text.y = element_blank(),
        # axis.line.y = element_blank(),#去除Y轴轴线
        axis.ticks.y = element_blank()
  )+ #去除Y轴刻度
  # scale_y_discrete(expand = c(0, 0.5)) + #调整起始图形距离X轴距离
  scale_fill_manual(values = ggsci::pal_nejm()(6)) +
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"), text = element_text(size = 14, colour = "black"), legend.title = element_blank(), axis.text = element_text(size = 14, colour = "black")) #+
ggsave("geom_density_ridges_order_trajectory_Group_0.pdf", ppsed3, width = 6, height = 5)

#仅显示候选细胞标示对应cluster颜色的图
col_set2 <- c("red", "orange", "purple", "blue")
names(col_set2) <- levels(Hep_scRNA$anno1)
names(col_set2)[length(col_set2)] <- "not_select"
plot_cells(cds, color_cells_by = "Trajectory2", label_cell_groups = FALSE, 
                   label_leaves = FALSE, trajectory_graph_color = "purple", cell_size = 0.8,
                   label_branch_points = FALSE, group_cells_by = "Trajectory2",
                   label_roots = F) + 
  theme(legend.position="right", legend.key.size =  unit(0.2, "inches"), legend.title = element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16), axis.text = element_text(size = 16)) + 
  # theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))  + 
  tidydr::theme_dr() + scale_color_manual(values = col_set2)

#拟时轨迹基因
modulated_genes2 <- graph_test(cds_sub33, neighbor_graph = "principal_graph", cores = 30)
genes <- subset(modulated_genes2, q_value < 0.05 & morans_I > 0.005)
saveRDS(genes, "Filtered_Trajectory_genes.Rds")

#拟时轨迹基因的拟时聚类和表达heatmap
mat <- pre_pseudotime_matrix(cds_obj = cds_sub33, gene_list = genes$gene_short_name, assays = "normalized")
saveRDS(mat, "Filtered_Trajectory_genes_exp_mat.Rds")
mat <- na.omit(mat)
library(ClusterGVis)
ck <- clusterData(exp = mat,
                  cluster.method = "kmeans",
                  cluster.num = 5)
ck_gene_res = ck[["long.res"]]
gene.cluster <- ck_gene_res[, c("cluster", "gene")]
colnames(gene.cluster) <- c("group", "gene")
gene.cluster$group[which(gene.cluster$group == "1")] <- "C3"
gene.cluster$group[which(gene.cluster$group == "2")] <- "C2"
gene.cluster$group[which(gene.cluster$group == "3")] <- "C2"
gene.cluster$group[which(gene.cluster$group == "4")] <- "C1"
gene.cluster$group[which(gene.cluster$group == "5")] <- "C1"
genes_traj <- c(rev(subset(gene.cluster, group == "C1")$gene), subset(gene.cluster, group == "C2")$gene, subset(gene.cluster, group == "C3")$gene)
pt.matrix2 <- exprs(cds_sub33)[match(genes_traj, rownames(rowData(cds_sub33))),order(pseudotime(cds_sub33))]
pt.matrix2 <- t(apply(pt.matrix2,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix2 <- t(apply(pt.matrix2,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix2) <- genes_traj;
#heatmap
library(circlize)
cols <- c("#ED0000FF","#42B540FF" , "#00468BFF")
row_anno <- data.frame(Cluster = gene.cluster$group, row.names = gene.cluster$gene)
ann_colors = list(Cluster = setNames(cols, paste0("C", 1:3)))
pheatmap::pheatmap(pt.matrix2, annotation_legend = F,
                           cluster_cols = F, cluster_rows = F, show_rownames = F,
                           clustering_method = "ward.D", show_colnames = F,
                           color = colorRampPalette(c("blue4","#33FFFf","yellow", "red3"))(100),
                           annotation_row = row_anno, annotation_colors = ann_colors,
                           main = "Pseudotime")
# 拟时变化基因富集分析
library(clusterProfiler)
library(enrichplot)
group <- data.frame(gene=gene.cluster$gene,group=gene.cluster$group)
Gene_ID <- bitr(gene.cluster$gene, fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb="org.Mm.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO_tra <- compareCluster(
  ENTREZID~group,
  data=data,
  fun="enrichGO",
  OrgDb="org.Mm.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
data_GO_sim_tra <- clusterProfiler::simplify(data_GO_tra,
                                             cutoff=0.5,
                                             by="p.adjust",
                                             select_fun=min)
data_GO_sim_tra <- setReadable(data_GO_sim_tra, OrgDb = "org.Mm.eg.db")
saveRDS(data_GO_tra, "Cluster3_Trajectory_genes_Kmeans_GO_all.Rds")

