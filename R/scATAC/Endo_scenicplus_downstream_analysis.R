# Endothelial network with scenicplus
################################ prepared data
# RNA
table(Endo_scRNA1$annotation6)
# Apold1+ Lsec     EC_other 
#                1939         5207 
library(SeuratDisk)
Idents(Endo_scRNA1) <- Endo_scRNA1$annotation6
# we used the SCT data
Endo_scRNA1@assays$RNA <- Endo_scRNA1@assays$SCT
SaveH5Seurat(Endo_scRNA1, filename = "Endo_scRNA1.h5Seurat")
Convert("Endo_scRNA1.h5Seurat", dest = "h5ad")
# Endo_scRNA1.h5ad

# ATAC
table(Endo_scATAC1$annotation6)
# Apold1+ Lsec     EC_other 
#                 193         4872 
Endo_scATAC1_meta <- Endo_scATAC1@meta.data[, c("annotation6", "seurat_clusters", "Group")]
Endo_scATAC1_meta$barcode <- rownames(Endo_scATAC1_meta)
Endo_scATAC1_meta <- Endo_scATAC1_meta[, c("barcode", "annotation6", "seurat_clusters", "Group")]
Endo_scATAC1_meta$barcode <- sub(".*_", "", Endo_scATAC1_meta$barcode)
Endo_scATAC1_meta <- Endo_scATAC1_meta[!duplicated(Endo_scATAC1_meta$barcode),]
write.table(Endo_scATAC1_meta, file = "cell_data.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


############################################# downstream analysis and plot
############################################# 内皮细胞调控网络下游分析
# eRegulon_metadata_result3：内皮细胞scenic+调控网络结果, 获取的方式在脚本Endo_get_network_table.ipynb
setwd("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/out_4/")
library(readr)
eRegulon_metadata_result3 <- read_csv("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/out_4/eRegulon_metadata_result3.csv")

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


############################################# CoveragePlot of Hbegf
consensus_regions <- read.delim("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/consensus_peak_calling/consensus_regions.bed", header=FALSE)
library(regioneR)
gr1 <- toGRanges(consensus_regions)
# Clusters
DefaultAssay(Endo_scATAC1) <- "ATAC"
pdf("Hbegf_cluster_CoveragePlot.pdf", width = 9, height = 4)
CoveragePlot(
  object = Endo_scATAC1,
  region = "Hbegf", group.by = "annotation6",
  ranges = gr1,
  features = "Hbegf",
  extend.downstream = 5000, extend.upstream = 2000,
  ranges.title = "MACS2"
)
dev.off()
# Group
pdf("Hbegf_Group_CoveragePlot.pdf", width = 9, height = 4)
CoveragePlot(
  object = Endo_scATAC1,
  region = "Hbegf", group.by = "Group",
  ranges = gr1,
  features = "Hbegf",
  extend.downstream = 5000, extend.upstream = 2000,
  ranges.title = "MACS2"
)
dev.off()


############################################# 扫描Hbegf启动子区域的模体
# Hbegf 基因的坐标：chr18-36504927-36515805
library(motifmatchr)
library(JASPAR2020) #JASPAR2020
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10) #BSgenome.Hsapiens.UCSC.hg19
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
peaks2 <- GRanges(seqnames = c("chr18"),
                  ranges = IRanges(start = c(36518891),end = c(36519391)), strand = "-")
# Get motif matches for example motifs in peaks
motif_ix <- matchMotifs(pfm, peaks2, genome = "mm10") 
motif_matches_matrix <- motifMatches(motif_ix) # Extract matches matrix from SummarizedExperiment result
mm <- as.matrix(motif_matches_matrix)
# motif的基因名字
pfm_motif_name <- Endo_scATAC1@assays[["ATAC"]]@motifs@motif.names
pfm_motif_name <- as.data.frame(pfm_motif_name)
pfm_motif_name <- t(pfm_motif_name)
pfm_motif_name <- as.data.frame(pfm_motif_name)
mm <- t(mm)
mm <- as.data.frame(mm)
MM <- subset(mm, V1 == "TRUE")
MM$motif_name <- pfm_motif_name$V1[match(rownames(MM), rownames(pfm_motif_name))]
head(MM)
write.table(MM, "Hbegf_peaks_36518891_36519391_motifmatchr_find_motif_result.txt")

motif_gene_n <- "FOSL1"
pdf(paste(motif_gene_n,"motif.pdf", sep = "_"),width = 3.5, height = 3)
MotifPlot(
  object = Endo_scATAC1,
  motifs = "MA0477.2"
)
dev.off()

motif_gene_n <- "FOSL2"
pdf(paste(motif_gene_n,"motif.pdf", sep = "_"),width = 3.5, height = 3)
MotifPlot(
  object = Endo_scATAC1,
  motifs = "MA0478.1"
)
dev.off()

motif_gene_n <- "JUN"
pdf(paste(motif_gene_n,"motif.pdf", sep = "_"),width = 4.5, height = 3)
MotifPlot(
  object = Endo_scATAC1,
  motifs = "MA0489.1"
)
dev.off()

### Footprint
Endo_scATAC1 <- Footprint(
  object = Endo_scATAC1, in.peaks = T, 
  motif.name = c("FOSL1", "FOSL2", "JUN"), 
  genome = BSgenome.Mmusculus.UCSC.mm10
)
Endo_scATAC1$Two_group <- factor(Endo_scATAC1$Two_group, levels = c("Apold1+ Lsec(Stress)", "Other"))
Idents(Endo_scATAC1) <- Endo_scATAC1$Two_group
PlotFootprint(Endo_scATAC1, features = c("FOSL1"), label = F)
PlotFootprint(Endo_scATAC1, features = c("FOSL2"), label = F)
PlotFootprint(Endo_scATAC1, features = c("JUN"), label = F)
# subset, Group
Apold1_Lsec_obj <- Footprint(
  object = subset(Endo_scATAC1, annotation6 == "Apold1_Lsec"), in.peaks = T, 
  motif.name = c("FOSL1", "FOSL2", "JUN"), 
  genome = BSgenome.Mmusculus.UCSC.mm10
)
PlotFootprint(Apold1_Lsec_obj, features = c("FOSL1"), label = F, group.by =  "Group")
PlotFootprint(Apold1_Lsec_obj, features = c("FOSL2"), label = F, group.by = "Group")
PlotFootprint(Apold1_Lsec_obj, features = c("JUN"), label = F, group.by = "Group")

### split vlnplot
VlnPlot(Endo_scRNA1, features = c("Fosl1"), group.by = "Clusters", split.by = "Group") + theme(axis.title.x = element_blank()) & geom_point(stat = 'summary',fun=mean, position = position_dodge(width = 0.9))
VlnPlot(Endo_scRNA1, features = c("Fosl2"), group.by = "Clusters", split.by = "Group") + theme(axis.title.x = element_blank()) & geom_point(stat = 'summary',fun=mean, position = position_dodge(width = 0.9))
VlnPlot(Endo_scRNA1, features = c("Jun"), group.by = "Clusters", split.by = "Group") + theme(axis.title.x = element_blank()) & geom_point(stat = 'summary',fun=mean, position = position_dodge(width = 0.9))


############################################# 计算相关性
avg1 <- AverageExpression(subset(Endo_scRNA1, Clusters == "Apold1+ Lsec(Stress)"), features = c("Hbegf","Fosl1", "Fosl2", "Jun"), assays = "SCT", group.by = "Group")$SCT
avg1_df <- as.data.frame(t(avg1))
cor(avg1_df)
#           Hbegf     Fosl1     Fosl2       Jun
# Hbegf 1.0000000 0.6843026 0.8487261 0.2646273
# Fosl1 0.6843026 1.0000000 0.7240405 0.0625012
# Fosl2 0.8487261 0.7240405 1.0000000 0.5318654
# Jun   0.2646273 0.0625012 0.5318654 1.0000000

# END！