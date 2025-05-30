########## IReNA构建肝细胞拟时变化过程的调控网络
# 使用scATAC-seq(bams.bam)作为一个ATAC-seq的整体筛选转录因子的印记。
# cds_sub33即choose_graph_segments函数选择的肝细胞（候选肝细胞）
# Hep_scRNA即肝细胞的scRNA-seq的seurat对象
# mo@colData@listData[["seurat_clusters"]] 即Hep_scRNA$anno1，即肝细胞的亚群注释结果
# gene.cluster即肝细胞拟时变化基因的聚类信息

library(IReNA)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(monocle3)
monocle_object <- cds_sub33
seurat_object <- Hep_scRNA[, rownames(colData(cds_sub33))]
monocle_object <- monocle_object[rownames(mat), ]
seurat_object <- seurat_object[rownames(mat), ]

Monocle3.add_pseudotime <- function(seurat_object, monocle_object){
  se <- seurat_object
  mo <- monocle_object
  se[["Pseudotime"]] <- mo@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  se[["State"]] <- mo@colData@listData[["seurat_clusters"]]
  return(se)
}
seurat_with_time <- Monocle3.add_pseudotime(seurat_object, monocle_object)

expression_profile <- get_SmoothByBin_PseudotimeExp(seurat_with_time, Bin = 200, FcType = "Q95", method = "Pseudotime")
expression_profile_filter <- expression_profile
K.set = 5
library(gplots)
clustering <- clustering_Kmeans(expression_profile_filter, K1=K.set)
clustering.copy <- clustering
clustering$KmeansGroup[match(subset(gene.cluster, group == "C1")$gene, rownames(clustering))] <- "C1"
clustering$KmeansGroup[match(subset(gene.cluster, group == "C2")$gene, rownames(clustering))] <- "C2"
clustering$KmeansGroup[match(subset(gene.cluster, group == "C3")$gene, rownames(clustering))] <- "C3"
Kmeans_clustering_ENS <- add_ENSID(clustering, Spec1='Mm')
# infers the potential regulation for each gene in the expression profile through GENIE3
# 有G个基因，将会以每个基因作为输出（靶基因表达量），其余基因作为输入（调控基因表达量）来构建G个随机森林模型。进而计算每一个模型中，每个调控基因对靶基因的重要性，从而得出它们之间的调控关系。然后根据所有模型的调控关系进行排序。
library("GENIE3")
library("doParallel")
weightMat <- GENIE3(as.matrix(seurat_with_time@assays$SCT@data), nCores = 15, treeMethod = "RF",
                     K = "sqrt",
                     nTrees = 500)
weightMat <- getLinkList(weightMat2)
regulation <- weightMat[weightMat[,3] > as.numeric(quantile(weightMat$weight)[3]),] #the third column is the weight of the link

### add regulation type for each gene pair
regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS, regulation)

motif1 <- Tranfac201803_Mm_MotifTFsF
motifTF <- c()
for (i in 1:nrow(motif1)) {
  TF <- strsplit(motif1[i,5],';')[[1]]
  motifTF <- c(motifTF, TF)
}
regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]


### 用Hep scatac筛选footprint对应的调控网络，与上面GENIC3调控网络的交集作为最后的调控网络
######################################## （1）differential_peaks.bed
DefaultAssay(Hep_scATAC) <- 'ATAC'
Idents(Hep_scATAC) <- Hep_scATAC$Group
# wilcox is the default option for test.use
da_peaks_1 <- FindMarkers(
  object = Hep_scATAC,
  ident.1 = "12h",
  ident.2 = "0h",
  test.use = 'wilcox', only.pos = T,
  assay = "ATAC",
  min.pct = 0.1
)
da_peaks_1_sig <- da_peaks_1[da_peaks_1$p_val_adj < 0.01 & da_peaks_1$avg_log2FC > 0.1, ]
da_peaks_2 <- FindMarkers(
  object = Hep_scATAC,
  ident.1 = "96h",
  ident.2 = "0h",
  test.use = 'wilcox', only.pos = T,
  assay = "ATAC",
  min.pct = 0.1
)
da_peaks_2_sig <- da_peaks_2[da_peaks_2$p_val_adj < 0.01 & da_peaks_2$avg_log2FC > 0.1, ]
da_peaks_3 <- FindMarkers(
  object = Hep_scATAC,
  ident.1 = "96h",
  ident.2 = "12h",
  test.use = 'wilcox', only.pos = T,
  assay = "ATAC",
  min.pct = 0.1
)
da_peaks_3_sig <- da_peaks_3[da_peaks_3$p_val_adj < 0.01 & da_peaks_3$avg_log2FC > 0.1, ]

Hep_DEG_peaks <- c(rownames(da_peaks_1_sig), rownames(da_peaks_2_sig), rownames(da_peaks_3_sig))
Hep_DEG_peaks_df <- matrix(sapply(Hep_DEG_peaks, function(x) strsplit(x, split = "-")[[1]]), byrow = T, ncol = 3)
Hep_DEG_peaks_df <- as.data.frame(Hep_DEG_peaks_df)
Hep_DEG_peaks_df <- Hep_DEG_peaks_df[-(grep("GL", Hep_DEG_peaks_df$V1)),] ## 去掉特殊染色体的片段，否则报错
# JH584304.1
Hep_DEG_peaks_df <- Hep_DEG_peaks_df[-(grep("JH", Hep_DEG_peaks_df$V1)),] ## 去掉特殊染色体的片段，否则报错
write.table(Hep_DEG_peaks_df, 'differential_peaks.bed', quote = F, row.name = F, col.names = F, sep = '\t')
dim(Hep_DEG_peaks_df)
# [1] 4427    3


################################### （2）Calculate footprints: Linux
### Install RGT
# conda create -n RGT
# conda activate RGT
# pip install --user RGT
# rgt-hint

### 配置genome
# cd ~/rgtdata
# python setupGenomicData.py --mm10
### footprinting 
# rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-prefix=Hep bams.bam differential_peaks.bed
# 结果保存在Hep.bed


########################################## （3）Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)
# footprint occupancy score (FOS) 
motif1 <- Tranfac201803_Mm_MotifTFsF
# obtaining the motif sequences, then use fimo software to identify binding motifs in the footprints. 
###merge footprints whose distance is less than 4
filtered_footprints <- read.table('Hep.bed',sep = '\t')
filtered_footprints <- filtered_footprints[,c(1,2,3,5)] ### 必须，否则报错
fastadir <- '/home/lijinlian/rgtdata/mm10/'
# merged_fasta <- get_merged_fasta(filtered_footprints,fastadir) ## 报错，没有结果出来
# get_merged_fasta报错，我换用下面merge_footprints实现将距离接近（<4bp）的footprints合并

if (!require(data.table)) {
  install.packages("data.table")
}
library(data.table)
setDT(filtered_footprints)

# 定义合并足迹的函数
merge_footprints <- function(dt) {
  # 按照染色体和起始位置排序
  dt <- dt[order(V1, V2)]
  # 初始化合并后的足迹列表
  merged <- data.table(V1 = character(), V2 = integer(), V3 = integer(), V5 = integer())
  # 遍历足迹
  i <- 1
  while (i <= nrow(dt)) {
    # 取出当前足迹
    current <- dt[i, ]
    # 查找需要合并的足迹
    j <- i + 1
    while (j <= nrow(dt) && abs(dt[j, V2] - current$V3) < 4) {
      current$V3 <- max(current$V3, dt[j, V3])
      current$V5 <- max(current$V5, dt[j, V5])
      j <- j + 1
    }
    # 将合并后的足迹添加到结果列表中
    merged <- rbind(merged, current)
    # 更新下一个足迹的索引
    i <- j
  }
  return(merged)
}
merged_footprints <- merge_footprints(filtered_footprints)
head(merged_footprints)
write.table(merged_footprints, "Hep_merged_footprints.bed", sep = "\t", row.names = F, col.names = F, quote = F)
# 获取footprint的序列
# bedtools getfasta -name -fi /home/lijinlian/rgtdata/mm10/genome_mm10.fa -bed Hep_merged_footprints.bed -fo Hep_merged_footprint.fasta
# Hep_merged_footprint.fasta


### Identify differentially expressed genes related motifs
motif1 <- motifs_select(Tranfac201803_Mm_MotifTFsF, rownames(Kmeans_clustering_ENS)) ###Kmeans_clustering_ENS was obtained in part1
### run find_motifs()
fimodir <- '/home/lijinlian/software/meme-5.4.1/src/fimo'
outputdir1 = "/data2/lijinlian/APAP_project/Hep_scRNA/monocle3_wjl/IReNA_bulkatac/fimo/fimo/" # indicating the output path of shell script
outputdir <- "/data2/lijinlian/APAP_project/Hep_scRNA/monocle3_wjl/IReNA_bulkatac/fimo/fimo/"
motifdir <- '/data2/lijinlian/APAP_project/Irena_motif/'
sequencedir <- './Hep_merged_footprint.fasta'
find_motifs(motif1, step = 20, fimodir, outputdir1, outputdir, motifdir, sequencedir)
### run fimo_all script in shell
shell_code <- paste0('sh ', outputdir,'Fimo_All.sh')
system(shell_code, wait=TRUE)
# ### delete shell scripts
shell_code2 <- paste0('rm ',outputdir1,'Fimo*.sh')
system(shell_code2, wait=TRUE)

###Combine all footprints of motifs
combined <- combine_footprints(outputdir1)
peaks <- read.delim('differential_peaks.bed')
overlapped <- overlap_footprints_peaks(combined, peaks)
saveRDS(overlapped, "Hep_footprints_overlapped_motifs.Rds")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
list1 <- get_related_genes(overlapped, txdb = txdb, motif=Tranfac201803_Mm_MotifTFsF, Species = 'Mm')
###Get candidate genes/TFs-related peaks
list2 <- get_related_peaks(list1, Kmeans_clustering_ENS)
### output filtered footprints
write.table(list2[[1]],'filtered_footprints.bed', quote = F, row.names = F, col.names = F, sep = '\t')

## 报错的话，就去linux上跑这个命令
shell_code1 <- 'samtools view -hb -L filtered_footprints.bed bams.bam > Hep_1_filter.bam'
system(shell_code1, wait=TRUE)

# count the cuts of each position in footprints by wig_track(), and use these cuts to calculate the FOS of footprints to identify enriched TFs which determine the regulatory relationship
### calculate cuts of each position in footprints
bamfilepath1 <- 'Hep_1_filter.bam'
### set parameter 'workers' to make this function run in parallel
cuts1 <- cal_footprint_cuts(bamfilepath = bamfilepath1, bedfile = list2[[1]], workers = 40, index_bam = T)
cut_list <- list(cuts1)
### get related genes of footprints with high FOS
potential_regulation <- Footprints_FOS(cut_list, list2[[2]], FOS_threshold = 0.5)
### Use information of footprints with high FOS to refine regulatory relationships
filtered_regulatory <- filter_ATAC(potential_regulation, regulatory_relationships)
saveRDS(filtered_regulatory, "filtered_regulatory_fos0.5.Rds")

### 超几何检验
TFs_list <- network_analysis(filtered_regulatory, Kmeans_clustering_ENS, TFFDR1 = -log10(0.05), TFFDR2 = -log10(0.05))
saveRDS(TFs_list, "Hep_irena_bulkatac_TFs_list.Rds")
saveRDS(filtered_regulatory, "Hep_irena_filtered_regulatory.0.05.Rds")

# 富集转录因子间的调控网络
pdf("plot_tf_network.pdf", width = 7, height = 7)
plot_tf_network(TFs_list)
dev.off()
write.csv(TFs_list[["TF_network"]], "Hep_irena_bulkatac_enrichTF_regulate_module.csv", quote = F)

### TF pseudotime热图
motif1 <- Tranfac201803_Mm_MotifTFsF
filtered_regulatory_relationships <- filtered_regulatory[filtered_regulatory$Weight > 0.0001, ]
Module_Tfs <- unique(c(as.character(unique(filtered_regulatory_relationships$TFSymbol)), as.character(unique(filtered_regulatory_relationships$TargetSymbol[filtered_regulatory_relationships$TargetSymbol %in% motif1$TFs]))))
pt.matrix2 <- exprs(cds_sub33)[match(Module_Tfs, rownames(rowData(cds_sub33))),order(pseudotime(cds_sub33))]
pt.matrix2 <- t(apply(pt.matrix2,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix2 <- t(apply(pt.matrix2,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix2) <- Module_Tfs
library(circlize)
cols <- c("#ED0000FF","#42B540FF" , "#00468BFF")
row_anno <- data.frame(Cluster = Kmeans_clustering_ENS$KmeansGroup[match(Module_Tfs, Kmeans_clustering_ENS$Symbol)], row.names = Module_Tfs)
row_anno$Cluster <- paste0("C", row_anno$Cluster)
table(row_anno$Cluster)
# C1 C2 C3
# 79 62 28
ann_colors = list(Cluster = setNames(cols, paste0("C", 1:3)))
pheatmap::pheatmap(pt.matrix2[c(rownames(subset(row_anno, Cluster == "C1")), rownames(subset(row_anno, Cluster == "C2")),rownames(subset(row_anno, Cluster == "C3"))),], 
                           cluster_cols = F, cluster_rows = F, show_rownames = F,
                           clustering_method = "ward.D", 
                           color = colorRampPalette(c("navy","white", "red3"))(100),
                           annotation_row = row_anno, annotation_colors = ann_colors,
                           main = "Pseudotime")

# END！