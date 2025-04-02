# scATAC, Endothelial, clustering and annotation
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(future)
plan("multisession", workers = 4)

scATAC.high.combine <- readRDS("scATAC.data.high.rmDoublets.Rds")
## 我们只选取不在其他细胞类型分布的内皮细胞作后续的聚类分析
##filter false Endothelial cells
scATAC.high.combine.En.fiter =subset(scATAC.high.combine,seurat_clusters %in% c(0))
Endo_scATAC1 =subset(scATAC.high.combine.En.fiter, predicted.id == "Endothelial")

## Endothelial peak re-cluster 
DefaultAssay(Endo_scATAC1) <- 'ATAC'
Endo_scATAC1 <- RunTFIDF(Endo_scATAC1)
Endo_scATAC1 <- FindTopFeatures(Endo_scATAC1, min.cutoff = 'q80')
Endo_scATAC1 <- RunSVD(Endo_scATAC1, features = VariableFeatures(Endo_scATAC1))
Endo_scATAC1 <- RunUMAP(Endo_scATAC1, reduction = 'lsi', dims = 2:10)
Endo_scATAC1 <- FindNeighbors(Endo_scATAC1, reduction = 'lsi', dims = 2:10)
Endo_scATAC1 <- FindClusters(Endo_scATAC1, verbose = FALSE, algorithm = 3, resolution = 0.3)
DimPlot(Endo_scATAC1, label = TRUE) +ggtitle("Endothelial cells")|DimPlot(object = Endo_scATAC1, label = T, group.by = "Group")

## harmony
DefaultAssay(Endo_scATAC1) <- 'ATAC'
Endo_scATAC1 <- RunTFIDF(Endo_scATAC1)
Endo_scATAC1 <- FindTopFeatures(Endo_scATAC1, min.cutoff = 'q25')
Endo_scATAC1 <- RunSVD(Endo_scATAC1, features = VariableFeatures(Endo_scATAC1))
Endo_scATAC2 <- RunHarmony(
  object = Endo_scATAC1,
  group.by.vars = 'Group',
  assay.use = 'ATAC',reduction = "lsi", project.dim = FALSE, verbose = FALSE
)
Endo_scATAC2 <- RunUMAP(Endo_scATAC2, reduction = "harmony", dims = 2:10, verbose = FALSE)
Endo_scATAC2 <- FindNeighbors(object = Endo_scATAC2, reduction = 'harmony', dims = 2:10)
Endo_scATAC2 <- FindClusters(object = Endo_scATAC2, resolution = 0.2, verbose = FALSE, algorithm = 3)
DimPlot(Endo_scATAC2, label = TRUE)|DimPlot(object =Endo_scATAC2, label = T, group.by = "Group")


###
Endo_scATAC2$Clusters <- Endo_scATAC2$seurat_clusters
Endo_scATAC2$Clusters <- as.character(Endo_scATAC2$Clusters)
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "1")] <- "Chst15+ Lsec"
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "2")] <- "Mki67+ Lsec(Proliferating)"
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "3")] <- "Irgm1+ Lsec"
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "4")] <- "Rspo3+ Lvec(Pericentral)"
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "5")] <- "Efnb1+ Lvec(Periportal)"
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "6")] <- "Alb+ Lsec"
Endo_scATAC2$Clusters[which(Endo_scATAC2$seurat_clusters == "7")] <- "Apold1+ Lsec(Stress)"
Endo_scATAC2$Clusters <- factor(Endo_scATAC2$Clusters, levels = c("Chst15+ Lsec", "Mki67+ Lsec(Proliferating)", "Irgm1+ Lsec", "Rspo3+ Lvec(Pericentral)", "Efnb1+ Lvec(Periportal)", "Alb+ Lsec", "Apold1+ Lsec(Stress)"))

Endo_scATAC2@meta.data %>% 
  group_by(Group, Clusters) %>% 
  dplyr::summarise(count = n()) %>%
  ggplot() + 
  geom_bar(aes(Group, count, fill=Clusters), stat = "identity", position = "fill", width = 0.7) + 
  theme(axis.text.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 16, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain")) + 
  # scale_y_break(c(0.85,0.95),space=0,expand = c(0,0))+#对坐标轴进行截断，并设置范围12-16
  theme(legend.title=element_blank(), legend.text=element_text(size=14), axis.title.x = element_blank(), axis.title.y = element_blank())+
  theme(panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  labs(x="",y="Percentage") 

# Stress signatures
jjDotPlot(object = Endo_scATAC2, col.min = 0,
               base_size = 14, textSize = 14, 
               bar.width = 3, segWidth = 0.3,
               gene = stress_genes_ec_scrna, 
               xtree = F, ytree = F, legend.position = "bottom",
               id = 'Clusters')+ ggtitle("Stress signatures")

# Egfr ligand
ligands_of_Egfr <- c("Hsp90aa1","Tgfb1", "Grn","Hbegf", "Tgfa", "Itgb1", "Icam1", "Epha2")
jjDotPlot(object = Endo_scATAC2, col.min = 0,
               base_size = 14, textSize = 14, 
               bar.width = 3, segWidth = 0.3,
               gene = c(ligands_of_Egfr), 
               xtree = F, ytree = F, legend.position = "bottom",
               id = 'Clusters') + ggtitle("Egfr ligands")

#### ##### DEG and GO enrichment analysis
DefaultAssay(Endo_scATAC2) <- "ACTIVITY"
Idents(Endo_scATAC2) <- Endo_scATAC2$Clusters
object.name = "Endo_scATAC_Clusters"
min.pct = 0.25
logfc.threshold = 0.25
diff.wilcox = FindAllMarkers(Endo_scATAC2, only.pos = T, min.pct = min.pct, logfc.threshold = logfc.threshold)
all.markers = diff.wilcox %>% filter(p_val_adj <= 0.01) %>% group_by(cluster)

library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
topn.markers = 100
cluster_markers <- all.markers %>% group_by(cluster) %>% top_n(topn.markers, wt = avg_log2FC) #tibble
cluster_markers <- as.data.frame(cluster_markers)
clusters <- unique(cluster_markers$cluster)
enrichment_results <- list()
for (cluster in clusters) {
  cluster_genes_symbol <- cluster_markers[cluster_markers$cluster == cluster, "gene"]
  print(head(cluster_genes_symbol))
  cluster_genes_entrez <- mapIds(org.Mm.eg.db, as.character(cluster_genes_symbol), column = "ENTREZID", keytype = "SYMBOL")
  cluster_genes_entrez <- cluster_genes_entrez[!is.na(cluster_genes_entrez)]
  if (length(cluster_genes_entrez) == 0) {
    next
  }
  res1 <- enrichGO(gene = cluster_genes_entrez,
                   OrgDb = org.Mm.eg.db, 
                   readable = TRUE,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
  enrichment_results[[as.character(cluster)]] <- res1
}
segment_plot_list <- list()
for (i in 1:length(enrichment_results)) {
  res1 <- enrichment_results[[i]]
  cluster_name <- names(enrichment_results[i])
  print(paste("cluster:", cluster_name, sep = ""))
  topn <- 10
  top_10_descriptions <- head(res1[order(res1$p.adjust), "Description"], topn)
  
  library(ggplot2)
  library(forcats)
  p <- ggplot(res1[res1$Description %in% top_10_descriptions, ], 
              aes(-1*log10(p.adjust), 
                  fct_reorder(factor(Description), -1*log10(p.adjust)))) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = -1*log10(p.adjust), size = Count)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    ylab(NULL) +
    xlab(bquote("-log"[10] ~ "p.adjust")) +
    labs(title = paste0("Cluster_", cluster_name)) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    theme(text = element_text(size = 16), axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16))
  
  segment_plot_list[[cluster_name]] <- p
}
plot.ncol = 2
segment_plot <- patchwork::wrap_plots(segment_plot_list, ncol = plot.ncol, )

# END!