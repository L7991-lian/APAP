# Endothelial:clustering, annotation
# All为所有细胞类型的scRNA对象

# Endothelial clustering
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(ggplot2)
library(tidyverse)
library(harmony)
library(readxl)

sce.normaize <- function(sce, variable.features.n = 2000, vars.to.regress = c("nCount_RNA"), return.only.var.genes = FALSE, PC = 20, resolution = 0.5){
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  sce <- SCTransform(sce, variable.features.n = variable.features.n, vars.to.regress = vars.to.regress, return.only.var.genes = return.only.var.genes)
  varigenes <- VariableFeatures(object = sce)
  # varigenes.d.mt <- grep("^mt-", varigenes, value = T, invert = T)
  sce <- RunPCA(sce, features = varigenes, npcs=60) 
  sce <- RunUMAP(sce, dims = 1:PC, assay = "SCT")
  sce <- FindNeighbors(sce, dims = 1:PC)
  sce <- FindClusters(sce, resolution = resolution) %>% identity()
  DimPlot(sce, reduction = "umap", label = T) + DimPlot(sce, reduction = "umap", label = F, group.by = "Group")
  ggsave("PCA_UMAP.pdf", p1, width = 14, height = 8)
  return(sce)
}
# Dims clustering 
harmony.clustering <- function(sce, Dims = seq(10, 20, 5), resolution = 0.05){
  # create dir
  Dims <- Dims
  namefc <- c()
  for (l in Dims) {
    namef <- paste("Dim", l, sep = "_")
    dir.create(namef)
    namefc <- c(namefc, namef)
  }
  # clustering
  pathwd <- getwd()
  for (m in 1:length(Dims)) {
    path <- paste(pathwd, namefc[m], sep = "/")
    setwd(path)
    sce <- sce %>%
      RunUMAP(reduction = "harmony", dims = 1:Dims[m]) %>%
      FindNeighbors(reduction = "harmony", dims = 1:Dims[m]) %>%
      FindClusters(reduction.type = "harmony", resolution = resolution) %>% identity()
    DimPlot(sce, reduction = "umap", label = T)+ggtitle(Dims[m])|DimPlot(sce, reduction = "umap", group.by = "Group", label = F)
    return(sce)
  }
}
## RUN
# Endothelial
sce <- subset(All, celltype == "Endothelial")
sce <- sce.normaize(sce, variable.features.n = 3000, vars.to.regress = NULL, return.only.var.genes = FALSE, PC = 20, resolution = 0.1)
# Harmony
library(harmony)
sce <- sce %>%
  RunHarmony("Group", assay.use="SCT", verbose = FALSE)
Endo_scrna_harmony_sce <- harmony.clustering(sce, Dims = c(22), resolution = 0.05)


###### annotation
Endo_scRNA1$Clusters <- Endo_scRNA1$seurat_clusters
Endo_scRNA1$Clusters <- as.numeric(Endo_scRNA1$Clusters)
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "1")] <- "Chst2+ Lsec"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "2")] <- "Apold1+ Lsec(Stress)"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "3")] <- "Alb+ Lsec"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "4")] <- "Rpl32+ EC"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "5")] <- "Apold1+ Lsec(Stress)"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "6")] <- "Rspo3+ Lvec(Pericentral)"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "7")] <- "Efnb1+ Lvec(Periportal)"
Endo_scRNA1$Clusters[which(Endo_scRNA1$Clusters == "8")] <- "Mki67+ EC(Proliferating)"
Endo_scRNA1$Clusters <- factor(Endo_scRNA1$Clusters, levels = c("Rspo3+ Lvec(Pericentral)", "Efnb1+ Lvec(Periportal)", "Chst2+ Lsec", "Apold1+ Lsec(Stress)", "Alb+ Lsec", "Rpl32+ EC", "Mki67+ EC(Proliferating)"))
DimPlot(Endo_scRNA1, group.by = "Clusters", cols = cols, pt.size = 0.8)


###### other 
# jjDotPlot
Endo_scRNA1_clusters_marker_genes <- c("Vwf", "Stab2", "Rspo3", "Efnb1", "Chst2", "Apold1", "Alb", "Rpl32", "Mki67")
jjDotPlot(object = Endo_scRNA1, col.min = 0,base_size = 14, textSize = 14, bar.width = 3, segWidth = 0.3,
               gene = Endo_scRNA1_clusters_marker_genes, xtree = F, ytree = F, legend.position = "bottom",
               id = 'Clusters')

# Roe_Analysis
Roe_Analysis <- function(train, tissue = "tissue", celltype = "finalannotate", plot_location = ".", data_location = ".", perfix = 'rna') {
  library(epitools)
  library(reshape2)
  library(tidyverse)
  data <- table(train@meta.data[, tissue], train@meta.data[, celltype]) %>%
    matrix(nrow = dim(table(train@meta.data[, tissue])))
  ei_data <- expected(data)
  data_roe <- data / ei_data
  colnames(data_roe) <- names(table(train@meta.data[, celltype]))
  rownames(data_roe) <- names(table(train@meta.data[, tissue]))
  data_roe_plot <- melt(data_roe)
  data_roe_plot$significant <- case_when(
    data_roe_plot$value < 1 ~ "±",
    data_roe_plot$value > 1 & data_roe_plot$value < 1.5 ~ "+",
    data_roe_plot$value > 1.5 & data_roe_plot$value < 2 ~ "++",
    data_roe_plot$value > 2 ~ "+++",
  )
  colnames(data_roe_plot)[colnames(data_roe_plot) == "Var1"] <- "tissue"
  colnames(data_roe_plot)[colnames(data_roe_plot) == "Var2"] <- "celltype"
  p <- ggplot(data_roe_plot, aes(x = tissue, y = celltype, label = significant)) +
    geom_tile(aes(fill = value), color="gray",size=0.5) +
    coord_equal() + 
    geom_text(color = "black") + 
    scale_fill_viridis_c(option = "C", direction = 1, alpha = 1) +
    theme(axis.title = element_blank(), text = element_text(size = 14, colour = "black"), legend.title = element_blank())
  pdf(paste0(plot_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".pdf"), width = 5, height = 6)
  print(p)
  dev.off()
  write.table(data_roe, paste0(data_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".txt"), sep = "\t", quote = FALSE)
}
Roe_Analysis(Endo_scRNA1, tissue = "Group", celltype = "Clusters", perfix = 'Endo_scRNA')

# Stress signatures
stress_genes_ec_scrna <- c("Atf3", "Atf4", "Atf6","Egr1", "Fosl1","Hif1a", "Hsp90aa1", "Hspa1a", "Hspa1b", "Hspa8","Hspa9", "Hspb1", "Hspd1", "Hsph1", "Bcl10", "Bcl2l1", "Cdkn1a", "Hmox1", "Tlr4")
jjDotPlot(object = Endo_scRNA1, col.min = 0,
               base_size = 14, textSize = 16, 
               bar.width = 3, segWidth = 0.3,
               gene = stress_genes_ec_scrna, 
               xtree = F, ytree = F, legend.position = "bottom",
               id = 'Clusters')+ ggtitle("Stress signatures")

###### Egfr ligands
# list_cellchat_object[[2]]为6小时那个样本的cellchat通讯分析结果
regeneration_proliferation <- c("Cd81","Egfr","Lrp1","Ptprf")
regeneration_proliferation_df <- findEnrichedSignaling(list_cellchat_object[[2]], features = regeneration_proliferation, pattern ="incoming", thresh = 0.05)
regeneration_proliferation_df_endo_to_hep <- subset(regeneration_proliferation_df, target == "Hepatocyte" & source == "Endothelial")
regeneration_proliferation_df_endo_to_hep_signalling <- regeneration_proliferation_df_endo_to_hep[, c(3:5)]

Egfr_signaling <- subset(regeneration_proliferation_df_endo_to_hep_signalling, receptor == "Egfr")
# ligands_of_Egfr <- c("Hsp90aa1","Tgfb1", "Grn","Hbegf", "Tgfa", "Itgb1", "Icam1", "Epha2")
# Egfr_signaling <- ligands_of_Egfr
jjDotPlot(object = Endo_scRNA1, col.min = 0,base_size = 14, textSize = 14, 
               gene = unique(Egfr_signaling$ligand), 
               xtree = T, ytree = F, legend.position = "bottom", bar.width = 3, segWidth = 0.3,
               id = 'Clusters') + ggtitle("Egfr ligands")


###### DEG and GO enrichment analysis
DefaultAssay(Endo_scRNA1) <- "SCT"
Idents(Endo_scRNA1) <- Endo_scRNA1$Clusters
object.name = "Endo_scRNA_Clusters"
min.pct = 0.25
logfc.threshold = 0.25
diff.wilcox = FindAllMarkers(Endo_scRNA1, only.pos = T, min.pct = min.pct, logfc.threshold = logfc.threshold)
all.markers = diff.wilcox %>% filter(p_val_adj <= 0.01) %>% group_by(cluster)
all.markers <- all.markers[-grep("mt-", all.markers$gene) , ]

library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
topn.markers = 100
cluster_markers <- all.markers %>% group_by(cluster) %>% top_n(topn.markers, wt = avg_log2FC)
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
plot.ncol = 3
patchwork::wrap_plots(segment_plot_list, ncol = plot.ncol, )

# END!