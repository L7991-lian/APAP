# scATAC, Hepatocyte, clustering and annotation
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(future)
plan("multisession", workers = 4)

# clustering
scATAC.high.combine <- readRDS("scATAC.data.high.rmDoublets.Rds")
##filter Hepatocyte cells
## 我们只选取不在其他细胞类型分布的肝细胞作后续的聚类分析
scATAC.high.combine.He.fiter =subset(scATAC.high.combine, seurat_clusters %in% c(2,7,11,13))
Hep_scATAC =subset(scATAC.high.combine.He.fiter, predicted.id == "Hepatocyte")

## Hepatocyte peak re-cluster
DefaultAssay(Hep_scATAC) <- 'ATAC'
Hep_scATAC <- RunTFIDF(Hep_scATAC)
Hep_scATAC <- FindTopFeatures(Hep_scATAC, min.cutoff = 'q80')
Hep_scATAC <- RunSVD(Hep_scATAC, features = VariableFeatures(Hep_scATAC))
Hep_scATAC <- RunUMAP(Hep_scATAC, reduction = 'lsi', dims = 2:10)
Hep_scATAC <- FindNeighbors(Hep_scATAC, reduction = 'lsi', dims = 2:10)
Hep_scATAC <- FindClusters(Hep_scATAC, verbose = FALSE, algorithm = 3, resolution = 0.4)
DimPlot(Hep_scATAC, label = TRUE)|DimPlot(object = Hep_scATAC, label = T, group.by = "Group")


# annotation
Hep_scATAC$anno1 <- Hep_scATAC$seurat_clusters
Hep_scATAC$anno1 <- as.character(Hep_scATAC$anno1)
Hep_scATAC$anno1[which(Hep_scATAC$anno1 == "1")] <- "Cps1+"
Hep_scATAC$anno1[which(Hep_scATAC$anno1 == "2")] <- "Cyp2e1+"
Hep_scATAC$anno1[which(Hep_scATAC$anno1 == "3")] <- "Cdh1+"
Hep_scATAC$anno1[which(Hep_scATAC$anno1 == "4")] <- "Mt1+"
Hep_scATAC$anno1[which(Hep_scATAC$anno1 == "5")] <- "Kdr+"
Hep_scATAC$anno1[which(Hep_scATAC$anno1 == "6")] <- "Hsph1+"
Hep_scATAC$anno1 <- factor(Hep_scATAC$anno1, levels = c("Cps1+", "Cyp2e1+", "Cdh1+", "Mt1+", "Kdr+", "Hsph1+"))
DimPlot(Hep_scATAC, group.by = "anno1", cols = gcolors, pt.size = 1)|DimPlot(Hep_scATAC, group.by = "Group", cols = unname(gcolors), pt.size = 1) 

## jjDotPlot
Hep_scATAC_anno1_markers <- c("Cps1", "Igfbp2", "Mfsd2a", "Cyp2e1","Cdh1","Mt1", "Kdr", "Hsph1")
DefaultAssay(Hep_scATAC) <- "ACTIVITY"
jjDotPlot(object = Hep_scATAC, col.min = 0,base_size = 14, textSize = 14, segWidth = 0.4, bar.width = 3,
               gene = Hep_scATAC_anno1_markers, xtree = F, ytree = F, legend.position = "top", 
               id = 'anno1')

## histogram
Hep_scATAC@meta.data %>% 
  group_by(Group, anno1) %>% 
  dplyr::summarise(count = n()) %>%
  ggplot() + 
  geom_bar(aes(Group, count, fill=anno1), stat = "identity", position = "fill", width = 0.7) + 
  theme(axis.text.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 16, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = .5, face = "plain")) + 
  # scale_y_break(c(0.85,0.95),space=0,expand = c(0,0))+#对坐标轴进行截断，并设置范围12-16
  theme(legend.title=element_blank(), legend.text=element_text(size=16), axis.title.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values = gcolors) +
  labs(x="",y="Percentage") +
  # theme_bw() + 
  theme(panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  labs(x="",y="Percentage") 

# AddModuleScore: regeneration receptors
regeneration_receptors <- c("Cd81","Egfr","Lrp1","Ptprf")
Hep_scATAC <- AddModuleScore(Hep_scATAC, features = list(regeneration_receptors), name = "regeneration_receptors")
library(scales)
library(RColorBrewer)
FeaturePlot(Hep_scATAC, features = "regeneration_receptors1", pt.size = 1.2, order = T) &
  scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),
                        guide = "colorbar",
                        limits = c(min(Hep_scATAC$regeneration_receptors1), max(Hep_scATAC$regeneration_receptors1))) &

Nebulosa::plot_density(Hep_scATAC, 
                             features = c("Igfbp2"),
                             size = 1.2, 
                             # pal = "magma",
                             joint=TRUE, slot = "data",combine = T, adjust = 1,
                             reduction = "umap")
Nebulosa::plot_density(Hep_scATAC, 
                             features = c("Mfsd2a"),
                             size = 1.2, 
                             # pal = "magma",
                             joint=TRUE, slot = "data",combine = T, adjust = 1,
                             reduction = "umap")

