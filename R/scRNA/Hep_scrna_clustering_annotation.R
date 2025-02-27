#Hepatocyte: clustering and annotation
####
# variable.features.n = 2000, vars.to.regress = NULL, Dim13, harmony
# HVG2000 
# ## S3 method for class 'Seurat', DEFAULT
# RunHarmony(
#   object,
#   group.by.vars,
#   reduction.use = "pca",
#   dims.use = NULL,
#   reduction.save = "harmony",
#   project.dim = TRUE,
#   ...
# )
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
# Hepatocyte
sce <- subset(SeuratObject, celltype == "Hepatocyte")
sce <- sce.normaize(sce, variable.features.n = 2000, vars.to.regress = NULL, return.only.var.genes = FALSE, PC = 20, resolution = 0.05)
# Harmony
library(harmony)
sce <- sce %>%
  RunHarmony("Group", assay.use="SCT", verbose = FALSE)
Hep_scrna_harmony_sce <- harmony.clustering(sce, Dims = seq(13), resolution = 0.05)


#### annotation
Endo_scRNA1_clusters_marker_genes <- c("Thrsp", "Cps1", "Cyp2e1", "Cdk1")
ncol = 4
FeaturePlot(Hep_scRNA, features = Endo_scRNA1_clusters_marker_genes, ncol = ncol, order = T, pt.size = 1.2)

Hep_scRNA$anno1 <- Hep_scRNA$seurat_clusters
Hep_scRNA$anno1 <- as.numeric(Hep_scRNA$anno1)
Hep_scRNA$anno1[which(Hep_scRNA$anno1 == "1")] <- "Thrsp+(Regenerating/Normal)"
Hep_scRNA$anno1[which(Hep_scRNA$anno1 == "2")] <- "Cps1+(Stress)"
Hep_scRNA$anno1[which(Hep_scRNA$anno1 == "3")] <- "Cyp2e1+(Pericentral)"
Hep_scRNA$anno1[which(Hep_scRNA$anno1 == "4")] <- "Cdk1+(Proliferating)"
Hep_scRNA$anno1 <- factor(Hep_scRNA$anno1, levels = c("Thrsp+(Regenerating/Normal)", "Cps1+(Stress)", "Cyp2e1+(Pericentral)", "Cdk1+(Proliferating)"))
DimPlot(Hep_scRNA, group.by = "anno1", pt.size = 0.5) 
DimPlot(Hep_scRNA, group.by = "Group", pt.size = 0.5) 


# jjDotPlot
Hep_scRNA_anno1_markers <- c("Thrsp", "Cps1", "Cyp2e1", "Cdk1")
jjDotPlot(object = Hep_scRNA, col.min = 0,base_size = 14, textSize = 14, segWidth = 0.4, bar.width = 3,
               gene = unique(Hep_scRNA_anno1_markers), xtree = F, ytree = F, legend.position = "top",plot.margin = c(0,1,0,0),
               id = 'anno1')

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
  p <- ggplot(data_roe_plot, aes(x = celltype, y = tissue, label = significant)) +
    geom_tile(aes(fill = value), color="gray",size=0.5) +
    coord_equal() + 
    geom_text(color = "black") + 
    scale_fill_viridis_c(option = "C", direction = 1, alpha = 1) +
    theme(axis.title = element_blank(),legend.title = element_blank(), legend.position = "right", axis.text.x = element_text(angle = 35, size = 14, colour = "black", hjust = 1, vjust = 1), axis.text.y = element_text(size = 14, colour = "black"))
  pdf(paste0(plot_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".pdf"), width = 4, height = 4.5)
  print(p)
  dev.off()
}
Roe_Analysis(Hep_scRNA, tissue = "Group", celltype = "anno1", perfix = 'Hep_scRNA_anno1')


##### AddModuleScore: regeneration_receptors 
regeneration_receptors <- c("Egfr","Lrp1","Cd81","Ptprf")
Hep_scRNA <- AddModuleScore(Hep_scRNA, features = list(regeneration_receptors), name = "regeneration_receptors")
library(scales)
library(RColorBrewer)
FeaturePlot(Hep_scRNA, features = "regeneration_receptors1", pt.size = 1.2, order = T) &
  scale_color_gradientn(colours =  rev(brewer.pal(n = 11, name = "RdBu")),
                        guide = "colorbar",
                        limits = c(min(Hep_scRNA$regeneration_receptors1), max(Hep_scRNA$regeneration_receptors1)))



# Hep_scrna_harmony_sce
# END!