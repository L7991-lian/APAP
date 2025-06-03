# scRNA-seq preprocess, celltype annotation
library(Seurat)
library(Matrix)
library(readr)
library(dplyr)
library(purrr)

# add prefix to cellbarcode for each sample
add_prefix_cellname <- function(base_path, time_points) {
  list_of_data <- list()
  for (time in time_points) {
    data_path <- paste0(base_path, time, "soupX_matrix/")
    data <- Read10X(data_path)
    colnames(data) <- paste0("APAP_", time, "_", colnames(data))
    list_of_data[[paste0("APAP_", time)]] <- data
  }
  return(list_of_data)
}
base_path <- "/data4/jinlianli/APAP/Data.analysis/2022.9.1.SoupX.0.1/SoupX/"
time_points <- c("0h", "6h", "12h", "24h", "48h", "96h")
all_data <- add_prefix_cellname(base_path, time_points)

# Create Seurat Object
APAP_mt0.2_list <- lapply(all_data, function(x) {
  object <- CreateSeuratObject(x, project = names(x), min.cells = 3) 
  # Mitochondrial gene ratio
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")
  meanfeature <- mean(object$nFeature_RNA)
  sdfeature <- sd(object$nFeature_RNA)
  meancount <- mean(object$nCount_RNA)
  sdcount <- sd(object$nCount_RNA)
  thresholdfeature <- meanfeature + 4*sdfeature
  thresholdcount <- meancount + 4*sdcount
  object1 <- subset(object, subset = nFeature_RNA >= 400 & nFeature_RNA < thresholdfeature & percent.mt < 20 & nCount_RNA >= 1000 & nCount_RNA < thresholdcount) #mt<20%
})
names(APAP_mt0.2_list) <- names(all_data)

# doublet 
library(clustree)
PCA_UMAP_list <- list()
set.resolutions <- seq(0.1, 0.8, by = 0.1)
for (l in 1:length(APAP_mt0.2_list)) {
  object <- APAP_mt0.2_list[[l]]
  object_name <- names(APAP_mt0.2_list[l])
  object <- SCTransform(object, vars.to.regress = c("nCount_RNA", "percent.mt"), verbose = FALSE, variable.features.n = 3000)
  object <- RunPCA(object, npcs = 100, verbose = FALSE)
  ElbowPlot <- ElbowPlot(object = object, ndims = 60)
  ggsave(str_c(object_name, "ElbowPlot.pdf"), ElbowPlot, width = 8, height = 8)
  object  <- FindNeighbors(object = object , dims = 1:60, verbose = FALSE)
  object  <- FindClusters(object = object , resolution = set.resolutions, verbose = FALSE) 
  clustree <- clustree(object)
  ggsave(str_c(object_name, "clustree.pdf"), clustree, width = 8, height = 8)
  object  <- RunUMAP(object , dims = 1:60)
  PCA_UMAP_list[[object_name]] <- object
  PCA_UMAP.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = object, reduction = 'umap',label = TRUE, group.by = paste0("SCT_snn_res.", x))
    ggsave(str_c(object_name, "UMAP.SCT_snn_res.", x, ".pdf"), p, width = 7, height = 8)
  })
}
source("./doubletDetect.R")
for (i in seq_along(names(PCA_UMAP_list))) {
  doublet_list[[names(PCA_UMAP_list[i])]] <- doubletDetect(Seurat.object = PCA_UMAP_list[[i]], PCs = 1:60, doublet.rate = 0.061, annotation = "SCT_snn_res.0.7", sct = T)
}
# caculate doublet rate 
# Cluster: resolution=2
for (i in 1:length(doublet_list)) {
  doublet_list[[i]] <- FindClusters(doublet_list[[i]], resolution = c(2))
}
cluster_doublet_rate_list <- list()
for (i in 1:length(doublet_list)) {
  object <- doublet_list[[i]]
  name <- names(doublet_list[i])
  dt <- table(doublet_list[[i]]$seurat_clusters, doublet_list[[i]]$Doublet)
  cluster <- as.data.frame(table(doublet_list[[i]]$seurat_clusters))
  Doublet_rate <- apply(dt, 2, function(x) x/cluster$Freq)
  cluster_doublet_rate_list[[name]] <- Doublet_rate
}
# delete clusters which doublet rate equal and greater than 0.5
cluster_doublet_rate0.5_list <- list()
for (i in 1:length(cluster_doublet_rate_list)) {
  object <- as.data.frame(cluster_doublet_rate_list[[i]][,1])
  colnames(object) <- "Doublet_rate"
  name <- names(cluster_doublet_rate_list[i])
  object <- subset(object, Doublet_rate >= 0.5)
  cluster_doublet_rate0.5_list[[name]] <- object
}

# CACIMAR
source("/data4/jinlianli/APAP/Data.analysis/scRNA.function/CACIMAR.annotation.R")
for (i in 1:length(doublet_list)) {
  object <- doublet_list[[i]]
  name <- names(doublet_list[i])
  doublet_list[[i]] <- CACIMAR.annotation(object, name)
}

# self dedoublet those doublet rate > 0.5 and annotation to two celltypes
# 0h: delect cluster22
doublet_list_self <- list()
doublet_list_self[[names(doublet_list[1])]] <- doublet_list[[1]][, -which(doublet_list[[1]]$SCT_snn_res.2 == "22")]
# 6h
doublet_list_self[[names(doublet_list[2])]] <- doublet_list[[2]]
# 12h: delect cluster23 and cluster25
doublet_list_self[[names(doublet_list[3])]] <- doublet_list[[3]][, -which(doublet_list[[3]]$SCT_snn_res.2 == "23" | doublet_list[[3]]$SCT_snn_res.2 == "25")]
# 24h: delect cluster19
doublet_list_self[[names(doublet_list[4])]] <- doublet_list[[4]][, -which(doublet_list[[4]]$SCT_snn_res.2 == "19")]
# 48h
doublet_list_self[[names(doublet_list[5])]] <- doublet_list[[5]]
# 96h
# delect cluster15
doublet_list_self[[names(doublet_list[6])]] <- doublet_list[[6]][, -which(doublet_list[[6]]$SCT_snn_res.2 == "15")]
doublet_list_selfedit <- doublet_list_self[c("APAP0h", "APAP6h", "APAP12h", "APAP24h", "APAP48h", "APAP96h")]

## left those cells labels as Singlet
Dedoublet_list <- list()
for (i in 1:length(doublet_list_selfedit)) {
  object <- subset(doublet_list_selfedit[[i]], subset = Doublet == "Singlet")
  name <- names(doublet_list_selfedit[i])
  Dedoublet_list[[i]] <- object
}
names(Dedoublet_list) <- names(doublet_list)

# merge and clustering
SeuratObject <- merge(Dedoublet_list[[1]], Dedoublet_list[2:length(Dedoublet_list)])
Group <- substr(SeuratObject$orig.ident, 5, length(SeuratObject$orig.ident))
SeuratObject$Group <- Group
SeuratObject$Group <- factor(SeuratObject$Group, levels = c("0h", "6h", "12h", "24h", "48h", "96h"))

merge.clustering <- function(SeuratObject, scale.var.features.n = 3000, harmony.var = "Group", assay="SCT",  Dims = 20, resolution = 0.1){
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  library(homologene)
  s.genes <- human2mouse(s.genes)
  g2m.genes <- human2mouse(g2m.genes)
  SeuratObject <- CellCycleScoring(SeuratObject, s.features = s.genes$mouseGene, g2m.features = g2m.genes$mouseGene, set.ident = F)
  # cycle Score regress 
  SeuratObject <- SCTransform(SeuratObject, variable.features.n = scale.var.features.n, vars.to.regress = c("S.Score", "G2M.Score"), return.only.var.genes = F)
  SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(SeuratObject), nfeatures.print = 10, npcs=60)
  library(harmony)
  SeuratObject <- SeuratObject %>% RunHarmony(harmony.var, assay.use=assay, verbose = FALSE)
  SeuratObject <- SeuratObject %>%
    RunUMAP(reduction = "harmony", dims = 1:Dims) %>%
    FindNeighbors(reduction = "harmony", dims = 1:Dims) %>%
    FindClusters(reduction.type = "harmony", resolution = resolution) %>% identity()
  return(SeuratObject)
}
SeuratObject <- merge.clustering(SeuratObject, scale.var.features.n = 2000, harmony.var = "Group", assay = "SCT", Dims = 22,  resolution = 0.05)

annotation.scRNA <- function(scRNA.data, new.cluster.ids){
  cell.markers2 <- c("Apoc1", "Apoc3", "Fabp1", "Cyp4a14", "Gstm1", "Alb", # Hepatocyte
                     "Clec4g", "Aqp1", "Bmp2", "Ptprb", "Dpp4", "Egfl7", "F8", "Igfbp7", # Endothelial
                     "Cd79a" ,"Cd79b", "Ebf1", "Blnk", "Bank1", "Fcmr", "Ighd", "Ms4a1", "Pax5", "Igkc", # B
                     "Cd3d", "Trac", "Il7r", "Txk", "Ms4a4b", "Nkg7", "Gzma", "Klrb1c", "Eomes", #T/NK
                     "Lyz2", "Fcer1g", "Clec4f", "Cd5l", "Ctss","Vsig4", # Macrophage
                     "Ifi30", "Cstb", "Cyba", "Cybb", "Aif1", "Cox6a2", # Dendritic
                     "Irf8", "Runx2", "Siglech", "Klk1", "Rnase6", "Ccr9",  # Dendritic
                     "Il1r2", "Il1b", "S100a8", "S100a9", "Lcn2", "Wfdc21",  # Neutrophil
                     "Hbb−bs", "Hba−a1", "Alas2", "Bpgm", # Erythrocyte
                     "Dcn", "Col3a1", "Col1a2", "Col14a1", "Sparc", # HSC
                     "Sox9", "Spp1", "Clu", "Tm4sf4", "Epcam", "Anxa4", "Pkhd1", "Onecut2", "Hnf1b", "Tstd1", "Krt19", "Pdzk1ip1" #Cholangiocyte
  )
  p0 <- DotPlot(scRNA.data, features = unique(cell.markers2)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave("Clusters.markers.Dotplot2.pdf", p0, width = 24, height = 9)
  # annotation
  names(new.cluster.ids) <- levels(scRNA.data$seurat_clusters)
  scRNA.data <- RenameIdents(scRNA.data, new.cluster.ids)
  scRNA.data$CellType <- scRNA.data@active.ident
  # UMAP
  p1 <- DimPlot(scRNA.data, label = T, reduction = "umap")|DimPlot(scRNA.data, label = F, reduction = "umap", group.by = "Group")
  p3 <- DimPlot(scRNA.data, label = T, reduction = "tsne")|DimPlot(scRNA.data, label = F, reduction = "tsne", group.by = "Group")
  ggsave("Celltype.Label.pdf", p1, width = 14, height = 8)
  ggsave("Celltype.tSNE.pdf", p3, width = 14, height = 8)
  return(scRNA.data)
}
## rename idents
new.cluster.ids <- c("Hepatocyte", 
                     "Endothelial",
                     "Macrophage/Dendritic", 
                     "T/NK",
                     "B",
                     "Neutrophil",
                     "Hepatocyte",
                     "Cholangiocyte",
                     "Macrophage/Dendritic",
                     "Erythrocyte",
                     "Hepatocyte", 
                     "HSC"
)
SeuratObject <- annotation.scRNA(SeuratObject, new.cluster.ids)
SeuratObject.delect.Erythrocyte <- subset(SeuratObject, celltype != "Erythrocyte")
SeuratObject.delect.Erythrocyte$celltype <- droplevels(SeuratObject.delect.Erythrocyte$celltype)
SeuratObject.delect.Erythrocyte@active.ident <- droplevels(SeuratObject.delect.Erythrocyte@active.ident)

Dims.Clustering <- function(SeuratObject, Dims = seq(10, 40, 5), resolution.set = 0.1){
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  # clustering
  namefc <- c()
  for (i in Dims) {
    namef <- paste0("Dim", i)
    dir.create(namef)
    namefc <- c(namefc, namef)
  }
  pathwd <- getwd()
  for (m in 1:length(Dims)) {
    setwd(pathwd)
    path <- paste(pathwd, namefc[m], sep = "/")
    setwd(path)
    SeuratObject <- SeuratObject %>%
      RunUMAP(reduction = "harmony", dims = 1:Dims[m]) %>%
      RunTSNE(reduction = "harmony", dims = 1:Dims[m]) %>%
      FindNeighbors(reduction = "harmony", dims = 1:Dims[m]) %>%
      FindClusters(reduction.type = "harmony", resolution = resolution.set) %>% identity()
    titlegg1 <- paste0("Dims", Dims[m], ".Resolution", resolution.set)
    titlegg2 <- paste0("Dims", Dims[m], ".Resolution", resolution.set)
    p11 <- DimPlot(SeuratObject, reduction = "umap", label = T) + ggtitle(titlegg1)|DimPlot(SeuratObject, reduction = "umap", group.by = "Group", label = F)
    p12 <- DimPlot(SeuratObject, reduction = "tsne", label = T) + ggtitle(titlegg2)|DimPlot(SeuratObject, reduction = "tsne", group.by = "Group", label = F)
    ggsavename1 <- paste0("Try.UMAP.", "Dims", Dims[m], "Resolution", resolution.set, ".pdf")
    ggsavename2 <- paste0("Try.tSNE.", "Dims", Dims[m], "Resolution", resolution.set, ".pdf")
    ggsave(ggsavename1, plot = p11, width = 14, height = 8)
    ggsave(ggsavename2, plot = p12, width = 14, height = 8)
  }
  return(SeuratObject)
}
SeuratObject <- Dims.Clustering(SeuratObject.delect.Erythrocyte, Dims = seq(28), resolution.set = 0.05)

# check celltype markers expression
cell.markers <- c("Apoc1", "Apoc3", "Fabp1", "Cyp4a14", "Gstm1", "Alb", # Hepatocyte
                  "Clec4g", "Aqp1", "Bmp2", "Ptprb", "Dpp4", "Egfl7", "F8", # Endothelial
                  "Cd79a" ,"Cd79b", "Ebf1", "Blnk", "Bank1", "Fcmr", "Ighd", "Ms4a1", "Pax5", "Igkc", # B
                  "Cd3d", "Cd3g", #T
                  "Klre1", "Klrb1c", "Eomes","Gzma",  #NK
                  "Lyz2", "Fcer1g", "Clec4f", "Cd5l", "Ctss","Vsig4", # Macrophage
                  "Ifi30", "Cstb", "Cyba", "Cybb", "Aif1", "Cox6a2", # Dendritic
                  "Irf8", "Runx2", "Siglech", "Klk1", "Rnase6", "Ccr9",  # Dendritic
                  "Il1r2", "Il1b", "S100a8", "S100a9", "Lcn2", "Wfdc21",  # Neutrophil
                  "Hbb−bs", "Hba−a1", "Alas2", "Bpgm", # Erythrocyte
                  "Dcn", "Col3a1", "Col1a2", "Col14a1", "Sparc", # HSC
                  "Sox9", "Spp1", "Clu", "Hnf1b", "Tstd1", "Krt19" #Cholangiocyte
)
p0 <- DotPlot(SeuratObject, features = unique(cell.markers)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p0
# annotation
new.cluster.ids <- c("Hepatocyte", #0
                     "Endothelial", #1
                     "Macrophage/Dendritic", #2
                     "T/NK", #3
                     "B", #4
                     "Neutrophil", #5
                     "NK", #6
                     "Cholangiocyte", #7
                     "Macrophage/Dendritic", #8
                     "Endothelial", #9
                     "Hepatocyte", #10
                     "HSC" #11
)
names(SeuratObject) <- levels(SeuratObject$seurat_clusters)
SeuratObject <- RenameIdents(SeuratObject, new.cluster.ids)
SeuratObject$celltype <- SeuratObject@active.ident
celltype.markers <- c("Alb", # Hepatocyte
                      "Ptprb", # Endothelial
                      "Cybb", # Macrophage
                      "Cd3d", #T
                      "Ebf1", # B
                      "S100a9", # Neutrophil
                      "Klre1", #NK
                      "Sox9", #Cholangiocyte
                      "Dcn" # HSC
)
DotPlot(SeuratObject, group.by = "celltype", features = celltype.markers) 
celltype.markers <- c("Alb", # Hepatocyte
                      "Ptprb", # Endothelial
                      "Cybb", # Macrophage
                      "Cd3d", #T
                      "Ebf1", # B
                      "S100a9", # Neutrophil
                      "Klre1", #NK
                      "Sox9", #Cholangiocyte
                      "Dcn" # HSC
)
DotPlot(SeuratObject, group.by = "celltype", features = celltype.markers) 
pdf("scRNA_celltype_markers_dotplot.pdf")
jjDotPlot(object = SeuratObject, 
          col.min = 0,
          anno = T,
          base_size = 14,
          textSize = 14, 
          id = 'celltype',
          gene = celltype.markers,
          xtree = F, ytree = F, 
          legend.position = "right")
dev.off()
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
pdf("scrna_UMAP.pdf")
DimPlot(SeuratObject, group.by = "celltype", cols = colors)+ theme(legend.title = element_blank(), legend.text = element_text(size = 14)) + ggtitle("")| DimPlot(All_obj, group.by = "Group", cols = colors)+ theme(legend.title = element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()
saveRDS(SeuratObject, "/data4/jinlianli/APAP/Data.analysis/Object/2.Merge/All.Rds")
