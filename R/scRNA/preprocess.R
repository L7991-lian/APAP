# scRNA-seq preprocess
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
base_path <- "/SoupX/"
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
cluster_doublet_rate_list <- list()
for (i in 1:length(doublet_list)) {
  object <- doublet_list[[i]]
  name <- names(doublet_list[i])
  dt <- table(doublet_list[[i]]$seurat_clusters, doublet_list[[i]]$Doublet)
  cluster <- as.data.frame(table(doublet_list[[i]]$seurat_clusters))
  Doublet_rate <- apply(dt, 2, function(x) x/cluster$Freq)
  cluster_doublet_rate_list[[name]] <- Doublet_rate
  csv.name <- paste0("1.DoubletFinder/",name, "_Cluster.Doublet.rate.csv")
  write.csv(Doublet_rate, csv.name)
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
doublet_list_selfedit <- list()
for (i in 1:length(doublet_list)) {
  object <- doublet_list[[i]]
  name <- names(doublet_list[i])
  nrows <- nrow(cluster_doublet_rate0.5_list[[i]])
  if (nrows == 0) {
    doublet_list_selfedit[[i]] <- doublet_list[[i]]
  }else{
    doublet_list_selfedit[[i]] <- doublet_list[[i]][, -which(doublet_list[[i]]$SCT_snn_res.2 == rownames(cluster_doublet_rate0.5_list[[i]]))]
  }
}
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
scRNA.data$celltype <- scRNA.data@active.ident

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

saveRDS(SeuratObject, "APAP_scRNA_sce.Rds")
