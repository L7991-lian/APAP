# SoupX: ambient RNA
library(SoupX)
#Load data and estimate soup profile
filter_0h <- "/0h/outs/filtered_feature_bc_matrix/"
raw_0h <- "/0h/outs/raw_feature_bc_matrix/"
filter_6h <- "/6h/outs/filtered_feature_bc_matrix/"
raw_6h <- "/6h/outs/raw_feature_bc_matrix/"
filter_12h <- "/12h/outs/filtered_feature_bc_matrix/"
raw_12h <- "/12h/outs/raw_feature_bc_matrix/"
filter_24h <- "/24h/outs/filtered_feature_bc_matrix/"
raw_24h <- "/24h/outs/raw_feature_bc_matrix/"
filter_48h <- "/48h/outs/filtered_feature_bc_matrix/"
raw_48h <- "/48h/outs/raw_feature_bc_matrix/"
filter_96h <- "/96h_2/outs/filtered_feature_bc_matrix/"
raw_96h <- "/96h_2/outs/raw_feature_bc_matrix/"

toc_filter <- c(filter_0h, filter_6h, filter_12h, filter_24h, filter_48h, filter_96h)
tod_filter <- c(raw_0h, raw_6h, raw_12h, raw_24h, raw_48h, raw_96h)
toc_tod_data.frame <- data.frame(toc = toc_filter, tod = tod_filter)

out_list <- list()
sc_list <- list()
run_soupx <- function(toc,tod,rho=NULL) {
  library(Seurat)
  toc <- Read10X(toc,gene.column=1)
  tod <- Read10X(tod,gene.column=1)
  tod <- tod[rownames(toc),]
  ##聚类
  all <- toc
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 100, verbose = F) 
  all <- FindNeighbors(all, dims = 1:50)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:50)
  matx <- all@meta.data
  #####
  sc = SoupChannel(tod, toc)
  sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
  if (is.null(rho)) {
    tryCatch(
      {sc = autoEstCont(sc)}, 
      error=function(e) {
        print("autoEstCont Error !")
        sc = setContaminationFraction(sc, 0.1)}
    )
  }else{
    sc = setContaminationFraction(sc, rho)
  }
  #校正矩阵
  out = adjustCounts(sc)
  #保存两个矩阵文件
  Rds_name <- paste0(names[i], ".Rds")
  out_list[names[i]] <- out
  sc_list[names[i]] <- sc
  saveRDS(sc, Rds_name)
  DropletUtils:::write10xCounts(paste0("./", names[i], "soupX_matrix"), out, version="3")
}
# RUN
for (i in 1:nrow(toc_tod_data.frame)) {
  names <- c("0h", "6h", "12h", "24h", "48h", "96h")
  run_soupx(toc_tod_data.frame[i, "toc"], toc_tod_data.frame[i, "tod"], rho = 0.1)
}
