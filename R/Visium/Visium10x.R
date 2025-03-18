# APAP, 10x Visium 空间转录组数据分析
# 10x Visium中每个spot长宽约50um，因此每个spot测序数据包含多个细胞。

# linux 实现
# 将 GSM8404656_APAP24h_barcodes.tsv.gz 重命名为 barcodes.tsv.gz。
# for file in GSM8404656_APAP24h_*.gz; do
# mv "$file" "${file#GSM8404656_APAP24h_}"
# done

# raw_feature_bc_matrix:
# 28K 2 20 2020 barcodes.tsv.gz
# 298K 2 20 2020 features.tsv.gz
# 18M 2 20 2020 matrix.mtx.gz
# spatial:
# 1.6M 2 20 2020 aligned_fiducials.jpg
# 1.8M 2 20 2020 detected_tissue_image.jpg
# 163B 2 20 2020 scalefactors_json.json
# 4.6M 2 20 2020 tissue_hires_image.png
# 496K 2 20 2020 tissue_lowres_image.png
# 180K 2 20 2020 tissue_positions_list.csv

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)
library(Matrix)
# library(qs)
library(BiocParallel)
# library(RCTD)
library(spacexr)

#### 0h
ct = Read10X("./raw_feature_bc_matrix/0h/") 
## 记得gunzip
img <- Read10X_Image(image.dir = file.path("./spatial/", "0h"), filter.matrix = TRUE) 
sceP <- CreateSeuratObject(counts = ct, assay = "Spatial", project = "0h") 
img <- img[Cells(x = sceP)] 
DefaultAssay(sceP = img) <- "Spatial" 
sceP[['slice_0h']] <- img
sceP <- subset(sceP, nCount_Spatial > 5000)
sceP <- SCTransform(sceP, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

#### 3h
ct_3h = Read10X("./raw_feature_bc_matrix/3h/") 
img_3h <- Read10X_Image(image.dir = file.path("./spatial/", "3h"), filter.matrix = TRUE) 
sceP_3h <- CreateSeuratObject(counts = ct_3h, assay = "Spatial", project = "3h") 
img_3h <- img_3h[Cells(x = sceP_3h)] 
DefaultAssay(sceP_3h = img_3h) <- "Spatial" 
sceP_3h[['slice_3h']] <- img_3h
sceP_3h<-subset(sceP_3h, nCount_Spatial > 5000)
sceP_3h <- SCTransform(sceP_3h, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

#### 6h
ct_6h = Read10X("./raw_feature_bc_matrix/6h/") 
img_6h <- Read10X_Image(image.dir = file.path("./spatial/", "6h"), filter.matrix = TRUE) 
sceP_6h <- CreateSeuratObject(counts = ct_6h, assay = "Spatial", project = "6h") 
img_6h <- img_6h[Cells(x = sceP_6h)] 
DefaultAssay(sceP_6h = img_6h) <- "Spatial" 
sceP_6h[['slice_6h']] <- img_6h
sceP_6h<-subset(sceP_6h, nCount_Spatial > 5000)
sceP_6h <- SCTransform(sceP_6h, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

#### 24h
ct_24h = Read10X("./raw_feature_bc_matrix/24h/") 
img_24h <- Read10X_Image(image.dir = file.path("./spatial/", "24h"), filter.matrix = TRUE) 
sceP_24h <- CreateSeuratObject(counts = ct_24h, assay = "Spatial", project = "24h") 
img_24h <- img_24h[Cells(x = sceP_24h)] 
DefaultAssay(sceP_24h = img_24h) <- "Spatial" 
sceP_24h[['slice_24h']] <- img_24h
sceP_24h <- subset(sceP_24h, nCount_Spatial > 5000)
sceP_24h <- SCTransform(sceP_24h, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

##### merge
APAP_merge <- merge(
  x = sceP,
  y = list(sceP_3h, sceP_6h, sceP_24h),
  add.cell.ids = c("0h", "3h", "6h", "24h")
)
DefaultAssay(APAP_merge) <- "SCT"
APAP_merge$orig.ident <- factor(APAP_merge$orig.ident, levels = c("0h", "3h", "6h", "24h"))

#### Hbegf_Egfr_score
genes=c("Hbegf", "Egfr")
Hbegf_Egfr_interaction_all_cell <- as.data.frame(t(as.matrix(APAP_merge@assays$SCT@counts[genes,])))
Hbegf_Egfr_interaction_all_cell$Hbegf_Egfr_score <- Hbegf_Egfr_interaction_all_cell$Hbegf*Hbegf_Egfr_interaction_all_cell$Egfr
APAP_merge$Hbegf_Egfr <- Hbegf_Egfr_interaction_all_cell$Hbegf_Egfr_score
SpatialFeaturePlot(APAP_merge, features = c("Hbegf_Egfr"), pt.size.factor = 2.5)

# DONE!
