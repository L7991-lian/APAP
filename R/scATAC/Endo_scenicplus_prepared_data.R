# scenicplus构建内皮细胞亚群调控网络，准备内皮细胞的scRNA-seq对象、scATAC-seq的细胞barcode metadata

# scRNA-seq对象
Endo_scRNA1$annotation6 <- "Lsec_Hbegf"
Endo_scRNA1$annotation6[which(Endo_scRNA1$Clusters != "Apold1+ Lsec(Stress)")] <- "Lsec_other"
library(SeuratDisk)
Idents(Endo_scRNA1) <- Endo_scRNA1$annotation6
Endo_scRNA1@assays$RNA <- Endo_scRNA1@assays$SCT
SaveH5Seurat(Endo_scRNA1, filename = "Endo_scRNA1.h5Seurat")
Convert("Endo_scRNA1.h5Seurat", dest = "h5ad")
# Endo_scRNA1.h5ad

# cell metadata
library(Seurat)
Endo_scATAC1$annotation6 <- "Lsec_Hbegf"
Endo_scATAC1$annotation6[which(Endo_scATAC1$seurat_clusters != "7")] <- "Lsec_other"
Endo_scATAC1_meta <- Endo_scATAC1@meta.data[, c("annotation6", "seurat_clusters", "Group")]
Endo_scATAC1_meta$barcode <- rownames(Endo_scATAC1_meta)
Endo_scATAC1_meta <- Endo_scATAC1_meta[, c("barcode", "annotation6", "seurat_clusters", "Group")]
Endo_scATAC1_meta$barcode <- sub(".*_", "", Endo_scATAC1_meta$barcode)
Endo_scATAC1_meta <- Endo_scATAC1_meta[!duplicated(Endo_scATAC1_meta$barcode),]
write.table(Endo_scATAC1_meta, file = "cell_data.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
