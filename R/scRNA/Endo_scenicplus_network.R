# Endothelial network with scenicplus

################################ prepared data
# RNA
table(Endo_scRNA1$annotation6)
# Apold1+ Lsec     EC_other 
# 1939         5207 
library(SeuratDisk)
Idents(Endo_scRNA1) <- Endo_scRNA1$annotation6
# we used the SCT data
Endo_scRNA1@assays$RNA <- Endo_scRNA1@assays$SCT
SaveH5Seurat(Endo_scRNA1, filename = "Endo_scRNA1.h5Seurat")
Convert("Endo_scRNA1.h5Seurat", dest = "h5ad")
# Endo_scRNA1.h5ad

# ATAC
table(Endo_scATAC1$annotation6)
# Apold1+ Lsec     EC_other 
#         193         4872 
Endo_scATAC1_meta <- Endo_scATAC1@meta.data[, c("annotation6", "seurat_clusters", "Group")]
Endo_scATAC1_meta$barcode <- rownames(Endo_scATAC1_meta)
Endo_scATAC1_meta <- Endo_scATAC1_meta[, c("barcode", "annotation6", "seurat_clusters", "Group")]
Endo_scATAC1_meta$barcode <- sub(".*_", "", Endo_scATAC1_meta$barcode)
Endo_scATAC1_meta <- Endo_scATAC1_meta[!duplicated(Endo_scATAC1_meta$barcode),]
write.table(Endo_scATAC1_meta, file = "cell_data.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


################################ downstream analysis and plot






