### Identify cell type : CACIMAR
# install.packages("devtools")
# devtools::install_github("jiang-junyao/CACIMAR")
# object: seurat object
# Marker: contain two column, Marker and Celltype

library(CACIMAR)
Marker <- read.table('cctmr_landscape.txt',header = T) #343

CACIMAR.annotation <- function(object, name=NULL){
  Marker <- Marker[rownames(Marker)%in%rownames(object),]
  Marker$Marker <- rownames(Marker)
  zfcelltype <- Identify_CellType(object, Marker)
  celltype_cacimar <- as.data.frame(zfcelltype)
  # save table
  celltype_cacimar$method <- "CACIMAR"
  celltype_cacimar$Group <- levels(object@active.ident)
  colnames(celltype_cacimar) <- c("Celltype","Power", "Predicted.Celltype", "method", "Group")
  write.table(celltype_cacimar,file=paste0(name,"_CACIMAR_celltype.txt"), quote=F, sep = "\t")
  # add metadata
  object[["CACIMAR"]] <- "Unknown"
  for (i in 1:length(celltype_cacimar)) {
    object@meta.data[which(object@meta.data$seurat_clusters == celltype_cacimar$Group[i]), "CACIMAR"] <- celltype_cacimar$Predicted.Celltype[i]
  }
  # 
  object[["CACIMAR.power"]] <- "NA"
  for (i in 1:length(celltype_cacimar)) {
    object@meta.data[which(object@meta.data$seurat_clusters == celltype_cacimar$Group[i]), "CACIMAR.power"] <- strsplit(celltype_cacimar$Power[i], ",")[[1]][1] 
  }
  # 
  object[["CACIMAR_celltype"]] <- "Unknown"
  for (i in 1:length(celltype_cacimar)) {
    object@meta.data[which(object@meta.data$seurat_clusters == celltype_cacimar$Group[i]), "CACIMAR_celltype"] <- celltype_cacimar$Celltype[i]
  }
  object[["CACIMAR_celltype_power"]] <- "NA"
  for (i in 1:length(celltype_cacimar)) {
    object@meta.data[which(object@meta.data$seurat_clusters == celltype_cacimar$Group[i]), "CACIMAR_celltype_power"] <- celltype_cacimar$Power[i] 
  }
  return(object)
}
