# 公共人的scRNA-seq进行细胞类型间通讯分析
# Human_scRNA为APAP诱导肝损伤的公共的人的scRNA-seq数据的所有细胞类型的seurat对象
# 我们使用整合多个软件的配受体数据库cellchat_Human_db.new.Rds作为本次分析实验的配受体数据库
# cellchat_Human_db.new.Rds是小鼠那个通过mouse2human函数转换的

db.new <- readRDS("/data2/lijinlian/APAP_project/cellchat/human2024_2/cellchat_Human_db.new.Rds")
Idents(Human_scRNA) <- Human_scRNA$anno
Human_scRNA$anno <- factor(Human_scRNA$anno, levels = c("Hepatocytes", "Endothelia", "MPs", "T cells", "B cells", "Plasma", "Cholangiocytes", "Mesenchyme", "ILCs", "Unidentified"))

library(CellChat)
library(Seurat)
library(ggplot2)
library(liana)
library(tidyverse)
run_cellchat <- function(a1, metadata, 
                         group, 
                         assay="RNA", 
                         # species ="hs",
                         min.cells = 20,
                         # search = c('Cell-Cell Contact',
                         #            # 'ECM-Receptor',
                         #            'Secreted Signaling'),
                         filename_pre){
  library(CellChat)
  cellchat <- createCellChat(object = a1, meta = metadata, group.by = group, assay = assay)
  # if (species =='mm') {
  #   db = CellChatDB.mouse
  # }
  # if (species =='hs') {
  #   db = CellChatDB.human
  # }
  # if (species =='zf') {
  #   db = CellChatDB.zebrafish
  # }
  # CellChatDB <- db
  # CellChatDB.use <- subsetDB(CellChatDB, search = search)
  # cellchat@DB <- CellChatDB.use
  cellchat@DB <- db.new
  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 20)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE, type = "truncatedMean", trim = 0)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat, thresh = 0.05)
  groupSize <- as.numeric(table(cellchat@idents))
  # par(mfrow = c(1,2), xpd=TRUE)
  pdf(paste(filename_pre,"netVisual_circle_count.pdf", sep = "_"), width = 6, height = 6)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", alpha.edge = 1)
  dev.off()
  pdf(paste(filename_pre, "netVisual_circle_weight.pdf", sep = "_"), width = 6, height = 6)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",alpha.edge = 1)
  dev.off()
  saveRDS(cellchat, paste(filename_pre, "cellchat_object.Rds", sep = "_"))
  return(cellchat)
}
groups_vector <- names(table(Human_scRNA$condition))
for (i in 1:length(groups_vector)) {
  obj <- subset(Human_scRNA, condition == groups_vector[i])
  print(table(obj$condition))
  filename <- groups_vector[i]
  run_cellchat(obj, 
               metadata = obj@meta.data, 
               group = "anno", 
               min.cells = 5, 
               filename_pre = filename)}


### 合并
list_cellchat_object_human2 <- list()
for (j in 1:length(groups_vector)) {
  nfile <- paste(groups_vector[j], "cellchat_object.Rds", sep = "_")
  cellchat_obj <- readRDS(nfile)
  cellchat_obj2 <- aggregateNet(cellchat_obj, thresh = 0.05)
  list_cellchat_object_human2[[groups_vector[j]]] <- cellchat_obj2
}
list_cellchat_object_human2 <- lapply(list_cellchat_object_human2, function(x){
  cellchat_obj <- netAnalysis_computeCentrality(x, slot.name = "netP")
})
cellchat_human2 <- mergeCellChat(list_cellchat_object_human2, add.names = names(list_cellchat_object_human2), cell.prefix = T)


### 肝脏非实质细胞作用于肝细胞的weight的heatmap
H_W1 <- Healthy_APAP_cellchat@net[["weight"]]/max(Healthy_APAP_cellchat@net[["weight"]])
A_W1 <- human_APAP_cellchat@net[["weight"]]/max(human_APAP_cellchat@net[["weight"]])
W_df <- rbind(H_W1[,1], A_W1[,1])
rownames(W_df) <- c("Healthy", "APAP")
pheatmap::pheatmap(W_df[,2:10], 
                   cluster_rows = F, cluster_cols = F, 
                   border_color = "white", main = "Interaction score to hepatocytes" ,
                   fontsize = 14, cellheight = 25)

# END!