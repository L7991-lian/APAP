# 实验室自测的scRNA-seq进行细胞类型间通讯分析
# All为所有细胞类型的scRNA对象
# 我们使用整合多个软件的配受体数据库CellChatDB.mouse_user.rda作为本次分析实验的配受体数据库

### cellchat
load("/data2/lijinlian/APAP_project/cellchat/allcelltype/CellChatDB.mouse_user.rda")
db.new <- CellChatDB.mouse_user
run_cellchat <- function(a1, metadata, 
                         group, 
                         assay="SCT", 
                         # species ="mm", 
                         min.cells = 20,
                         # search = c('Cell-Cell Contact','ECM-Receptor',
                         #            'Secreted Signaling'),
                         color.use = color.use,
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
  cellchat <- projectData(cellchat, PPI.mouse) # raw.use = TRUE 表示使用raw数据，而不使用上一步projectData后的结果
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE, type = "truncatedMean", trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat, thresh = 0.05)
  groupSize <- as.numeric(table(cellchat@idents))
  # par(mfrow = c(1,2), xpd=TRUE)
  pdf(paste(filename_pre,"netVisual_circle_count.pdf", sep = "_"), width = 6, height = 6)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use = color.use)
  dev.off()
  pdf(paste(filename_pre, "netVisual_circle_weight.pdf", sep = "_"), width = 6, height = 6)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = color.use)
  dev.off()
  saveRDS(cellchat, paste(filename_pre, "cellchat_object.Rds", sep = "_"))
  return(cellchat)
}
Celltype_col = c(
  "Hepatocyte" = "#F7931E",
  "Endothelial" = "#F8736A",
  "Macrophage" = "#A9A9A9",
  "T" = "#96CEB4",
  "B" = "#A3A500",
  "Neutrophil" = "#C0C130",
  "NK" = "#9D73C2",
  "Cholangiocyte" = "#B74CAB",
  "HSC" = "#E686C9"
)
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
groups_vector <- names(table(All$Group))
for (i in 1:length(groups_vector)) {
  obj <- subset(All, Group == groups_vector[i])
  filename <- groups_vector[i]
  run_cellchat(obj, 
               metadata = obj@meta.data, 
               group = "celltype", 
               color.use = colors, 
               min.cells = 20, 
               filename_pre = filename)}
# 合并
groups_vector <- c("0h", "6h", "12h", "24h", "48h", "96h")
list_cellchat_object <- list()
for (j in 1:length(groups_vector)) {
  nfile <- paste(groups_vector[j], "cellchat_object.Rds", sep = "_")
  cellchat_obj <- readRDS(nfile)
  cellchat_obj2 <- aggregateNet(cellchat_obj, thresh = 0.05)
  list_cellchat_object[[groups_vector[j]]] <- cellchat_obj2
}
list_cellchat_object <- lapply(list_cellchat_object, function(x){
  cellchat_obj <- netAnalysis_computeCentrality(x, slot.name = "netP")
})
cellchat <- mergeCellChat(list_cellchat_object, add.names = names(list_cellchat_object))


### heatmap: weight to hepatocyte
weight_list <- lapply(c("0h", "6h", "12h", "24h", "48h", "96h"), function(x) {
  cellchat@net[[x]][["weight"]]/max(cellchat@net[[x]][["weight"]])
})
names(weight_list) <- c("0h", "6h", "12h", "24h", "48h", "96h")
w_to_hep_0h <- weight_list[[1]][,1]
w_to_hep_6h <- weight_list[[2]][,1]
w_to_hep_12h <- weight_list[[3]][,1]
w_to_hep_24h <- weight_list[[4]][,1]
w_to_hep_48h <- weight_list[[5]][,1]
w_to_hep_96h <- weight_list[[6]][,1]
w_to_hep_l <- rbind(w_to_hep_0h,w_to_hep_6h,w_to_hep_12h,w_to_hep_24h,w_to_hep_48h,w_to_hep_96h)
rownames(w_to_hep_l) <- c("0h", "6h", "12h", "24h", "48h", "96h")
pheatmap::pheatmap(w_to_hep_l, cluster_rows = F, cluster_cols = F, border_color = "white", fontsize = 14)

### GO enrichment analysis: 6h, endo to hep, signaling
endo_to_hep_6h <- netVisual_bubble(list_cellchat_object[[2]], targets.use = 1, sources.use = 2, remove.isolate = FALSE, return.data = T)
gene.cluster <- data.frame("gene" = c(unique(endo_to_hep_6h$communication$ligand), unique(endo_to_hep_6h$communication$receptor)), "cluster" = c(rep("Ligands", length(unique(endo_to_hep_6h$communication$ligand))), rep("Receptors", length(unique(endo_to_hep_6h$communication$receptor)))))
library(clusterProfiler)
library(enrichplot)
group <- data.frame(gene=gene.cluster$gene,group=gene.cluster$cluster)
Gene_ID <- bitr(group$gene, fromType="SYMBOL",
                toType="ENTREZID",
                OrgDb="org.Mm.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group,
  data=data,
  fun="enrichGO",
  OrgDb="org.Mm.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
data_GO_sim <- clusterProfiler::simplify(data_GO,
                                         cutoff=0.5,
                                         by="p.adjust",
                                         select_fun=min)
data_GO_sim <- setReadable(data_GO_sim, OrgDb = "org.Mm.eg.db")
data_GO_sim_incoming_df_receptors <- data_GO_sim@compareClusterResult
# terms barplot
select_term <- c(
  "GO:0018108", "GO:0050678", "GO:0031589", "GO:0045787", "GO:0019216",
  "GO:0034329", "GO:0031099", "GO:0070508", "GO:0045740", "GO:0045834"
)
select_term <- rev(select_term)
select_term_id <- c()
for (i in select_term) {
  id_t <- grep(paste0("^", i, "$"), data_GO_sim_incoming_df_receptors$ID)
  select_term_id <- c(select_term_id, id_t)
}
select_Go_receptors_terms <- data_GO_sim_incoming_df_receptors[select_term_id, ]
select_Go_receptors_terms$`-log10(p.adjust)` <- -log10(select_Go_receptors_terms$p.adjust)

which(select_Go_receptors_terms$Description == "peptidyl-tyrosine phosphorylation")
select_Go_receptors_terms$Description[10] <- "peptidyl-tyrosine phosphorylation1"
select_Go_receptors_terms <- select_Go_receptors_terms[c(1:4, 11, 5:10),]
select_Go_receptors_terms$Description <- factor(select_Go_receptors_terms$Description, levels =unique(select_Go_receptors_terms$Description))
ggplot(select_Go_receptors_terms,
              aes(x=Description,y=-log10(p.adjust), fill=Cluster)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = c("orange", "purple")) +
  coord_flip() +
  xlab("GO terms") +
  ylab("-log10(p.adjust)") +
  # labs(title = "")+  #设置标题
  theme_bw()+theme(text = element_text(size = 14), axis.text = element_text(size = 16, colour = "black"))


##### 内皮细胞作用于肝细胞的信号中，作用于"Cd81","Egfr","Lrp1","Ptprf"的信号
regeneration_proliferation <- c("Cd81","Egfr","Lrp1","Ptprf")
regeneration_proliferation_df <- findEnrichedSignaling(list_cellchat_object[[2]], features = regeneration_proliferation, pattern ="incoming", thresh = 0.05)
regeneration_proliferation_df_endo_to_hep <- subset(regeneration_proliferation_df, target == "Hepatocyte" & source == "Endothelial")
regeneration_proliferation_df_endo_to_hep_signalling <- regeneration_proliferation_df_endo_to_hep[, c(3:5)]
saveRDS(regeneration_proliferation_df_endo_to_hep_signalling, "regeneration_proliferation_df_endo_to_hep_6h_signalling.Rds")
library(tidyverse)
library(circlize)
circos.par(start.degree = 90)
chordDiagram(
  regeneration_proliferation_df_endo_to_hep_signalling, directional = 1,
  annotationTrack = "grid", link.arr.width = 0.1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  big.gap = 18, diffHeight = -0.01,
  small.gap = 6,
  preAllocateTracks = list(
    track.height = max(strwidth(unlist(chordDiagram_data$source)))
  )
)
circos.track(
  track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$ylim[1], 
      CELL_META$sector.index,  facing = "clockwise", 
      niceFacing = TRUE, adj = c(0, 0.3), cex = 1.15,
    )
  }, bg.border = NA
)
title("Signalling of regeneration receptors from Endothelial to Hepatocyte")
circos.clear()

# END!