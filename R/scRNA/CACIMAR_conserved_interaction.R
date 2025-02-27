# CACIMAR计算细胞间相互作用的保守性得分
# list_cellchat_object_human2[[2]]是公共人的数据的损伤后的scRNA-seq样本，list_cellchat_object[[2]]是实验室自测的6小时的小鼠scRNA-seq样本

### 提取出两个物种的cellchat的细胞间配受体互作结果，用于下面CACIMAR计算保守性的输入
# human
APAP_human_all <- netVisual_bubble(list_cellchat_object_human2[[2]], 
                                   sources.use = c(1:5, 7,8), targets.use = c(1:5, 7,8),  
                                   angle.x = 45, return.data = T, 
                                   remove.isolate = FALSE)
APAP_human_all_signalling <- APAP_human_all[["communication"]]
APAP_human_all_signalling_u <- APAP_human_all_signalling[, 1:5]
APAP_human_all_signalling_u <- APAP_human_all_signalling_u[, c(3,4,1,2,5)]
write.table(APAP_human_all_signalling_u, "APAP_human_signalling.txt")

# mouse wj
APAP_mousewj_all <- netVisual_bubble(list_cellchat_object[[2]], 
                                     sources.use = c(1:5, 8,9), targets.use = c(1:5, 8,9),  
                                     angle.x = 45, return.data = T, 
                                     remove.isolate = FALSE)
APAP_mousewj_all_signalling <- APAP_mousewj_all[["communication"]]
APAP_mousewj_all_signalling_u <- APAP_mousewj_all_signalling[, 1:5]
APAP_mousewj_all_signalling_u <- APAP_mousewj_all_signalling_u[, c(3,4,1,2,5)]
write.table(APAP_mousewj_all_signalling_u, "APAP_mousewj_signalling.txt")

### 计算保守性
# 输入保守的细胞类型（因为肝脏主要细胞类型在人和小鼠中，比较明确，这里直接指定保守细胞类型。如果需要，有可以使用CACIMAR包的函数计算保守的细胞类型）
library(CACIMAR)
conserved_cell_types_mm_zf <- data.frame("mm" = c(names(table(list_cellchat_object[[2]]@idents))[-c(6,7)]), "hs" = c(names(table(list_cellchat_object_human2[[2]]@idents))[c(1:5, 7,8)]))
conserved_cell_types_mm_zf$mm <- paste0("mm", conserved_cell_types_mm_zf$mm)
conserved_cell_types_mm_zf$hs <- paste0("hs", conserved_cell_types_mm_zf$hs)
conserved_cell_types_mm_zf
#               mm               hs
# 1    mmHepatocyte    hsHepatocytes
# 2   mmEndothelial     hsEndothelia
# 3    mmMacrophage            hsMPs
# 4             mmT        hsT cells
# 5             mmB        hsB cells
# 6 mmCholangiocyte hsCholangiocytes
# 7           mmHSC     hsMesenchyme

load(system.file("extdata", "HOM_matrix.rda", package = "CACIMAR"))
load(system.file("extdata", "species_names_ref.rda",package = "CACIMAR"))
species_names_ref$sy_name[which(species_names_ref$sy_name == "hu")] <- "hs"

# colnames 可以不同，但是colname的顺序应该与参考数据一致
colnames(APAP_mousewj_all_signalling_u) <- colnames(sp1_ccc)
colnames(APAP_human_all_signalling_u) <- colnames(sp1_ccc)
# 将每个物种的细胞间相互作用除以对应物种总的相互作用，将相互作用得分作一个标准化
APAP_mousewj_all_signalling_u$LRscore_scale <- APAP_mousewj_all_signalling_u$prob/sum(na.omit(APAP_mousewj_all_signalling_u$prob))  ####na.omit 很重要
APAP_human_all_signalling_u$LRscore_scale <- APAP_human_all_signalling_u$prob/sum(na.omit(APAP_human_all_signalling_u$prob))

# 鉴定每对相互作用的细胞的保守的细胞间相互作用，并计算保守性得分
conserved_result_mm_zf <- CACIMAR::Identify_Conserved_CCI1(species1_cci=APAP_mousewj_all_signalling_u,
                                                  species2_cci=APAP_human_all_signalling_u,
                                                  conserved_cell_types_df=conserved_cell_types_mm_zf,
                                                  species_name1 = "mm",
                                                  species_name2 = "hs", 
                                                  HOM_matrix = HOM_matrix,
                                                  species_names_ref = species_names_ref)
cci_conserved_Weights_table_mm_zf <- CACIMAR::Caculate_cell_pair_cci_score(conserved_result_df=conserved_result_mm_zf,
                                                                  species1_cci = APAP_mousewj_all_signalling_u,
                                                                  species2_cci = APAP_human_all_signalling_u,
                                                                  conserved_cell_types_df = conserved_cell_types_mm_zf,
                                                                  species_name1 = "mm",
                                                                  species_name2 = "hs")
# 将上面计算的保守性得分，转换成符合ChordDiagram绘图输入数据需要的格式，然后绘制保守性得分ChordDiagram
cci_data_mm_zf <- Make_ChordDiagram_data(cci_conserved_Weights_table = cci_conserved_Weights_table_mm_zf, species_name1 = "mm", species_name2 = "hs")
cci_data_mm_zf_to_hep <- subset(cci_data_mm_zf, rep_target == "Hepatocyte/Hepatocytes")
cci_data_mm_zf_to_hep <- subset(cci_data_mm_zf_to_hep, rep_source != "Hepatocyte/Hepatocytes")
ChordDiagram(net = cci_data_mm_zf_to_hep, 
             filename = paste("mm", "hs", "_Conserved_interaction_score_to_Hepatocyte_between_mouse_and_human_APAP_chordDiagram.pdf.pdf", sep = "_"), 
             brewer_pal_used = "Dark2", link_colors_threshold = 0.9, Score_factor_size = 150, picture_width = 6, picture_height = 6)

######### 提取出保守相互作用信号LR，然后将作用于肝细胞的再生受体的信号绘制在ChordDiagram展示
cci_conserved_result_in_hs <- conserved_result_mm_zf[["sp2_cci_conserved_df"]]
cci_conserved_result_in_hs_endo_to_hep <- subset(cci_conserved_result_in_hs, c_pair == "hsEndothelia_hsHepatocytes")

cci_conserved_result_in_hs_endo_to_hep_chordD <- cci_conserved_result_in_hs_endo_to_hep[, c("ligand","receptor", "prob")]
cci_conserved_regeneration_receptors_result_in_hs_endo_to_hep_chordD <- subset(cci_conserved_result_in_hs_endo_to_hep_chordD, receptor %in% c("CD81","EGFR","LRP1","PTPRF"))
chordDiagram(
  cci_conserved_regeneration_receptors_result_in_hs_endo_to_hep_chordD, directional = 1,
  # grid.col = grid.col, 
  annotationTrack = "grid", link.arr.width = 0.1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  # self.link = "1",
  # link.auto = F,
  # group = group, 
  big.gap = 15,
  diffHeight = -0.01,
  # small.gap = 4,
  preAllocateTracks = list(
    track.height = max(strwidth(unlist(chordDiagram_data$source)))
  )
)
# 在空轨迹中放置文本标签
circos.track(
  track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$ylim[1], 
      CELL_META$sector.index,  facing = "clockwise", 
      niceFacing = TRUE, adj = c(0, 0.3), cex = 1.15,
    )
  }, bg.border = NA
)
title("Conserved regeneration receptors signallings from endothelial to hepatocyte in human data")
circos.clear()


