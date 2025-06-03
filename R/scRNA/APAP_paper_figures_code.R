# 此代码是，APAP与肝再生研究项目，截至2025年6月，绘制的所有分析结果的图片代码。现整理于本人GitHub仓库中，公开给大家参考。
# 统一说明：data2下的文件夹在新服务器，data4下的文件夹在旧服务器

# Figure1
# B
All_scrna <- readRDS("/data2/lijinlian/APAP/All.Rds")
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
pdf("scrna_UMAP.pdf")
DimPlot(All_scrna, group.by = "celltype", cols = colors)+ theme(legend.title = element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()
library(scRNAtoolVis)
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
pdf("scRNA_celltype_markers_dotplot.pdf")
jjDotPlot(object = All_scrna, 
          col.min = 0,
          anno = T,
          base_size = 14,
          textSize = 14, 
          id = 'celltype',
          gene = celltype.markers,
          xtree = F, ytree = F, 
          legend.position = "right")
dev.off()

# C
All_scatac <- readRDS("/data2/lijinlian/APAP/scATAC.data_scale.Rds")
All_scatac$predicted.id <- as.character(All_scatac$predicted.id)
All_scatac$predicted.id <- factor(All_scatac$predicted.id, c("Hepatocyte", "Endothelial", "Macrophage", "T", "B", "Neutrophil", "NK", "Cholangiocyte", "HSC"))         
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
pdf("scatac_UMAP.pdf")
DimPlot(All_scatac, group.by = "predicted.id", cols = colors)+ theme(legend.title = element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
dev.off()
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
pdf("scATAC_celltype_markers_dotplot.pdf")
jjDotPlot(object = All_scatac, 
          col.min = 0,
          anno = T,
          base_size = 14,
          textSize = 14, 
          id = 'predicted.id',
          gene = celltype.markers,
          xtree = F, ytree = F, 
          legend.position = "right")
dev.off()

# D
a <- readRDS("/home/lijinlian/scRNA_celltype_Group_DEG/scRNA_celltype_Correlation_pheatmap_annobar.Rds")
p <- pheatmap::pheatmap(a$cor,
                        cluster_rows = F, 
                        cluster_cols = F, 
                        border_color = "white",
                        annotation_col = a$annotation_col, 
                        annotation_row = a$annotation_col,
                        annotation_colors = a$ann_colors,
                        show_rownames = F, show_colnames = F, 
                        cellwidth = 15, cellheight = 15,
                        treeheight_row = 0, 
                        main = "Correlation between samples of celltypes",
                        color = colorRampPalette(colors = c("blue","white","red3"))(100))
ggsave("scRNA_celltype_Correlation_pheatmap_annobar.pdf", p, width = 15, height = 15)

result_df_sig <- read.table("hsc_hep_endo_maf_chol_6h_sigp.adjust0.05_terms_uncluster.txt")
# terms pheatmap
set_terms1 <- c(
  #endo
  # "GO:0006282", #regulation of DNA repair
  "GO:0006914", #autophagy
  "GO:0001667", # ameboidal-type cell migration
  "GO:0010631", #epithelial cell migration
  "GO:0033002", #muscle cell proliferation
  "GO:0035966", #response to topologically incorrect protein
  "GO:0071559", #response to transforming growth factor beta
  "GO:0198738", #cell-cell signaling by wnt
  "GO:0042058", #regulation of epidermal growth factor receptor signaling pathwa
  "GO:0031099", #regeneration
  "GO:0042246", #tissue regeneration
  # "GO:0050678", #regulation of epithelial cell proliferation
  "GO:0072575", #epithelial cell proliferation involved in liver morphogenesis
  "GO:0072574", #hepatocyte proliferation
  # "GO:0032675", #regulation of interleukin-6 production
  #HSC
  # "GO:0001666", #response to hypoxia
  # "GO:0007263", # nitric oxide mediated signal transduction
  # "GO:0043114", #regulation of vascular permeability
  # "GO:0010632", #regulation of epithelial cell migration
  #Hep
  "GO:0019216", #regulation of lipid metabolic process
  "GO:0009062", #fatty acid catabolic process
  "GO:0006635", #fatty acid beta-oxidation
  "GO:0006006", #glucose metabolic process
  "GO:0030301", #cholesterol transport
  "GO:0050817", # coagulation
  "GO:0006956", #complement activation
  #Mac
  "GO:0140236", #translation at presynapse
  "GO:0046651", #lymphocyte proliferation
  "GO:0032943", #mononuclear cell proliferation
  "GO:0050863", #regulation of T cell activation
  "GO:0042113", #B cell activation
  "GO:0007159", #leukocyte cell-cell adhesion
  "GO:0035924", #cellular response to vascular endothelial growth factor stimulus
  "GO:0032623", #interleukin-2 production
  #Chol
  "GO:0042273", #ribosomal large subunit biogenesis
  "GO:0008406", #gonad development
  "GO:0045682", #regulation of epidermis development
  # "GO:0045444", #fat cell differentiation
  "GO:0048011" #neurotrophin TRK receptor signaling pathway
)
term_ids <- c()
for (i in seq_along(set_terms1)) {
  id <- grep(set_terms1[i], result_df_sig$ID)
  term_ids <- c(term_ids, id)
}
hsc_hep_endo_maf_6h_set_terms_df <- result_df_sig[term_ids,]
hsc_hep_endo_maf_6h_set_terms_df <- hsc_hep_endo_maf_6h_set_terms_df[,c("Description","p.adjust","cluster")]
hsc_hep_endo_maf_6h_set_terms_df$log10_p.adjust <- -log10(hsc_hep_endo_maf_6h_set_terms_df$p.adjust)
hsc_hep_endo_maf_6h_set_terms_df <- hsc_hep_endo_maf_6h_set_terms_df[, c("cluster", "log10_p.adjust", "Description")]
df_wide <- hsc_hep_endo_maf_6h_set_terms_df %>%
  pivot_wider(names_from = cluster, values_from = c(log10_p.adjust))
rownames(df_wide) <- df_wide$Description
rn <- df_wide$Description
df_wide <- df_wide[, -1]
rownames(df_wide) <- rn
ph3 <- pheatmap::pheatmap(df_wide, border_color = "white", fontsize = 12, cluster_rows = F, cluster_cols = F, cellwidth = 30)
ggsave("hsc_hep_endo_maf_Chol_notcluster_6h_set_terms_pheatmap.pdf", ph3, width = 8, height = 6.5)

# E
library(ggplot2)
cellcycle_score_df_normalized <- read.table("/data2/lijinlian/APAP_project/coembed_scrna_scatac/S.cellcycle_score_df_normalized_width_table.txt")
cellcycle_score_df_long <- cellcycle_score_df_normalized %>%
  pivot_longer(cols = colnames(cellcycle_score_df_normalized)[1:9], names_to = "Celltype", values_to = "S.Score")
cellcycle_score_df_long$Celltype <- factor(cellcycle_score_df_long$Celltype, levels = levels(All$celltype))
cellcycle_score_df_long$Group <- factor(cellcycle_score_df_long$Group, levels = levels(All$Group))
cellcycle_score_picture <- ggplot(cellcycle_score_df_long, aes(x = Group, y = S.Score, group = Celltype, color = Celltype)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Group", y = "S.Score") +
  theme_minimal() +
  scale_color_manual(values =  as.vector(Celltype_col))+
  theme(text = element_text(size = 14), axis.title.x = element_blank())+
  facet_wrap(~Celltype, scales = "free_y")
ggsave("All_celltype_S.Score_line_plot.pdf", cellcycle_score_picture, width = 9, height = 7)

cellcycle_score_df_normalized <- read.table("/data2/lijinlian/APAP_project/coembed_scrna_scatac/G2M.cellcycle_score_df_normalized_width_table.txt")
cellcycle_score_df_long <- cellcycle_score_df_normalized %>%
  pivot_longer(cols = colnames(cellcycle_score_df_normalized)[1:9], names_to = "Celltype", values_to = "G2M.Score")
cellcycle_score_df_long$Celltype <- factor(cellcycle_score_df_long$Celltype, levels = levels(All$celltype))
cellcycle_score_df_long$Group <- factor(cellcycle_score_df_long$Group, levels = levels(All$Group))
cellcycle_score_picture <- ggplot(cellcycle_score_df_long, aes(x = Group, y = G2M.Score, group = Celltype, color = Celltype)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Group", y = "G2M.Score") +
  theme_minimal() +
  scale_color_manual(values =  as.vector(Celltype_col))+
  theme(text = element_text(size = 14), axis.title.x = element_blank())+
  facet_wrap(~Celltype, scales = "free_y")
ggsave("All_celltype_G2M.Score_line_plot.pdf", cellcycle_score_picture, width = 9, height = 7)
# phx
S.cellcycle_score_df_raw_width_table <- read.csv("/data4/jinlianli/public_liver_hep/S.cellcycle_score_df_raw_width_table.txt", sep="")
S.cellcycle_score_df_raw_width_table_scale <- S.cellcycle_score_df_raw_width_table[,2:7]
S.cellcycle_score_df_raw_width_table_scale_normalized <- apply(S.cellcycle_score_df_raw_width_table_scale, 1, function(x) x - x[1])
library(ggplot2)
S.cellcycle_score_df_raw_width_table_scale_normalized <- as.data.frame(S.cellcycle_score_df_raw_width_table_scale_normalized)
S.cellcycle_score_df_raw_width_table_scale_normalized$group <- rownames(S.cellcycle_score_df_raw_width_table_scale_normalized)
S.cellcycle_score_df_long <- S.cellcycle_score_df_raw_width_table_scale_normalized %>%
  pivot_longer(cols = colnames(S.cellcycle_score_df_raw_width_table_scale_normalized)[1:13], names_to = "Celltype", values_to = "S.Score")
S.cellcycle_score_df_long$Celltype <- factor(S.cellcycle_score_df_long$Celltype, levels = c("Hepatocyte", "Endothelial", "Macrophage", "HSC", "Cholangio", "Basophil", "Neutr", "T", "NKT", "B", "NK", "Dendritic", "ILC"))
S.cellcycle_score_df_long$group <- factor(S.cellcycle_score_df_long$group, levels = c("D0","X8h", "D1", "D2", "D3", "D7"))
Celltype_col <- setNames(c("tomato1", '#4363d8', "dodgerblue4", "darkmagenta", "brown4", "hotpink" ,"burlywood1",  "cadetblue1", "cadetblue4", "chartreuse3", #T
                           
                           "chartreuse4", "slateblue1" ,"slateblue4"), c("Hepatocyte", "Endothelial", "Macrophage", "HSC", "Cholangio", "Basophil", "Neutr", "T", "NKT", "B", "NK", "Dendritic", "ILC"))
S.cellcycle_score_picture <- ggplot(S.cellcycle_score_df_long, aes(x = group, y = S.Score, group = Celltype, color = Celltype)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Group", y = "S.Score") +
  theme_minimal() +
  scale_color_manual(values =  as.vector(Celltype_col))+
  theme(text = element_text(size = 14), axis.title.x = element_blank())+
  facet_wrap(~Celltype, scales = "free_y", ncol = 3)
ggsave("PHX_All_celltype_S.Score_line_plot.pdf", S.cellcycle_score_picture, width = 9, height = 9)

G2M.cellcycle_score_df_raw_width_table <- read.csv("/data4/jinlianli/public_liver_hep/G2M.cellcycle_score_df_raw_width_table.txt", sep="")
rownames(G2M.cellcycle_score_df_raw_width_table) <- G2M.cellcycle_score_df_raw_width_table$Celltype
G2M.cellcycle_score_df_raw_width_table_scale <- G2M.cellcycle_score_df_raw_width_table[,2:7]
G2M.cellcycle_score_df_raw_width_table_scale_normalized <- apply(G2M.cellcycle_score_df_raw_width_table_scale, 1, function(x) x - x[1])
G2M.cellcycle_score_df_raw_width_table_scale_normalized <- as.data.frame(G2M.cellcycle_score_df_raw_width_table_scale_normalized)
G2M.cellcycle_score_df_raw_width_table_scale_normalized$group <- rownames(G2M.cellcycle_score_df_raw_width_table_scale_normalized)
G2M.cellcycle_score_df_long <- G2M.cellcycle_score_df_raw_width_table_scale_normalized %>%
  pivot_longer(cols = colnames(G2M.cellcycle_score_df_raw_width_table_scale_normalized)[1:13], names_to = "Celltype", values_to = "S.Score")
G2M.cellcycle_score_df_long$Celltype <- factor(G2M.cellcycle_score_df_long$Celltype, levels = c("Hepatocyte", "Endothelial", "Macrophage", "HSC", "Cholangio", "Basophil", "Neutr", "T", "NKT", "B", "NK", "Dendritic", "ILC"))
G2M.cellcycle_score_df_long$group <- factor(S.cellcycle_score_df_long$group, levels = c("D0","X8h", "D1", "D2", "D3", "D7"))
Celltype_col <- setNames(c("tomato1", '#4363d8', "dodgerblue4", "darkmagenta", "brown4", "hotpink" ,"burlywood1",  "cadetblue1", "cadetblue4", "chartreuse3", #T
                           
                           "chartreuse4", "slateblue1" ,"slateblue4"), c("Hepatocyte", "Endothelial", "Macrophage", "HSC", "Cholangio", "Basophil", "Neutr", "T", "NKT", "B", "NK", "Dendritic", "ILC"))
G2M.cellcycle_score_picture <- ggplot(G2M.cellcycle_score_df_long, aes(x = group, y = S.Score, group = Celltype, color = Celltype)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Group", y = "G2M.Score") +
  theme_minimal() +
  scale_color_manual(values =  as.vector(Celltype_col))+
  theme(text = element_text(size = 14), axis.title.x = element_blank())+
  facet_wrap(~Celltype, scales = "free_y", ncol = 3)
ggsave("PHX_All_celltype_G2M.Score_line_plot.pdf", G2M.cellcycle_score_picture, width = 9, height = 9)


# Figure-S1
# D
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
DimPlot(All_scrna, group.by = "Group", cols = colors)+ theme(legend.title = element_blank(), legend.text = element_text(size = 14)) + ggtitle("")
DimPlot(All_scatac, group.by = "Group", cols = colors)+ theme(legend.title = element_blank(), legend.text = element_text(size = 14)) + ggtitle("")

# E
SeuratObject <- readRDS("/data4/jinlianli/APAP/Data.analysis/Object/2.Merge/All.Rds")
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
celltype.group.table <- table(SeuratObject$Group, SeuratObject$celltype)
celltype.group.table <- as.data.frame(round(apply(celltype.group.table, 1, function(x) x/sum(x)), 3))
celltype.group.table$Total <- table(SeuratObject$celltype)
celltype.group.table_t <- as.data.frame(t(celltype.group.table))
celltype.group.table_t$Total <- c(table(SeuratObject$Group), sum(celltype.group.table$Total))
library(tidyverse)
library(RColorBrewer)
p <- SeuratObject@meta.data %>% 
  group_by(Group, celltype) %>% 
  summarise(count = n()) %>%
  ggplot() + geom_bar(aes(Group, 
                          count, 
                          fill=celltype), 
                      stat = "identity", 
                      position = "fill") + 
  labs(x="Group",y="Ratio")+ 
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"),
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain")) + 
  theme(plot.title = element_text(size = 15, face = "bold"), legend.title=element_text(size=15), legend.text=element_text(size=15)) + scale_fill_manual(values = colors)
ggsave("scRNA_Celltype.histogram.pdf", p, width = 8, height = 6)

All_scatac <- readRDS("/data4/jinlianli/APAP/Data.analysis/Object/scATAC.data.Rds")
colors <- c(rgb(247/255,147/255,30/255), rgb(248/255,115/255,106/255), rgb(169/255,169/255,169/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(192/255,193/255,48/255), rgb(157/255,115/255,194/255), rgb(183/255,76/255,171/255), rgb(230/255,134/255,201/255), rgb(140/255,198/255,63/255), rgb(255/255,191/255,15/255),
            rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(129/255,70/255,58/255), rgb(0/255,114/255,189/255), rgb(74/255,76/255,191/255))
celltype.group.table <- table(All_scatac$Group, All_scatac$celltype)
celltype.group.table <- as.data.frame(round(apply(celltype.group.table, 1, function(x) x/sum(x)), 3))
celltype.group.table$Total <- table(All_scatac$celltype)
celltype.group.table_t <- as.data.frame(t(celltype.group.table))
celltype.group.table_t$Total <- c(table(All_scatac$Group), sum(celltype.group.table$Total))
library(tidyverse)
library(RColorBrewer)
p <- All_scatac@meta.data %>% 
  group_by(Group, celltype) %>% 
  summarise(count = n()) %>%
  ggplot() + geom_bar(aes(Group, 
                          count, 
                          fill=celltype), 
                      stat = "identity", 
                      position = "fill") + 
  labs(x="Group",y="Ratio")+ 
  theme(axis.text.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 0.5, vjust = 0.5, face = "plain"),
        axis.title.x = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = 0.5, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain")) + 
  theme(plot.title = element_text(size = 15, face = "bold"), legend.title=element_text(size=15), legend.text=element_text(size=15)) + scale_fill_manual(values = colors)
ggsave("scATAC_Celltype.histogram.pdf", p, width = 4, height = 6)

# F
hsc_ec_hep_maf_Bas_8h_scRNA_avg <- read.table("hsc_ec_hep_maf_Bas_8h_scRNA_avg.txt")
ph <- pheatmap::pheatmap(hsc_ec_hep_maf_Bas_8h_scRNA_avg, scale = "row", 
                         cluster_rows = T, cluster_cols = T, show_rownames = F, 
                         # cutree_rows =4,
                         treeheight_row = 0,
                         fontsize = 14,
                         treeheight_col = 0,clustering_method = "average", cellwidth = 30,
                         border_color = "white")
ggsave("PHx_target_celltype_DEG_pheatmap3.pdf",ph,width = 8, height = 6)

hsc_ec_hep_maf_Bas_8h_DEG_GO_result <- read.csv("/data4/jinlianli/public_liver_hep/hsc_ec_hep_maf_Bas_8h_DEG_GO_result.txt", sep="")
set_terms1 <- c(
  # Hepatocyte
  "GO:0010038", #"response to metal ion"
  "GO:0006953", #"acute-phase response"
  # "GO:0072577", #"endothelial cell apoptotic process"
  "GO:0002526", #"acute inflammatory response"
  "GO:0050817", # coagulation
  # "GO:0042730", #"fibrinolysis"
  # "GO:0031639", #"plasminogen activation"
  # Basophil
  "GO:0010631", #"epithelial cell migration"
  "GO:0043620", #"regulation of DNA-templated transcription in response to stress"
  # "GO:0050678", #"regulation of epithelial cell proliferation"
  # "GO:0010586", #"miRNA metabolic process"
  "GO:0071559", #"response to transforming growth factor beta"
  # "GO:0035296", #"regulation of tube diameter"
  "GO:0030509", #"BMP signaling pathway"
  "GO:0032640", #"tumor necrosis factor production"
  # "GO:0019216", #regulation of lipid metabolic process
  # "GO:0009062", #fatty acid catabolic process
  # "GO:0006635", #fatty acid beta-oxidation
  # "GO:0006006", #glucose metabolic process
  # "GO:0030301", #cholesterol transport
  # Mac
  "GO:0140236", #translation at presynapse
  "GO:0046651", #lymphocyte proliferation
  "GO:0032943", #mononuclear cell proliferation
  "GO:0050863", #regulation of T cell activation
  "GO:0042113", #B cell activation
  # "GO:0007159", #leukocyte cell-cell adhesion
  # "GO:0035924", #cellular response to vascular endothelial growth factor stimulus
  # "GO:0042058", #"regulation of epidermal growth factor receptor signaling pathway"
  # "GO:0032728", #"positive regulation of interferon-beta production"
  "GO:0032623", #interleukin-2 production
  # "GO:0032661", #"regulation of interleukin-18 production" 
  # "GO:0036120", # "cellular response to platelet-derived growth factor stimulus"
  # "GO:0032674", #"regulation of interleukin-5 production"
  # "GO:0043123", #"positive regulation of canonical NF-kappaB signal transduction"
  "GO:0048012", #"hepatocyte growth factor receptor signaling pathway"
  # Endo
  # "GO:0006914", #autophagy
  # "GO:0001667", # ameboidal-type cell migration
  # "GO:0001935", #"endothelial cell proliferation"
  # "GO:0010631", #epithelial cell migration
  "GO:0033002", #muscle cell proliferation
  "GO:0010632", #regulation of epithelial cell migration
  # "GO:0035966", #response to topologically incorrect protein
  # "GO:0071559", #response to transforming growth factor beta
  # "GO:0198738", #cell-cell signaling by wnt
  "GO:0060070", #"canonical Wnt signaling pathway"
  # "GO:0042058", #regulation of epidermal growth factor receptor signaling pathwa
  "GO:0038084", #"vascular endothelial growth factor signaling pathway" 
  "GO:0031099", #regeneration
  "GO:0042246", #tissue regeneration
  "GO:0031100", #"animal organ regeneration"
  # "GO:0050678", #regulation of epithelial cell proliferation
  "GO:0072575", #epithelial cell proliferation involved in liver morphogenesis
  "GO:0003382" #epithelial cell morphogenesis
  # "GO:0072574", #hepatocyte proliferation
  # "GO:0032675", #regulation of interleukin-6 production
  #HSC
  # "GO:0001666", #response to hypoxia
  # "GO:0007263", # nitric oxide mediated signal transduction
  # "GO:0043114", #regulation of vascular permeability
  # "GO:0006979" #"response to oxidative stress"
)
result_df_sig <- hsc_ec_hep_maf_Bas_8h_DEG_GO_result
term_ids <- c()
for (i in seq_along(set_terms1)) {
  id <- grep(set_terms1[i], result_df_sig$ID)
  term_ids <- c(term_ids, id)
}
hsc_hep_endo_maf_6h_set_terms_df <- result_df_sig[term_ids,]
hsc_hep_endo_maf_6h_set_terms_df <- hsc_hep_endo_maf_6h_set_terms_df[,c("Description","p.adjust","cluster")]
hsc_hep_endo_maf_6h_set_terms_df$log10_p.adjust <- -log10(hsc_hep_endo_maf_6h_set_terms_df$p.adjust)
hsc_hep_endo_maf_6h_set_terms_df <- hsc_hep_endo_maf_6h_set_terms_df[, c("cluster", "log10_p.adjust", "Description")]

df_wide <- hsc_hep_endo_maf_6h_set_terms_df %>%
  pivot_wider(names_from = cluster, values_from = c(log10_p.adjust))
rownames(df_wide) <- df_wide$Description
rn <- df_wide$Description
head(df_wide)
df_wide <- as.data.frame(df_wide)
df_wide <- df_wide[, -1]
rownames(df_wide) <- rn
df_wide <- df_wide[, c("Hepatocyte", "Basophil", "Macrophage", "Endothelial", "HSC")]
ph3 <- pheatmap::pheatmap(df_wide, border_color = "white", fontsize = 12, cluster_rows = F, cluster_cols = F, cellwidth = 30, cellheight = 17)
ggsave("hsc_ec_hep_maf_Bas_8h_DEG_GO_select_terms_pheatmap3.pdf", ph3, width = 8, height = 8)


# Figure2
# A
library(CellChat)
setwd("/data2/lijinlian/APAP_project/cellchat/allcelltype")
groups_vector <- c("0h", "6h", "12h", "24h", "48h", "96h")
list_cellchat_object_apap <- list()
for (j in 1:length(groups_vector)) {
  nfile <- paste(groups_vector[j], "cellchat_object.Rds", sep = "_")
  cellchat_obj <- readRDS(nfile)
  cellchat_obj2 <- aggregateNet(cellchat_obj, thresh = 0.05)
  list_cellchat_object_apap[[groups_vector[j]]] <- cellchat_obj2
}
list_cellchat_object_apap <- lapply(list_cellchat_object_apap, function(x){
  cellchat_obj <- netAnalysis_computeCentrality(x, slot.name = "netP")
})
cellchat <- mergeCellChat(list_cellchat_object, add.names = names(list_cellchat_object))
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
pdf("cellchat_allcelltype_weight_to_hep_pheatmap.pdf", width = 4, height = 4)
pheatmap::pheatmap(w_to_hep_l, cluster_rows = F, cluster_cols = F, border_color = "white", fontsize = 14)
dev.off()

# B
data_GO_sim <- read.table("/data2/lijinlian/APAP_project/cellchat/allcelltype_2/cellchat_to_hep2/Endo_Incoming_to_hepatocytes_data_GO_sim.txt")
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
pgo <- ggplot(select_Go_receptors_terms,
              aes(x=Description,y=-log10(p.adjust), fill=Cluster)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = c("orange", "purple")) +
  coord_flip() +
  xlab("GO terms") +
  ylab("-log10(p.adjust)") +
  # labs(title = "")+  #设置标题
  theme_bw()+theme(text = element_text(size = 14), axis.text = element_text(size = 16, colour = "black"))
ggsave("GO_barplot_Endo_Incoming_to_hepatocytes_selected_terms.pdf", pgo, width = 7.6, height = 5)

# C
p <- jjDotPlot(object = Hep_scRNA, col.min = 0, base_size = 12, textSize = 14, segWidth = 0.3, 
               bar.width = 2.5,
               gene = c("Egfr", "Lrp1", "Cd81", "Ptprf", "Jak1", "Stat3", "Ccnd1",  "Cdk1", "Mki67", "Stmn1", "Top2a"), 
               xtree = F, ytree = F, legend.position = "top",
               plot.margin = c(0,1,0,0),
               id = 'Group')
ggsave("Hep_scRNA_regeneration_receptor_Group_jjDotPlot4.pdf", p, width =4, height = 4)

# D
W_df <- read.table("/data2/lijinlian/APAP_project/cellchat/human2024_allcelltype_cell10/HUMAN_cellchat_allcelltype_weight_to_hep_weight.txt")
pdf("HUMAN_cellchat_allcelltype_weight_to_hep_pheatmap.pdf", width = 4, height = 3)
pheatmap::pheatmap(W_df, 
                   cluster_rows = F, cluster_cols = F, 
                   border_color = "white", main = "Interaction score to hepatocytes",
                   fontsize = 14, cellheight = 25)
dev.off()

# E
cci_data_mm_zf <- read.table("/data2/lijinlian/APAP_project/cellchat/human2024_allcelltype_cell10/APAP_cci_conserved_Weights_data_mm_hs.txt")
cci_data_mm_zf_to_hep <- subset(cci_data_mm_zf, rep_target == "Hepatocyte/Hepatocytes")
cci_data_mm_zf_to_hep <- subset(cci_data_mm_zf_to_hep, rep_source != "Hepatocyte/Hepatocytes")
ChordDiagram(net = cci_data_mm_zf_to_hep, 
             filename = paste("mm", "hs", "_Conserved_interaction_score_to_Hepatocyte_between_mouse_and_human_APAP_chordDiagram.pdf.pdf", sep = "_"), 
             brewer_pal_used = "Dark2", link_colors_threshold = 0.9, Score_factor_size = 150, picture_width = 6, picture_height = 6)

# F
conserved_result_mm_zf <- readRDS("/data2/lijinlian/APAP_project/cellchat/human2024_allcelltype_cell10/APAP_cci_conserved_result_mm_hs.Rds")
cci_conserved_result_in_hs <- conserved_result_mm_zf[["sp2_cci_conserved_df"]]
cci_conserved_result_in_hs_endo_to_hep <- subset(cci_conserved_result_in_hs, c_pair == "hsEndothelia_hsHepatocytes")
cci_conserved_result_in_hs_endo_to_hep_chordD <- cci_conserved_result_in_hs_endo_to_hep[, c("ligand","receptor", "prob")]
cci_conserved_regeneration_receptors_result_in_hs_endo_to_hep_chordD <- subset(cci_conserved_result_in_hs_endo_to_hep_chordD, receptor %in% c("CD81","EGFR","LRP1","PTPRF"))
pdf("Conserved_Signalling_of_regeneration_receptors_from_Endothelial_to_Hepatocyte_human_APAP_chordDiagram.pdf", width = 7, height = 7)
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
dev.off()
circos.clear()


# Figure-S2
# A
APAP_merge <- readRDS("/data4/jinlianli/Spatial_apap/APAP_merge_241212.Rds")
genes <- c("Egfr", "Hbegf")
Hmox1_Hbegf_interaction_all_cell <- as.data.frame(t(as.matrix(APAP_merge@assays$SCT@counts[genes,])))
Hmox1_Hbegf_interaction_all_cell$Hmox1_Hbegf <- Hmox1_Hbegf_interaction_all_cell$Egfr*Hmox1_Hbegf_interaction_all_cell$Hbegf
APAP_merge$Hbegf_Egfr <- Hmox1_Hbegf_interaction_all_cell$Hmox1_Hbegf
pdf("APAP_merge_Hbegf_Egfr_interaction_score_SpatialDimPlot_250408.pdf", width = 12, height = 4)
SpatialFeaturePlot(APAP_merge, features = c("Hbegf_Egfr"), pt.size.factor = 2.5)
dev.off()

# B
w_to_hep_l <- read.table("/home/lijinlian/test/PHX_scRNA_miguel/cellchat_2/PHx_celltype_to_hepatocyte_interaction_weight_table.txt")
ph2 <- pheatmap::pheatmap(w_to_hep_l, cluster_rows = F, cluster_cols = F, border_color = "white", fontsize = 14)
ggsave("PHx_celltype_to_hepatocyte_interaction_weight_pheatmap2.pdf",ph2,width = 4.5, height = 3.5)

# C
Egfr_endo_to_hep_signalling <- read.table("/data4/jinlianli/public_liver_hep/PHx_8h_Signalling_of_Egfr_from_Endothelial_to_Hepatocyte_chordDiagram_data.txt")
library(tidyverse)
library(circlize)
circos.par(start.degree = 90)
pdf("PHx_8h_Signalling_of_Egfr_from_Endothelial_to_Hepatocyte_chordDiagram.pdf", width = 8, height = 8)
chordDiagram(
  Egfr_endo_to_hep_signalling, 
  directional = 1,
  annotationTrack = "grid", link.arr.width = 0.1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  big.gap = 8, diffHeight = -0.01,
  small.gap = 3,
  preAllocateTracks = list(
    track.height = 0.25)
)
circos.track(
  track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$ylim[1], 
      CELL_META$sector.index,  facing = "clockwise", 
      niceFacing = TRUE, adj = c(0, 0.3), cex = 1.25,
    )
  }, bg.border = NA
)
title("Signalling of regeneration receptors from Endothelial to Hepatocyte")
dev.off()
circos.clear()

# D
# 该图直接从数据库获取https://db.cngb.org/stomics/lista/


# Figure3
# A
Hep_scRNA <- readRDS("/data2/lijinlian/APAP_project/Hep_scRNA.Rds")
colors_set_cluster <- c("red", "orange", "purple", "blue")
names(colors_set_cluster) <- levels(Hep_scRNA$anno1)
dimp1 <- DimPlot(Hep_scRNA, group.by = "anno1", cols = colors_set_cluster, pt.size = 0.5) + tidydr::theme_dr() + theme(title = element_blank(), text = element_text(size = 14))
ggsave("Hep_scRNA_clusters_anno1_Dimplot.pdf", dimp1, width = 6, height = 4)

# B
library(scales)
library(RColorBrewer)
pdf("mouse_Hep_scRNA_regeneration_receptors1.pdf", width = 4, height = 4)
FeaturePlot(Hep_scRNA, features = "regeneration_receptors1", pt.size = 1.2, order = T) &
  scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),
                        # values = rescale(c(min(Hepatocyte$Pericentral.signature1),0,max(Hepatocyte$Pericentral.signature1))),
                        guide = "colorbar",
                        limits = c(min(Hep_scRNA$regeneration_receptors1), max(Hep_scRNA$regeneration_receptors1))) &
  tidydr::theme_dr()
# theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) + theme(legend.position = "right") 
dev.off()

# C
Roe_Analysis <- function(train, tissue = "tissue", celltype = "finalannotate", plot_location = ".", data_location = ".", perfix = 'rna') {
  library(epitools)
  library(reshape2)
  library(tidyverse)
  data <- table(train@meta.data[, tissue], train@meta.data[, celltype]) %>%
    matrix(nrow = dim(table(train@meta.data[, tissue])))
  ei_data <- expected(data)
  data_roe <- data / ei_data
  colnames(data_roe) <- names(table(train@meta.data[, celltype]))
  rownames(data_roe) <- names(table(train@meta.data[, tissue]))
  data_roe_plot <- melt(data_roe)
  data_roe_plot$significant <- case_when(
    data_roe_plot$value < 1 ~ "±",
    data_roe_plot$value > 1 & data_roe_plot$value < 1.5 ~ "+",
    data_roe_plot$value > 1.5 & data_roe_plot$value < 2 ~ "++",
    data_roe_plot$value > 2 ~ "+++",
  )
  colnames(data_roe_plot)[colnames(data_roe_plot) == "Var1"] <- "tissue"
  colnames(data_roe_plot)[colnames(data_roe_plot) == "Var2"] <- "celltype"
  p <- ggplot(data_roe_plot, aes(x = celltype, y = tissue, label = significant)) +
    geom_tile(aes(fill = value), color="gray",size=0.5) +
    coord_equal() + 
    geom_text(color = "black") + 
    scale_fill_viridis_c(option = "C", direction = 1, alpha = 1) +
    theme(axis.title = element_blank(),legend.title = element_blank(), legend.position = "right", axis.text.x = element_text(angle = 35, size = 14, colour = "black", hjust = 1, vjust = 1), axis.text.y = element_text(size = 14, colour = "black"))
  pdf(paste0(plot_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".pdf"), width = 4, height = 4.5)
  print(p)
  dev.off()
  write.table(data_roe, paste0(data_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".txt"), sep = "\t", quote = FALSE)
}
Roe_Analysis(Hep_scRNA, tissue = "Group", celltype = "anno1", perfix = 'Hep_scRNA_anno1')

# D
Mfsd2a_Igfbp2_exp <- as.data.frame(t(as.data.frame(Hep_scRNA@assays$RNA@data[c("Igfbp2", "Mfsd2a"),])))
Mfsd2a_Igfbp2_exp <- as.data.frame(Mfsd2a_Igfbp2_exp)
Mfsd2a_Igfbp2_exp$Group <- as.data.frame(Hep_scRNA$Group)[,1]
Mfsd2a_Igfbp2_exp$coexp <- Mfsd2a_Igfbp2_exp$Igfbp2*Mfsd2a_Igfbp2_exp$Mfsd2a
p2 <- ggplot(Mfsd2a_Igfbp2_exp, aes(x = Igfbp2, y = Mfsd2a, color = coexp)) +
  geom_point(alpha = 0.6, size = 2) +  
  scale_color_viridis_c()+
  labs(title = "Correlation of Igfbp2 and Mfsd2a",
       x = "Igfbp2",
       y = "Mfsd2a") +
  theme_minimal() + 
  facet_grid(.~Group) +
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 12,colour = 'black',face = 'plain'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=12,hjust=0.5),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"))
ggsave("Mfsd2a_Igfbp2_exp_scatterplot_split3_sct.pdf", p2, width = 24, height = 4)

# E
Hep_scATAC <- readRDS("/data2/lijinlian/APAP_project/Hep_scATAC.Rds")
gcolors <- ggsci::pal_igv()(6)
names(gcolors) <- levels(Hep_scATAC$Clusters)
dimp1 <- DimPlot(Hep_scATAC, group.by = "Clusters", cols = gcolors, pt.size = 1) + tidydr::theme_dr()
ggsave("Hep_scATAC_clusters_Dimplot.pdf", dimp1, width = 4, height = 4)

Hep_scATAC_anno1_markers <- c("Igfbp2", "Mfsd2a", "Cyp2e1","Cdh1","Mt1", "Kdr", "Hsph1")
DefaultAssay(Hep_scATAC) <- "ACTIVITY"
p <- jjDotPlot(object = Hep_scATAC, col.min = 0,base_size = 14, textSize = 14, segWidth = 0.3, bar.width = 2.5,
               gene = Hep_scATAC_anno1_markers, xtree = F, ytree = F, legend.position = "top",plot.margin = c(0,1,0,0),
               id = 'anno1')
ggsave("Hep_scATAC_anno1_markers_anno1_jjDotPlot2.pdf", p, width =4, height = 4.5)

# F
pdf("mouse_Hep_scATAC_regeneration_receptors1.pdf", width = 4, height = 4)
FeaturePlot(Hep_scATAC, features = "regeneration_receptors1", pt.size = 1.2, order = T) &
  scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),
                        # values = rescale(c(min(Hepatocyte$Pericentral.signature1),0,max(Hepatocyte$Pericentral.signature1))),
                        guide = "colorbar",
                        limits = c(min(Hep_scATAC$regeneration_receptors1), max(Hep_scATAC$regeneration_receptors1))) &
  tidydr::theme_dr()
# theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) + theme(legend.position = "right") 
dev.off()


# Figure-S3
# A
colors_set <- c("green4", "red2", "orange2", "purple", "pink2", "springgreen", "blue","sienna4")
dimp2 <- DimPlot(Hep_scRNA, group.by = "Group", cols = unname(colors_set), pt.size = 0.5) + tidydr::theme_dr()+ theme(title = element_blank(), text = element_text(size = 14))
ggsave("Hep_scRNA_clusters_Group_Dimplot.pdf", dimp2, width = 4, height = 4, limitsize = F)
Hep_scRNA_anno1_markers <- c("Thrsp", "Igfbp2", "Cyp2e1", "Cdk1")
p <- jjDotPlot(object = Hep_scRNA, col.min = 0,base_size = 14, textSize = 14, segWidth = 0.4, bar.width = 3,
               gene = unique(Hep_scRNA_anno1_markers), xtree = F, ytree = F, legend.position = "top",plot.margin = c(0,1,0,0),
               id = 'anno1')
ggsave("Hep_scRNA_anno1_markers_jjDotPlot.pdf", p, width =4, height = 4)

# B
data_GO_sim <- readRDS("/data2/lijinlian/APAP_project/Hep_scRNA/Hep_4/Hep_scRNA_cluster_data_GO_sim.Rds")
data_GO_sim_Hep_scrna <- data_GO_sim@compareClusterResult
terms_list <- list(data_GO_sim_Hep_scrna$Description[c(5:10,18,19,21,22, 
                                                       38,42,43,46,47,48,52,55,81,83,
                                                       275:279,285, 286, 289, 
                                                       465:468, 471,472, 475:478
)])
dp <- dotplot(data_GO_sim, showCategory=terms_list[[1]]) + scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=65)) & theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12))
ggsave(paste0(object.name, "_dotplot3.pdf"), dp, width = 7, height = 9, limitsize = F)

# C
pdf("Igfbp2_Hep_scRNA_split_FeaturePlot.pdf", width = 24, height = 4)
FeaturePlot(
  Hep_scRNA,
  features = c("Igfbp2"),
  max.cutoff = 3, 
  order = T, split.by = "Group",
  pt.size = 1.2  
) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
dev.off()
pdf("Mfsd2a_Hep_scRNA_split_FeaturePlot.pdf", width = 24, height = 4)
FeaturePlot(
  Hep_scRNA,
  features = c("Mfsd2a"),
  max.cutoff = 3, 
  order = T, split.by = "Group",
  pt.size = 1.2  
) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
dev.off()

# D
DefaultAssay(Hep_scATAC) <- "ACTIVITY"
pdf("Hep_scATAC_Igfbp2_plot_density.pdf",width = 4, height = 4)
d1 <- Nebulosa::plot_density(Hep_scATAC, 
                             features = c("Igfbp2"),
                             size = 1.2, 
                             # pal = "magma",
                             joint=TRUE, slot = "data",combine = T, adjust = 1,
                             reduction = "umap")& tidydr::theme_dr()
d1
dev.off()
pdf("Hep_scATAC_Mfsd2a_plot_density.pdf",width = 4, height = 4)
d1 <- Nebulosa::plot_density(Hep_scATAC, 
                             features = c("Mfsd2a"),
                             size = 1.2, 
                             # pal = "magma",
                             joint=TRUE, slot = "data",combine = T, adjust = 1,
                             reduction = "umap")& tidydr::theme_dr()
d1
dev.off()

# E
snRNA.apap2023.mouse.Hepatocyte <- readRDS("/data4/jinlianli/APAP/ANALYSIS/Nature2024_mouse_APAP_scRNA_Hep/snRNA.apap2023.mouse.Hepatocyte.Rds")
colors_set <- ggsci::pal_lancet()(9)
pdf("Mouse_Hep_2024scRNA_seurat_clusters_umap.pdf", width = 5, height = 5)
DimPlot(snRNA.apap2024.mouse.Hepatocyte,cols =  cb_palette, group.by = "timepoint", pt.size = 0.5) + tidydr::theme_dr() 
DimPlot(snRNA.apap2024.mouse.Hepatocyte, cols = colors_set, group.by = "seurat_clusters", pt.size = 0.5) + tidydr::theme_dr()
dev.off()

regeneration_receptors <- c("Egfr", "Lrp1", "Cd81", "Ptprf")
snRNA.apap2024.mouse.Hepatocyte <- AddModuleScore(snRNA.apap2024.mouse.Hepatocyte, features = list(regeneration_receptors), name = "regeneration_receptors")
library(scales)
library(RColorBrewer)
pdf("mouse.Hep_regeneration_receptors_signature.pdf", width = 4, height =4)
FeaturePlot(snRNA.apap2024.mouse.Hepatocyte, features = "regeneration_receptors1", pt.size = 1.2, order = T) &
  scale_color_gradientn(colours =  rev(brewer.pal(n = 11, name = "RdBu")), 
                        # values = rescale(c(min(Hepatocyte$Pericentral.signature1),0,max(Hepatocyte$Pericentral.signature1))),
                        guide = "colorbar",
                        limits = c(min(snRNA.apap2024.mouse.Hepatocyte$regeneration_receptors1), max(snRNA.apap2024.mouse.Hepatocyte$regeneration_receptors1))) &
  tidydr::theme_dr()
# theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) + theme(legend.position = "right") 
dev.off()

# F
Hepatocyte_Hypoxia <- c("MFSD2A","IGFBP2")
f.names <- "Human.Hep_IGFBP2_MFSD2A"
ncol = 2
fp <- FeaturePlot(human_hep, features = Hepatocyte_Hypoxia, ncol = ncol, order = T, pt.size = 1.2) & viridis::scale_colour_viridis(option = "plasma") & theme(axis.title.x = element_blank()) & tidydr::theme_dr()
ggsave(paste(f.names, "ncol_FeaturePlot.pdf", sep = "_"), fp, width = ncol*4, height = ceiling(length(unique(Hepatocyte_Hypoxia))/ncol)*4, limitsize = F)

# G
Hep_human_snRNA_apap <- readRDS("/data4/jinlianli/APAP/ANALYSIS/Nature2024_human_APAP_scRNA_Hep/Hep_human_snRNA_apap.Rds")
colors_set <- ggsci::pal_lancet()(9)
pdf("Human_Hep_2024scRNA_seurat_clusters_umap.pdf", width = 5, height = 5)
DimPlot(Hep_human_snRNA_apap, cols = colors_set, group.by = "condition", pt.size = 0.5) + tidydr::theme_dr() 
DimPlot(Hep_human_snRNA_apap, cols = colors_set, group.by = "seurat_clusters", pt.size = 0.5) + tidydr::theme_dr()
dev.off()

regeneration_receptors <- c("EGFR",  "LRP1", "CD81", "PTPRF")
human_hep <- AddModuleScore(human_hep, features = list(regeneration_receptors), name = "regeneration_receptors")
library(scales)
library(RColorBrewer)
pdf("human_hep_regeneration_receptors_signature.pdf", width = 4, height =4)
FeaturePlot(human_hep, features = "regeneration_receptors1", pt.size = 1.2, order = T) &
  scale_color_gradientn(colours =  rev(brewer.pal(n = 11, name = "RdBu")), 
                        # values = rescale(c(min(Hepatocyte$Pericentral.signature1),0,max(Hepatocyte$Pericentral.signature1))),
                        guide = "colorbar",
                        limits = c(min(human_hep$regeneration_receptors1), max(human_hep$regeneration_receptors1))) &
  tidydr::theme_dr()
# theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid")) + theme(legend.position = "right") 
dev.off()

# H
Hepatocyte_Hypoxia <- c("MFSD2A","IGFBP2")
f.names <- "Human.Hep_IGFBP2_MFSD2A"
ncol = 2
fp <- FeaturePlot(human_hep, features = Hepatocyte_Hypoxia, ncol = ncol, order = T, pt.size = 1.2) & viridis::scale_colour_viridis(option = "plasma") & theme(axis.title.x = element_blank()) & tidydr::theme_dr()
ggsave(paste(f.names, "ncol_FeaturePlot.pdf", sep = "_"), fp, width = ncol*4, height = ceiling(length(unique(Hepatocyte_Hypoxia))/ncol)*4, limitsize = F)


# Figure4
# A
cds <- readRDS("/data2/lijinlian/APAP_project/Hep_scRNA/monocle3_wjl/cds_obj.Rds")
root.name <- "Hep_3"
pt1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
                  label_leaves = FALSE, trajectory_graph_color = "purple", cell_size = 0.8,
                  label_branch_points = FALSE, 
                  label_roots = F) + 
  theme(legend.position="right", 
        # legend.key.size =  unit(0.2, "inches"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 16), 
        text = element_text(size = 16), axis.text = element_text(size = 16)) + 
  # theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))  + 
  tidydr::theme_dr()+
  scale_color_gradientn(colours = c("#0088C3","#6FBCDB","#BCE5F1","#F2EFDA","#FFC482","#FF7B4C","#FF5831") # c("blue4", "#2b7ab5", "#a9d0e0","#b2e091", "#fa8d8f", "red3")
  )
ggsave(paste0("Root.", root.name, ".Trajectory3.pdf"), pt1, width = 5, height = 5)

cds_sub33 <- readRDS("/data2/lijinlian/APAP_project/Hep_scRNA/monocle3_wjl/cds_choose_graph_segments.Rds")
id_s2 <- as.vector(sapply(rownames(colData(cds_sub33)), function(x) grep(x,rownames(colData(cds)))))
cds@colData$seurat_clusters <- Hep_scRNA$anno1
Trajectory2 <- as.character(cds@colData$seurat_clusters)
Trajectory2[setdiff(1:16017, id_s2)] <- "not_select"
cds@colData$Trajectory2 <- Trajectory2
sce_sub33 <- Hep_scRNA[, rownames(colData(cds_sub33))]
sce_sub33$pseudotime <- cds_sub33@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pseudotime_df1 <- as.data.frame(sce_sub33$pseudotime)
colnames(pseudotime_df1) <- "pseudotime"
pseudotime_df1$annotation <- sce_sub33$anno1
pseudotime_df1$Group <- sce_sub33$Group
colors_set_cluster <- c("red", "orange", "purple", "blue")
names(colors_set_cluster) <- levels(Hep_scRNA$anno1)
cds@colData$Trajectory2 <- factor(cds@colData$Trajectory2, levels = c(levels(pseudotime_df1$annotation), "not_select"))
root.name <- "Trajectory2"
col_set2 <- c(colors_set_cluster[levels(pseudotime_df1$annotation)],"white")
names(col_set2)[length(col_set2)] <- "not_select"
pt11 <- plot_cells(cds, color_cells_by = "Trajectory2", label_cell_groups = FALSE, 
                   label_leaves = FALSE, trajectory_graph_color = "purple", cell_size = 0.8,
                   label_branch_points = FALSE, group_cells_by = "Trajectory2",
                   label_roots = F) + 
  theme(legend.position="right", legend.key.size =  unit(0.2, "inches"), legend.title = element_blank(), legend.text = element_text(size = 16), text = element_text(size = 16), axis.text = element_text(size = 16)) + 
  # theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"))  + 
  tidydr::theme_dr() + scale_color_manual(values = col_set2)
ggsave(paste0("Root.", root.name, "_select_colors_annotaion.pdf"), pt11, width = 5.2, height = 5)

# B
ppsed3 <- ggplot(pseudotime_df1, aes(x = pseudotime, y = Group, fill = Group)) +#数据
  geom_density_ridges(stat = "density_ridges", from = 0) +#绘制山脊图
  theme(legend.position = "right", #图例去除
        axis.text.y = element_blank(),
        # axis.line.y = element_blank(),#去除Y轴轴线
        axis.ticks.y = element_blank()
  )+ #去除Y轴刻度
  # scale_y_discrete(expand = c(0, 0.5)) + #调整起始图形距离X轴距离
  scale_fill_manual(values = ggsci::pal_nejm()(6)) +
  theme_bw()+
  theme(panel.border = element_rect(fill=NA,color= "black", size=1, linetype= "solid"), text = element_text(size = 14, colour = "black"), legend.title = element_blank(), axis.text = element_text(size = 14, colour = "black")) #+
ggsave("geom_density_ridges_order_trajectory_Group_0.pdf", ppsed3, width = 6, height = 5)

# C
pt.matrix2 <- read.table("/data2/lijinlian/APAP_project/cellchat/allcelltype_2/cellchat_to_hep2/Cluster3_pseudotime_gene_pheatmap_table.txt")
library(circlize)
Kmeans_clustering_ENS_re <- read.table("/data2/lijinlian/APAP_project/cellchat/allcelltype_2/cellchat_to_hep2/Kmeans_clustering_ENS_3_Clusters.txt")
cols <- c("#ED0000FF","#42B540FF" , "#00468BFF")
row_anno <- data.frame(Cluster = Kmeans_clustering_ENS_re$KmeansGroup, row.names = Kmeans_clustering_ENS_re$Symbol)
ann_colors = list(Cluster = setNames(cols, paste0("C", 1:3)))
pdf("trajectory_genes_heatmap.pdf")
pheatmap::pheatmap(pt.matrix2, annotation_legend = F,
                           cluster_cols = F, cluster_rows = F, show_rownames = F,
                           clustering_method = "ward.D", show_colnames = F,
                           color = colorRampPalette(c("blue4","#33FFFf","yellow", "red3"))(100),
                           annotation_row = row_anno,
                           annotation_colors = ann_colors,
                           main = "Pseudotime")
dev.off()

# D
pt.matrix3 <- readRDS("/data2/lijinlian/APAP_project/Hep_scRNA/monocle3_wjl/IReNA_bulkatac_3M/Module_Tfs.heatmap.data_169tfs.Rds")
cols <- c("#ED0000FF","#42B540FF" , "#00468BFF")
row_anno <- data.frame(Cluster = Kmeans_clustering_ENS_re$KmeansGroup[match(rownames(pt.matrix3), Kmeans_clustering_ENS_re$Symbol)], row.names = Kmeans_clustering_ENS_re$Symbol[match(rownames(pt.matrix3), Kmeans_clustering_ENS_re$Symbol)])
ann_colors = list(Cluster = setNames(cols, paste0("C", 1:3)))
pdf("TFs_pseudotime_pheatmap_annorow2.pdf", width = 8, height = 8)
pheatmap::pheatmap(pt.matrix3, 
                           cluster_cols = F, cluster_rows = F, show_rownames = F,
                           clustering_method = "ward.D", 
                           color = colorRampPalette(c("navy","white", "red3"))(100),
                           annotation_row = row_anno, annotation_colors = ann_colors,
                           main = "")
dev.off()

# E
# 5T移动硬盘里面，原始pdf图片的内皮细胞对应的文件夹里面Hep_irena_bulk_enrich_tfs_regulate_module.cys是调控网络画图文件，可以自己打开修改

# F
library(readxl)
gene_exp_data <- read_excel("/data4/jinlianli/GSE210842_WTliver_vs_HNF4AKOliver.gene_exp.diff.xlsx")
gene_exp_data_sig <- gene_exp_data[gene_exp_data$significant == "yes",]
gene_exp_data_sig <- as.data.frame(gene_exp_data_sig)
rownames(gene_exp_data_sig) <- gene_exp_data_sig$gene
interest_genes <- c("Hnf4a", "Egfr")
gene_exp_data_sig_interest <- gene_exp_data_sig[match(interest_genes, gene_exp_data_sig$gene),]
gene_exp_data_sig_interest_df <- gene_exp_data_sig_interest[,c("WTliver", "HNF4AKOliver")]

pdf("Hnf4a_ko_target_exp_heatmap.pdf", width = 3, height = 3)
pheatmap::pheatmap(gene_exp_data_sig_interest_df, cluster_rows = F,cluster_cols = F, show_rownames = T, fontsize = 16,show_colnames = T)
dev.off()


# Figure5
# A
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
cols <- ggsci::pal_ucscgb()(8)
cols[3] <- "blue"
names(cols) <- levels(Endo_scRNA1$Clusters)
dimp1 <- DimPlot(Endo_scRNA1, group.by = "Clusters", cols = cols, pt.size = 0.8) + tidydr::theme_dr() + theme(text = element_text(size = 16), title = element_blank())
ggsave("Endo_scRNA1_clusters_Dimplot.pdf", dimp1, width = 6.4, height = 4.5)

Endo_scRNA1_clusters_marker_genes <- c("Vwf", "Stab2", "Rspo3", "Efnb1", "Chst2", "Apold1", "Alb", "Rpl32", "Mki67")
p <- jjDotPlot(object = Endo_scRNA1, col.min = 0,base_size = 14, textSize = 14, bar.width = 3, segWidth = 0.3,
               gene = Endo_scRNA1_clusters_marker_genes, xtree = F, ytree = F, legend.position = "bottom",
               id = 'Clusters')
ggsave("Endo_scRNA1_clusters_markers_jjDotPlot.pdf", p, width = 6, height = 5)

# B
Roe_Analysis <- function(train, tissue = "tissue", celltype = "finalannotate", plot_location = ".", data_location = ".", perfix = 'rna') {
  library(epitools)
  library(reshape2)
  library(tidyverse)
  data <- table(train@meta.data[, tissue], train@meta.data[, celltype]) %>%
    matrix(nrow = dim(table(train@meta.data[, tissue])))
  ei_data <- expected(data)
  data_roe <- data / ei_data
  colnames(data_roe) <- names(table(train@meta.data[, celltype]))
  rownames(data_roe) <- names(table(train@meta.data[, tissue]))
  data_roe_plot <- melt(data_roe)
  data_roe_plot$significant <- case_when(
    data_roe_plot$value < 1 ~ "±",
    data_roe_plot$value > 1 & data_roe_plot$value < 1.5 ~ "+",
    data_roe_plot$value > 1.5 & data_roe_plot$value < 2 ~ "++",
    data_roe_plot$value > 2 ~ "+++",
  )
  colnames(data_roe_plot)[colnames(data_roe_plot) == "Var1"] <- "tissue"
  colnames(data_roe_plot)[colnames(data_roe_plot) == "Var2"] <- "celltype"
  p <- ggplot(data_roe_plot, aes(x = tissue, y = celltype, label = significant)) +
    geom_tile(aes(fill = value), color="gray",size=0.5) +
    coord_equal() + 
    geom_text(color = "black") + 
    scale_fill_viridis_c(option = "C", direction = 1, alpha = 1) +
    theme(axis.title = element_blank(), text = element_text(size = 14, colour = "black"), legend.title = element_blank())
  pdf(paste0(plot_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".pdf"), width = 5, height = 6)
  print(p)
  dev.off()
  write.table(data_roe, paste0(data_location, "/data_roe_", tissue, "_", celltype, "_", perfix, ".txt"), sep = "\t", quote = FALSE)
}
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
Roe_Analysis(Endo_scRNA1, tissue = "Group", celltype = "Clusters", perfix = 'Endo_scRNA')

# C
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
cluster2_genes <- c("Hbegf")
f.names <- "Endo_scRNA_Hbegf"
ncol = 1
fp <- FeaturePlot(Endo_scRNA1, features = cluster2_genes, ncol = ncol, order = T, pt.size = 1.2) & viridis::scale_colour_viridis(option = "plasma") & theme(axis.title.x = element_blank()) & tidydr::theme_dr()
ggsave(paste(f.names, "FeaturePlot.pdf", sep = "_"), fp, width = ncol*4, height = ceiling(length(unique(cluster2_genes))/ncol)*4, limitsize = F)

# D
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
stress_genes_ec_scrna <- c("Atf3", "Atf4", "Atf6","Egr1", "Fosl1","Hif1a", "Hsp90aa1", "Hspa1a", "Hspa1b", "Hspa8","Hspa9", "Hspb1", "Hspd1", "Hsph1", "Bcl10", "Bcl2l1", "Cdkn1a", "Hmox1", "Tlr4")
p <- jjDotPlot(object = Endo_scRNA1, col.min = 0,
          base_size = 14, textSize = 16, 
          bar.width = 3, segWidth = 0.3,
          gene = stress_genes_ec_scrna, 
          xtree = F, ytree = F, legend.position = "bottom",
          id = 'Clusters')+ ggtitle("Stress signatures")
ggsave("stress_genes_ec_scrna_jjDotPlot3.pdf", p, width =6, height = 5)

# E
Endo_scATAC1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scATAC1.Rds")
g_cols_endo_scatac <- c(ggsci::pal_futurama()(7)[1:4], c("green", "blue", "pink"))
p11 <- DimPlot(Endo_scATAC1, cols = g_cols_endo_scatac, group.by = "Clusters", pt.size = 1.2) + tidydr::theme_dr()+ theme(legend.title = element_blank(), legend.text = element_text(size = 14), text = element_text(size = 14))
ggsave("Endo_scATAC_UMAP_annotation.pdf", p11, width = 6.4, height = 5, limitsize = F)

# F
Endo_scATAC1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scATAC1.Rds")
stress_genes_ec_scrna <- c("Atf3", "Atf4", "Atf6","Egr1", "Fosl1","Hif1a", "Hsp90aa1", "Hspa1a", "Hspa1b", "Hspa8","Hspa9", "Hspb1", "Hspd1", "Hsph1", "Bcl10", "Bcl2l1", "Cdkn1a", "Hmox1", "Tlr4")
p <- jjDotPlot(object = Endo_scRNA1, col.min = 0,
          base_size = 14, textSize = 16, 
          bar.width = 3, segWidth = 0.3,
          gene = stress_genes_ec_scrna, 
          xtree = F, ytree = F, legend.position = "bottom",
          id = 'Clusters')+ ggtitle("Stress signatures")
ggsave("stress_genes_ec_scATAC_jjDotPlot3.pdf", p, width =6, height = 5)


# Figure-S5
# A
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
ligands_of_Egfr <- c("Hsp90aa1","Tgfb1", "Grn","Hbegf", "Tgfa", "Itgb1", "Icam1", "Epha2")
p <- jjDotPlot(object = Endo_scRNA1, col.min = 0,base_size = 14, textSize = 14, 
               gene = unique(ligands_of_Egfr), 
               xtree = T, ytree = F, legend.position = "bottom", bar.width = 3, segWidth = 0.3,
               id = 'Clusters') + ggtitle("Egfr ligands")
ggsave("Egfr_signaling_Cluster_Endo_scRNA_jjDotPlot1.pdf", p, width =8, height =6)

# B
Endo_scATAC1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scATAC1.Rds")
ligands_of_Egfr <- c("Hsp90aa1","Tgfb1", "Grn","Hbegf", "Tgfa", "Itgb1", "Icam1", "Epha2")
p <- jjDotPlot(object = Endo_scATAC1, col.min = 0,base_size = 14, textSize = 14, 
               gene = unique(ligands_of_Egfr), 
               xtree = T, ytree = F, legend.position = "bottom", bar.width = 3, segWidth = 0.3,
               id = 'Clusters') + ggtitle("Egfr ligands")
ggsave("Egfr_signaling_Cluster_Endo_scATAC_jjDotPlot1.pdf", p, width =8, height =6)

# C
snRNA.apap2023.mouse.Endothelial <- readRDS("/data4/jinlianli/APAP/ANALYSIS/Nature2024_mouse_APAP_scRNA_Endo/snRNA.apap2023.mouse.Endothelial.Rds")
snRNA.apap2023.mouse.Endothelial$seurat_clusters <- as.numeric(snRNA.apap2023.mouse.Endothelial$seurat_clusters)
cols <- ggsci::pal_lancet()(7)
names(cols) <- levels(snRNA.apap2023.mouse.Endothelial$seurat_clusters)
dimp1 <- DimPlot(snRNA.apap2023.mouse.Endothelial, group.by = "seurat_clusters", cols = cols, pt.size = 0.8) + tidydr::theme_dr() + theme(element_text(size = 16), title = element_blank())
ggsave("snRNA.apap2023.mouse.Endothelial_clusters_Dimplot.pdf", dimp1, width = 5, height = 4.5)

MARKERS <- c("Apold1")
f.names <- "Endo_mouse_snRNA_apap_Apold1"
ncol = 1
fp <- FeaturePlot(snRNA.apap2023.mouse.Endothelial, features = MARKERS, ncol = ncol, order = T, pt.size = 0.8) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
ggsave(paste(f.names, "FeaturePlot.pdf", sep = "."), fp, width = ncol*4, height = ceiling(length(MARKERS)/ncol)*4)
MARKERS <- c("Hbegf")
f.names <- "Endo_mouse_snRNA_apap_Hbegf"
ncol = 1
fp <- FeaturePlot(snRNA.apap2023.mouse.Endothelial, features = MARKERS, ncol = ncol, order = T, pt.size = 0.8) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
ggsave(paste(f.names, "FeaturePlot.pdf", sep = "."), fp, width = ncol*4, height = ceiling(length(MARKERS)/ncol)*4)

# D
# 数据来源于Matchett et al., 2024, Nature
Endo_human_snRNA_apap <- readRDS("/data4/jinlianli/APAP/ANALYSIS/Endothelial_snRNA_human_2023/APAP_human.endothelial_raw.Rds")
Endo_human_snRNA_apap$seurat_clusters <- as.numeric(Endo_human_snRNA_apap$seurat_clusters)
Endo_human_snRNA_apap$seurat_clusters2 <- as.numeric(Endo_human_snRNA_apap$seurat_clusters)
Endo_human_snRNA_apap$seurat_clusters2[which(Endo_human_snRNA_apap$seurat_clusters2 == "7")] <- "5"
Endo_human_snRNA_apap$seurat_clusters2[which(Endo_human_snRNA_apap$seurat_clusters2 == "8")] <- "7"
Endo_human_snRNA_apap$seurat_clusters2[which(Endo_human_snRNA_apap$seurat_clusters2 == "9")] <- "8"
pdf("Human_APAP_2023_seurat_cluster2_UMAP.pdf", width = 5, height = 5)
DimPlot(Endo_human_snRNA_apap, cols = ggsci::pal_lancet()(9), group.by = "seurat_clusters2") & tidydr::theme_dr() + theme(title = element_blank())
dev.off()

MARKERS <- c("APOLD1")
f.names <- "Endo_Human_snRNA_apap_APOLD1"
ncol = 1
fp <- FeaturePlot(Endo_human_snRNA_apap, features = MARKERS, ncol = ncol, order = T, pt.size = 0.8) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
ggsave(paste(f.names, "FeaturePlot.pdf", sep = "."), fp, width = ncol*4, height = ceiling(length(MARKERS)/ncol)*4)
MARKERS <- c("HBEGF")
f.names <- "Endo_human_snRNA_apap_Hbegf"
ncol = 1
fp <- FeaturePlot(Endo_human_snRNA_apap, features = MARKERS, ncol = ncol, order = T, pt.size = 0.8) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
ggsave(paste(f.names, "FeaturePlot.pdf", sep = "."), fp, width = ncol*4, height = ceiling(length(MARKERS)/ncol)*4)

# E
Phx_Endo <- readRDS("/home/lijinlian/test/PHX_scRNA_miguel/Phx_Endo_SCT_miguel.Rds")
Apold1_Hbegf_exp <- Phx_Endo@assays$originalexp$data[c("Apold1", "Hbegf"),]
Apold1_Hbegf_exp <- as.data.frame(t(as.data.frame(Apold1_Hbegf_exp)))
Apold1_Hbegf_exp$coexp <- Apold1_Hbegf_exp$Apold1*Apold1_Hbegf_exp$Hbegf
Apold1_Hbegf_exp$Group <- as.data.frame(Phx_Endo$group)[,1]
p1 <- ggplot(Apold1_Hbegf_exp, aes(x = Apold1, y = Hbegf, color = coexp)) +
  geom_point(alpha = 0.6, size = 2) +  
  scale_color_viridis_c()+
  labs(title = "Correlation of Apold1 and Hbegf",
       x = "Apold1",
       y = "Hbegf") +
  theme_minimal() + 
  facet_grid(.~Group) +
  theme(axis.text = element_text(size = 12,colour = 'black'),
        axis.title = element_text(size = 12,colour = 'black',face = 'plain'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=12,hjust=0.5),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"))
ggsave("PHx_Endo_Apold1_Hbegf_exp_scatterplot_split.pdf", p1, width = 24, height = 4)


# Figure6
# A
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
eRegulon_metadata_result3 <- read_csv("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/out_4/eRegulon_metadata_result3.csv")
Endo_TFs <- unique(eRegulon_metadata_result3$TF)
Endo_scRNA1$Two_group <- as.character(Endo_scRNA1$Two_group)
Endo_scRNA1$Two_group[which(Endo_scRNA1$Two_group == "Cdkn1a+ Lsec")] <- "Apold1+ Lsec(Stress)"
Endo_scRNA1$Two_group[which(Endo_scRNA1$Two_group == "Other")] <- "Other"
Endo_scRNA1$Two_group <- factor(Endo_scRNA1$Two_group, levels = c("Apold1+ Lsec(Stress)", "Other"))
library(scRNAtoolVis)
id_obj = 'Two_group'
p_tfscrna <- jjDotPlot(object = Endo_scRNA1,
                       gene = rev(Endo_TFs), 
                       col.min = 0, 
                       base_size = 14, textSize = 14, 
                       bar.width = 2,segWidth = 0.3, 
                       x.text.angle = 45,  
                       x.text.hjust = 1,x.text.vjust = 1,
                       xtree = F, ytree = F, 
                       legend.position = "right",
                       id = id_obj) + coord_flip()
a=p_tfscrna$data
a$id <- as.character(a$id)
gene_df <- data.frame("gene" = unique(a$gene))
seq_c <- seq(1, nrow(a),2) #2是两次重复，两个组
for (i in seq_along(seq_c)) {
  value_avg_exp <- a$avg.exp[seq_c[i]]-a$avg.exp[c(seq_c[i]+1)]
  gene_df$cluster[i] <- ifelse(value_avg_exp > 0, a$id[seq_c[i]], a$id[c(seq_c[i]+1)])
}
# cluster select
id_v <- unique(gene_df$cluster)
id_1_df <- subset(a, id %in% c(id_v[1]))
id_2_df <- subset(a, id %in% c(id_v[2]))
# gene select
id_1_df <- subset(id_1_df, gene %in% subset(gene_df, cluster == id_v[1])$gene)
id_2_df <- subset(id_2_df, gene %in% subset(gene_df, cluster == id_v[2])$gene)
# order
id_1_df <- id_1_df[order(id_1_df$pct.exp,decreasing = T),]
id_2_df <- id_2_df[order(id_2_df$pct.exp,decreasing = T),]
Lsec_Hbegf_specific_tfs_df <- id_1_df
jfile_pre = "Endo_scRNA_Tfs2"
id_obj = 'Two_group'
p_tfscrna <- jjDotPlot(object = Endo_scRNA1,
                       gene = c(id_1_df$gene,id_2_df$gene), 
                       col.min = 0, 
                       base_size = 14, textSize = 14, 
                       segWidth = 0.3, bar.width = 2.5,
                       plot.margin = c(0,0.5,0,0),
                       x.text.angle = 45,  
                       x.text.hjust = 1,x.text.vjust = 1,
                       xtree = F, ytree = F, 
                       legend.position = "top",
                       id = id_obj)
ggsave(paste0(jfile_pre, "_jjDotPlot.pdf"), p_tfscrna, width = 8.7, height =3)

# B
# 在5t硬盘中原始pdf的文件夹下内皮细胞中的Endo_netework_25_4.cys，可自己修改图

# C
Endo_scATAC1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scATAC1.Rds")
consensus_regions <- read.delim("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/consensus_peak_calling/consensus_regions.bed", header=FALSE)
library(regioneR)
gr1 <- toGRanges(consensus_regions)
DefaultAssay(Endo_scATAC1) <- "ATAC"
pdf("Hbegf_cluster_CoveragePlot.pdf", width = 9, height = 4)
CoveragePlot(
  object = Endo_scATAC1,
  region = "Hbegf", group.by = "annotation6",
  ranges = gr1,
  features = "Hbegf",
  extend.downstream = 5000, extend.upstream = 2000,
  ranges.title = "MACS2"
)
dev.off()
pdf(paste(Coverage_gene,"Apold1_Lsec_Group_CoveragePlot_macs2.pdf", sep = "_"), width = 9, height = 4.5)
CoveragePlot(
  object = subset(Endo_scATAC1, annotation6 == "Apold1_Lsec"),
  region.highlight = regions_highlight,
  region = "Hbegf",
  ranges = gr1, peaks = F,ranges.title = "Macs2",
  group.by = "Group",
  features = "Hbegf",
  extend.upstream = 2000,
  extend.downstream = 5000, 
  ncol = 1
) + scale_fill_manual(values = c("red4", "green4"))
dev.off()

# D
Endo_scATAC1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scATAC1.Rds")
library(JASPAR2020) #JASPAR2020
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10) #BSgenome.Hsapiens.UCSC.hg19
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
Endo_scATAC1 <- AddMotifs(
  object = Endo_scATAC1,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
motif_gene_n <- "FOSL1"
pdf(paste(motif_gene_n,"motif.pdf", sep = "_"),width = 3.5, height = 3)
MotifPlot(
  object = Endo_scATAC1,
  motifs = "MA0477.2"
)
dev.off()

# E
Endo_scATAC1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scATAC1.Rds")
Endo_scATAC1 <- Footprint(
  object = Endo_scATAC1, in.peaks = T, 
  motif.name = c("FOSL1"), 
  genome = BSgenome.Mmusculus.UCSC.mm10
)
Endo_scATAC1$Two_group <- factor(Endo_scATAC1$Two_group, levels = c("Apold1+ Lsec(Stress)", "Other"))
Idents(Endo_scATAC1) <- Endo_scATAC1$Two_group
pdf("Endo_scATAC1_FOSL1_PlotFootprint.pdf", width = 5.6, height = 4)
PlotFootprint(Endo_scATAC1, features = c("FOSL1"), label = F)
dev.off()
Apold1_Lsec_obj <- Footprint(
  object = subset(Endo_scATAC1, annotation6 == "Apold1_Lsec"), in.peaks = T, 
  motif.name = c("FOSL1"), 
  genome = BSgenome.Mmusculus.UCSC.mm10
)
pdf("Endo_scATAC1_Apold1_Lsec_FOSL1_Group_split_PlotFootprint.pdf", width = 5, height = 4)
PlotFootprint(Apold1_Lsec_obj, features = c("FOSL1"), label = F, group.by =  "Group")
dev.off()

# F
Endo_scRNA1 <- readRDS("/data2/lijinlian/APAP_project/Endo_scRNA1.Rds")
gene <- "Fosl1"
fp <- FeaturePlot(Endo_scRNA1, features = gene, ncol = 1, order = T, pt.size = 1.2) & viridis::scale_colour_viridis(option = "plasma")& tidydr::theme_dr()
ggsave(paste0(gene,"_FeaturePlot.pdf"), fp, width = 4, height = 4)


# author:lijinlian