# Regeneration_hepatic_cell_scRNAseq为公共的部分肝切除损伤的小鼠scRNA-seq的seurat object

table(Regeneration_hepatic_cell_scRNAseq$annotation, Regeneration_hepatic_cell_scRNAseq$group)
#                 D0    8h    D1    D2    D3    D7
# PC hep         1921  1223  7571  7870  6644  1895
# PP hep        14398 20069 86135 53725 74220 31938
# Cyc. hep         94   119  1858 33474  5795   294
# HSC             227    47   112   124   310   111
# LSEC           5766 17038  5149  3461  2545  2376
# LVEC            457  1004   359   239   153   268
# Cyc. LSEC        25   111    14    79   465    15
# Cholangio        76   214    85    48    34    28
# Kupffer         533    86    31    56   257    34
# Capsule macro    73    34    22    35    43    53
# Mono            266   356   294   214   262   172
# Neutr           201   395   107   170   905    71
# NK              990    99   246   181    74   303
# NKT             334   168   163    61    43    56
# pDC             741   384   584   307   504   525
# Migra DC         73    21    14     5     7    33
# cDC             103   118   101    56    55   106
# ILC             608   528   399   199   119   120
# BC             3655   864  1364  1268  1834  4067
# Basophil          3     7    11    71    23    50
# Naive Cd4+ TC   703   129   186    97   111   295
# Naive Cd8+ TC  1038   230   317   174   197   480
# Cd4+ TC        5778  2850  2208   852   610  1296
# Cd8+ TC         440   189   410    77    56   192


### Hep
Phx_hep <- subset(Regeneration_hepatic_cell_scRNAseq, annotation == c("PC hep") | annotation ==  "PP hep" | annotation == "Cyc. hep")
options(future.globals.maxSize= 10000000*1024^2)
Phx_hep <- SCTransform(Phx_hep, assay = "originalexp", new.assay.name = "SCT", return.only.var.genes = FALSE, verbose = FALSE, variable.features.n = 5000)
Phx_hep <- RunPCA(Phx_hep, features = VariableFeatures(Phx_hep))
Phx_hep <- RunHarmony(Phx_hep, group.by.vars = c("group"))
Phx_hep <- RunUMAP(Phx_hep, dims = 1:25, reduction = "harmony")
# gene df
Igfbp2_Mfsd2a_exp <- Phx_hep@assays$originalexp$data[c("Igfbp2", "Mfsd2a"),]
Igfbp2_Mfsd2a_exp <- as.data.frame(t(as.data.frame(Igfbp2_Mfsd2a_exp)))
Igfbp2_Mfsd2a_exp$coexp <- Igfbp2_Mfsd2a_exp$Igfbp2*Igfbp2_Mfsd2a_exp$Mfsd2a
Igfbp2_Mfsd2a_exp$Group <- as.data.frame(Phx_hep$group)[,1]
#correlation and pvalue
mat_cor1 <- matrix(ncol = 2, nrow = 6)
for (i in seq_along(unique(Igfbp2_Mfsd2a_exp$Group))) {
  a <- subset(Igfbp2_Mfsd2a_exp, Group == unique(Igfbp2_Mfsd2a_exp$Group)[i])
  cor_value = cor(a$Igfbp2, a$Mfsd2a, method = "pearson", use = "complete.obs")
  p_value = cor.test(a$Igfbp2, a$Mfsd2a,method = "pearson")$p.value
  print(paste(i, c(cor_value, p_value), sep = "_"))
  mat_cor1[i,] <- c(cor_value, p_value)
}
rownames(mat_cor1) <- unique(Igfbp2_Mfsd2a_exp$Group)
colnames(mat_cor1) <- c("correlation", "p_value")
mat_cor1
# correlation      p_value
# D0  0.01087621 1.635226e-01
# 8h  0.06012403 1.307368e-18
# D1  0.03945627 3.037425e-34
# D2  0.03083092 1.937810e-21
# D3  0.01114069 1.039402e-03
# D7  0.01498295 5.641394e-03
# co-ratio
unname(table(Igfbp2_Mfsd2a_exp[Igfbp2_Mfsd2a_exp$coexp > 0,]$Group))/unname(table(Igfbp2_Mfsd2a_exp$Group))
# 0.028209346 0.128158423 0.010861831 0.005290894 0.008446901 0.004131626
# Igfbp2_Mfsd2a中表达再生相关受体的细胞比例
a <- subset(Igfbp2_Mfsd2a_exp,Igfbp2_Mfsd2a_exp$Igfbp2 > 0 & Igfbp2_Mfsd2a_exp$Mfsd2a > 0)
a1 <- a[a$regeneration_score > 0,]
nrow(a1)/nrow(a)
# [1] 0.7343889
#regeneration receptor score
object.name <- "Hep_Miguel_"
Phx_hep <- AddModuleScore(Phx_hep, features = list(unique(c("Egfr", "Cd81", "Ptprf", "Lrp1"))), name = "regeneration_receptor")
Igfbp2_Mfsd2a_exp$regeneration_score <- as.data.frame(Phx_hep$regeneration_receptor1)[,1]
#scatterplot
ggplot(Igfbp2_Mfsd2a_exp, aes(x = Igfbp2, y = Mfsd2a, color = regeneration_score)) +
  geom_point(alpha = 0.6, size = 2) +  
  scale_color_viridis_c()+
  labs(
    x = "Igfbp2 expression",
    y = "Mfsd2a expression") +
  theme_minimal() + 
  # facet_grid(.~Group) +
  facet_wrap(~Group, ncol = 3) + 
  theme(axis.text = element_text(size = 10,colour = 'black'),
        axis.title = element_text(size = 10,colour = 'black',face = 'plain'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=12,hjust=0.5),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"))


### Endo
Phx_Endo <- subset(Regeneration_hepatic_cell_scRNAseq, annotation == c("LSEC") | annotation ==  "Cyc. LSEC" | annotation == "LVEC")
Phx_Endo <- SCTransform(Phx_Endo, assay = "originalexp", new.assay.name = "SCT", return.only.var.genes = FALSE, verbose = FALSE, variable.features.n = 5000)
Phx_Endo <- RunPCA(Phx_Endo, features = VariableFeatures(Phx_Endo))
Phx_Endo <- RunHarmony(Phx_Endo, group.by.vars = "group")
Phx_Endo <- RunUMAP(Phx_Endo, dims = 1:25, reduction = "harmony")
#gene df
Apold1_Hbegf_exp <- Phx_Endo@assays$originalexp$data[c("Apold1", "Hbegf"),]
Apold1_Hbegf_exp <- as.data.frame(t(as.data.frame(Apold1_Hbegf_exp)))
Apold1_Hbegf_exp$coexp <- Apold1_Hbegf_exp$Apold1*Apold1_Hbegf_exp$Hbegf
Apold1_Hbegf_exp$Group <- as.data.frame(Phx_Endo$group)[,1]
Stab2_exp <- Phx_Endo@assays$originalexp$data[c("Stab2"),]
Stab2_exp <- as.data.frame(t(Stab2_exp))
Apold1_Hbegf_exp$Stab2 <- Stab2_exp[,1]
#correlation and pvalue
mat_cor2 <- matrix(ncol = 2, nrow = 6)
for (i in seq_along(unique(Apold1_Hbegf_exp$Group))) {
  a <- subset(Apold1_Hbegf_exp, Group == unique(Apold1_Hbegf_exp$Group)[i])
  cor_value = cor(a$Apold1, a$Hbegf, method = "pearson", use = "complete.obs")
  p_value = cor.test(a$Apold1, a$Hbegf,method = "pearson")$p.value
  print(paste(i, c(cor_value, p_value), sep = "_"))
  mat_cor2[i,] <- c(cor_value, p_value)
}
rownames(mat_cor2) <- unique(Apold1_Hbegf_exp$Group)
colnames(mat_cor2) <- c("correlation", "p_value")
write.table(mat_cor2, "Apold1_Hbegf_exp_correlation_table.txt")
mat_cor2
#   correlation      p_value
# D0  0.01978391 1.179008e-01
# 8h  0.08018794 2.748209e-27
# D1  0.00593834 6.590805e-01
# D2  0.01761374 2.790280e-01
# D3  0.01913924 2.818940e-01
# D7  0.01553199 4.233711e-01
unname(table(Apold1_Hbegf_exp[Apold1_Hbegf_exp$coexp > 0,]$Group))/unname(table(Apold1_Hbegf_exp$Group))
# [1] 0.008802817 0.051671900 0.016298443 0.017200318 0.028137844 0.018427980
#scatterplot
ggplot(Apold1_Hbegf_exp, aes(x = Apold1, y = Hbegf, color = Stab2)) +
  geom_point(alpha = 0.6, size = 2) +  
  scale_color_viridis_c()+
  labs(
    x = "Apold1 expression",
    y = "Hbegf expression") +
  theme_minimal() + 
  # facet_grid(.~Group) +
  facet_wrap(~Group, ncol = 3) + 
  theme(axis.text = element_text(size = 10,colour = 'black'),
        axis.title = element_text(size = 12,colour = 'black',face = 'plain'),
        axis.ticks = element_line(size=1),
        plot.title = element_text(size=12,hjust=0.5),
        panel.grid.major = element_line(color = "#cacfd2", linetype = "dashed"))


### celltype samples' correlation 
Regeneration_hepatic_cell_scRNAseq$Celltype <- Regeneration_hepatic_cell_scRNAseq$annotation
Regeneration_hepatic_cell_scRNAseq$Celltype <- as.character(Regeneration_hepatic_cell_scRNAseq$Celltype)
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "PC hep")] <- "Hepatocyte"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "PP hep")] <- "Hepatocyte"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Cyc. hep")] <- "Hepatocyte"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "LSEC")] <- "Endothelial"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "LVEC")] <- "Endothelial"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Cyc. LSEC")] <- "Endothelial"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Kupffer")] <- "Macrophage"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Capsule macro")] <- "Macrophage"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Mono")] <- "Macrophage"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "pDC")] <- "Dendritic"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Migra DC")] <- "Dendritic"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "cDC")] <- "Dendritic"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Naive Cd4+ TC")] <- "T"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Naive Cd8+ TC")] <- "T"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Cd4+ TC")] <- "T"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "Cd8+ TC")] <- "T"
Regeneration_hepatic_cell_scRNAseq$Celltype[which(Regeneration_hepatic_cell_scRNAseq$annotation == "BC")] <- "B"
table(Regeneration_hepatic_cell_scRNAseq$Celltype)
#     B    Basophil   Cholangio   Dendritic Endothelial  Hepatocyte         HSC         ILC  Macrophage 
# 13052         165         485        3737       39524      349243         931        1973        2821 
# Neutr          NK         NKT           T 
# 1849        1893         825       18915 
table(Regeneration_hepatic_cell_scRNAseq$Celltype, Regeneration_hepatic_cell_scRNAseq$group)
#                D0    8h    D1    D2    D3    D7
# B            3655   864  1364  1268  1834  4067
# Basophil        3     7    11    71    23    50
# Cholangio      76   214    85    48    34    28
# Dendritic     917   523   699   368   566   664
# Endothelial  6248 18153  5522  3779  3163  2659
# Hepatocyte  16413 21411 95564 95069 86659 34127
# HSC           227    47   112   124   310   111
# ILC           608   528   399   199   119   120
# Macrophage    872   476   347   305   562   259
# Neutr         201   395   107   170   905    71
# NK            990    99   246   181    74   303
# NKT           334   168   163    61    43    56
# T            7959  3398  3121  1200   974  2263

celltype_id <- unique(Regeneration_hepatic_cell_scRNAseq$Celltype)
celltype_id <- celltype_id[-which(celltype_id  == "Basophil")]
marker_list <- list()
library(dplyr)
for (m in seq_along(celltype_id)) {
  obj <- subset(Regeneration_hepatic_cell_scRNAseq, Celltype == celltype_id[m])
  Idents(obj) <- obj$group
  min.pct = 0.1 
  logfc.threshold = 0.1
  diff.wilcox = FindAllMarkers(obj, only.pos = T, min.pct = min.pct, logfc.threshold = logfc.threshold)
  diff.wilcox_sig <- diff.wilcox %>% filter(p_val_adj <= 0.01)
  marker_list[[celltype_id[m]]] <- diff.wilcox_sig
}
marker_df_scRNA <- do.call(rbind, marker_list)
table(marker_df_scRNA$cluster)
#    D0    8h    D1    D2    D3    D7 
# 14824  3923  1098  2471  3774  2913 

marker_scRNA_avg <- AverageExpression(Regeneration_hepatic_cell_scRNAseq, features = unique(marker_df_scRNA$gene), group.by = c("Celltype", "group"), assays = "originalexp")$originalexp
marker_scRNA_avg <- as.data.frame(marker_scRNA_avg)
a2 <- marker_scRNA_avg[,c(31:36, 25:30, 49:54, 37:42, 7:12, 13:18, 55:60, 73:78, 67:72, 1:6, 61:66, 19:24, 43:48)]

annotation_col = data.frame(Group = rep(c("D0","8h", "D1", "D2", "D3", "D7"), 13), CellType = sub("^(.*)_.*", "\\1", colnames(a2)))
rownames(annotation_col) = colnames(a2)
annotation_col$CellType <- factor(annotation_col$CellType, levels = unique(sub("^(.*)_.*", "\\1", colnames(a2))))
annotation_col$Group <- factor(annotation_col$Group, levels = c("D0","8h", "D1", "D2", "D3", "D7"))
group_col <- setNames(ggsci::pal_bmj()(6), c("D0","8h", "D1", "D2", "D3", "D7"))
Celltype_col <- setNames(c("tomato1", '#4363d8', "dodgerblue4", "darkmagenta","hotpink" , "brown4", "burlywood1",  "cadetblue1", "cadetblue4", "chartreuse3",  "chartreuse4", "slateblue1" ,"slateblue4"), unique(sub("^(.*)_.*", "\\1", colnames(a2))))
ann_colors = list(CellType = Celltype_col, Group = group_col)
pheatmap::pheatmap(cor(a2),
                        cluster_rows = F, 
                        cluster_cols = F, 
                        border_color = "white",
                        annotation_col = annotation_col,
                        annotation_row = annotation_col,
                        annotation_colors = ann_colors,
                        show_rownames = F, show_colnames = F, 
                        # cellwidth = 15, cellheight = 15,
                        treeheight_row = 0, main = "Correlation between samples of celltypes after PHx",
                        color = colorRampPalette(colors = c("blue","white","red3"))(100))

#### 肝脏驻留细胞类型的6h差异基因热图和功能富集分析
# DEG pheatmap
hsc_ec_hep_maf_Bas_8h_intersect <- c(subset(marker_list[[6]], cluster == "8h")$gene, subset(marker_list[[5]], cluster == "8h")$gene, subset(marker_list[[9]], cluster == "8h")$gene, subset(marker_list[[7]], cluster == "8h")$gene, subset(marker_list[[2]], cluster == "8h")$gene)
hsc_ec_hep_maf_Bas_8h_scRNA_avg <- marker_scRNA_avg[hsc_ec_hep_maf_Bas_8h_intersect, c(match(c("Hepatocyte_8h", "Endothelial_8h", "Macrophage_8h", "HSC_8h", "Basophil_8h"), colnames(marker_scRNA_avg)))]
colnames(hsc_ec_hep_maf_Bas_8h_scRNA_avg) <- sub("_8h", "", colnames(hsc_ec_hep_maf_Bas_8h_scRNA_avg))
pheatmap::pheatmap(hsc_ec_hep_maf_Bas_8h_scRNA_avg, scale = "row", 
                       cluster_rows = T, cluster_cols = T, show_rownames = F, 
                       # cutree_rows =4,
                       treeheight_row = 0,
                       fontsize = 14,
                       treeheight_col = 0,clustering_method = "average", cellwidth = 30,
                       border_color = "white")

# GO
hsc_ec_hep_maf_Bas_8h_nocluter_df <- data.frame(gene = c(subset(marker_list[[6]], cluster == "8h")$gene, subset(marker_list[[5]], cluster == "8h")$gene, subset(marker_list[[9]], cluster == "8h")$gene, subset(marker_list[[7]], cluster == "8h")$gene, subset(marker_list[[2]], cluster == "8h")$gene), cluster = c(rep("Hepatocyte", length(c(subset(marker_list[[6]], cluster == "8h")$gene))), rep("Endothelial", length(c(subset(marker_list[[5]], cluster == "8h")$gene))), rep("Macrophage", length(c(subset(marker_list[[9]], cluster == "8h")$gene))), rep("HSC", length(c(subset(marker_list[[7]], cluster == "8h")$gene))), rep("Basophil", length(c(subset(marker_list[[2]], cluster == "8h")$gene))))) 
hsc_ec_hep_maf_Bas_8h_nocluter_df_unique <- hsc_ec_hep_maf_Bas_8h_nocluter_df %>% distinct()
cluster_markers <- hsc_ec_hep_maf_Bas_8h_nocluter_df_unique
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
hsc_hep_endo_maf_6h <- "hsc_ec_hep_maf_Bas_8h_deg_notcluster"
clusters <- unique(cluster_markers$cluster)
enrichment_results <- list()
for (cluster in clusters) {
  cluster_genes_symbol <- cluster_markers[cluster_markers$cluster == cluster, "gene"]
  print(head(cluster_genes_symbol))
  cluster_genes_entrez <- mapIds(org.Mm.eg.db, as.character(cluster_genes_symbol), column = "ENTREZID", keytype = "SYMBOL")
  cluster_genes_entrez <- cluster_genes_entrez[!is.na(cluster_genes_entrez)]
  if (length(cluster_genes_entrez) == 0) {
    next
  }
  res1 <- enrichGO(gene = cluster_genes_entrez,
                   OrgDb = org.Mm.eg.db, 
                   readable = TRUE,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
  enrichment_results[[as.character(cluster)]] <- res1
}
result_df <- data.frame(cluster = character(),
                        description = character(),
                        pvalue = numeric(),
                        qvalue = numeric(),
                        stringsAsFactors = FALSE)
for (i in 1:length(enrichment_results)) {
  cluster <- clusters[i]
  enr <- enrichment_results[[i]]
  
  cluster_result <- enr@result
  cluster_result$cluster <- cluster
  
  result_df <- rbind(result_df, cluster_result)
}
result_df_sig <- subset(result_df, p.adjust < 0.05)
set_terms1 <- c(
  # Hepatocyte
  "GO:0010038", #"response to metal ion"
  "GO:0006953", #"acute-phase response"
  "GO:0072577", #"endothelial cell apoptotic process"
  "GO:0002526", #"acute inflammatory response"
  "GO:0050817", # coagulation
  "GO:0042730", #"fibrinolysis"
  "GO:0031639", #"plasminogen activation"
  
  # Basophil
  "GO:0010631", #"epithelial cell migration"
  "GO:0043620", #"regulation of DNA-templated transcription in response to stress"
  # "GO:0050678", #"regulation of epithelial cell proliferation"
  # "GO:0010586", #"miRNA metabolic process"
  "GO:0071559", #"response to transforming growth factor beta"
  "GO:0035296", #"regulation of tube diameter"
  "GO:0030509", #"BMP signaling pathway"
  "GO:0032640", #"tumor necrosis factor production"
  "GO:0019216", #regulation of lipid metabolic process
  "GO:0009062", #fatty acid catabolic process
  "GO:0006635", #fatty acid beta-oxidation
  "GO:0006006", #glucose metabolic process
  "GO:0030301", #cholesterol transport
  
  
  # Mac
  "GO:0140236", #translation at presynapse
  "GO:0046651", #lymphocyte proliferation
  "GO:0032943", #mononuclear cell proliferation
  "GO:0050863", #regulation of T cell activation
  "GO:0042113", #B cell activation
  "GO:0007159", #leukocyte cell-cell adhesion
  "GO:0035924", #cellular response to vascular endothelial growth factor stimulus
  "GO:0042058", #"regulation of epidermal growth factor receptor signaling pathway"
  "GO:0032728", #"positive regulation of interferon-beta production"
  "GO:0032623", #interleukin-2 production
  "GO:0032661", #"regulation of interleukin-18 production" 
  "GO:0036120", # "cellular response to platelet-derived growth factor stimulus"
  "GO:0032674", #"regulation of interleukin-5 production"
  "GO:0043123", #"positive regulation of canonical NF-kappaB signal transduction"
  "GO:0048012", #"hepatocyte growth factor receptor signaling pathway"
  
  # Endo
  "GO:0006914", #autophagy
  "GO:0001667", # ameboidal-type cell migration
  "GO:0001935", #"endothelial cell proliferation"
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
  "GO:0072574", #hepatocyte proliferation
  # "GO:0032675", #regulation of interleukin-6 production
  #HSC
  "GO:0001666", #response to hypoxia
  "GO:0007263", # nitric oxide mediated signal transduction
  "GO:0043114", #regulation of vascular permeability
  "GO:0006979" #"response to oxidative stress"
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
head(df_wide)
df_wide <- df_wide[, -1]
rownames(df_wide) <- rn
df_wide <- df_wide[, c("Hepatocyte", "Basophil", "Macrophage", "Endothelial", "HSC")]
pheatmap::pheatmap(df_wide, border_color = "white", fontsize = 12, cluster_rows = F, cluster_cols = F, cellwidth = 30, cellheight = 20)


### cell cycle genes Score, line plot, each cell types
Regeneration_hepatic_cell_scRNAseq <- CellCycleScoring(Regeneration_hepatic_cell_scRNAseq,
                                                       s.features = human2mouse(Seurat::cc.genes$s.genes)$mouseGene,
                                                       g2m.features = human2mouse(Seurat::cc.genes$g2m.genes)$mouseGene,
                                                       set.ident = T)
cellcycle_score_df <- Regeneration_hepatic_cell_scRNAseq@meta.data[, c("Celltype", "S.Score", "G2M.Score", "group")]  
library(dplyr)
library(tidyr)
# 使用dplyr的group_by和summarise函数来计算每个组的S.Score均值
cellcycle_score_df_grouped <- cellcycle_score_df %>%
  group_by(Celltype, group) %>%
  summarise(S.Score.Mean = mean(S.Score, na.rm = TRUE))
cellcycle_score_df_grouped_width_df <- cellcycle_score_df_grouped %>% pivot_wider(names_from = group, values_from = S.Score.Mean)
cellcycle_score_df_grouped_width_df <- as.data.frame(cellcycle_score_df_grouped_width_df)
rownames(cellcycle_score_df_grouped_width_df) <- cellcycle_score_df_grouped_width_df$Celltype
cellcycle_score_df_grouped_width_scale <- cellcycle_score_df_grouped_width_df[,2:7]  #第一列是celltype，后面6个组group
cellcycle_score_df_normalized <- apply(cellcycle_score_df_grouped_width_scale, 1, function(x) x - x[1]) #对0小时归一化
#S.Score line plot for each celltype
library(ggplot2)
cellcycle_score_df_normalized <- as.data.frame(cellcycle_score_df_normalized)
cellcycle_score_df_normalized$group <- rownames(cellcycle_score_df_normalized)
cellcycle_score_df_long <- cellcycle_score_df_normalized %>%
  pivot_longer(cols = colnames(cellcycle_score_df_normalized)[1:13], names_to = "Celltype", values_to = "S.Score")
# cellcycle_score_df_long$Celltype <- factor(cellcycle_score_df_long$Celltype, levels = levels(All$Celltype))
cellcycle_score_df_long$group <- factor(cellcycle_score_df_long$group, levels = c("D0","8h", "D1", "D2", "D3", "D7"))
ggplot(cellcycle_score_df_long, aes(x = group, y = S.Score, group = Celltype, color = Celltype)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Group", y = "S.Score") +
  theme_minimal() +
  scale_color_manual(values =  as.vector(Celltype_col))+
  theme(text = element_text(size = 14), axis.title.x = element_blank())+
  facet_wrap(~Celltype, scales = "free_y")

#G2M.Score
cellcycle_score_df_grouped <- cellcycle_score_df %>%
  group_by(Celltype, group) %>%
  summarise(G2M.Score.Mean = mean(G2M.Score, na.rm = TRUE))
cellcycle_score_df_grouped_width_df <- cellcycle_score_df_grouped %>% pivot_wider(names_from = group, values_from = G2M.Score.Mean)
cellcycle_score_df_grouped_width_df <- as.data.frame(cellcycle_score_df_grouped_width_df)
rownames(cellcycle_score_df_grouped_width_df) <- cellcycle_score_df_grouped_width_df$celltype
cellcycle_score_df_grouped_width_scale <- cellcycle_score_df_grouped_width_df[,2:7]
cellcycle_score_df_normalized <- apply(cellcycle_score_df_grouped_width_scale, 1, function(x) x - x[1])

library(ggplot2)
cellcycle_score_df_normalized <- as.data.frame(cellcycle_score_df_normalized)
cellcycle_score_df_normalized$group <- rownames(cellcycle_score_df_normalized)
cellcycle_score_df_long <- cellcycle_score_df_normalized %>%
  pivot_longer(cols = colnames(cellcycle_score_df_normalized)[1:13], names_to = "Celltype", values_to = "G2M.Score")
cellcycle_score_df_long$group <- factor(cellcycle_score_df_long$group, levels = c("D0","8h", "D1", "D2", "D3", "D7"))
#S.Score line plot for each celltype
ggplot(cellcycle_score_df_long, aes(x = group, y = G2M.Score, group = Celltype, color = Celltype)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Group", y = "G2M.Score") +
  theme_minimal() +
  scale_color_manual(values =  as.vector(Celltype_col))+
  theme(text = element_text(size = 14), axis.title.x = element_blank())+
  facet_wrap(~Celltype, scales = "free_y")

###cellchat
load("/data2/lijinlian/APAP_project/cellchat/allcelltype/CellChatDB.mouse_user.rda")
db.new <- CellChatDB.mouse_user
write.table(hsc_ec_hep_maf_Bas_8h_scRNA_avg, "hsc_ec_hep_maf_Bas_8h_scRNA_avg.txt")
setwd("./cellchat_2")
run_cellchat <- function(a1, metadata, 
                         group, 
                         assay="originalexp", 
                         # species ="mm", 
                         min.cells = 20,
                         # search = c('Cell-Cell Contact','ECM-Receptor',
                         #            'Secreted Signaling'),
                         # color.use = color.use,
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
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE, type = "truncatedMean", trim = 0)
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat, thresh = 0.05)
  groupSize <- as.numeric(table(cellchat@idents))
  # par(mfrow = c(1,2), xpd=TRUE)
  pdf(paste(filename_pre,"netVisual_circle_count.pdf", sep = "_"), width = 8, height = 8)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  pdf(paste(filename_pre, "netVisual_circle_weight.pdf", sep = "_"), width = 8, height = 8)
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  # saveRDS(cellchat, paste(filename_pre, "cellchat_object.Rds", sep = "_"))
  return(cellchat)
}
Idents(Regeneration_hepatic_cell_scRNAseq) <- Regeneration_hepatic_cell_scRNAseq$Celltype
groups_vector <- c("D0","8h", "D1", "D2", "D3", "D7")
list_cellchat_object <- list()
for (i in 1:length(groups_vector)) {
  obj <- subset(Regeneration_hepatic_cell_scRNAseq, group == groups_vector[i])
  filename <- groups_vector[i]
  ccobj <- run_cellchat(obj, 
                        metadata = obj@meta.data, 
                        group = "Celltype", 
                        # color.use = colors,
                        min.cells = 20, 
                        filename_pre = filename)
  ccobj <- aggregateNet(ccobj, thresh = 0.05)
  ccobj <- netAnalysis_computeCentrality(ccobj, slot.name = "netP")
  list_cellchat_object[[groups_vector[i]]] <- ccobj
}
# 合并
groups_vector <- c("D0","8h", "D1", "D2", "D3", "D7")
cellchat <- mergeCellChat(list_cellchat_object, add.names = names(list_cellchat_object))

### heatmap: weight to hepatocyte
weight_list <- lapply(c("D0","8h", "D1", "D2", "D3", "D7"), function(x) {
  cellchat@net[[x]][["weight"]]/max(cellchat@net[[x]][["weight"]])
})
names(weight_list) <- c("D0","8h", "D1", "D2", "D3", "D7")
w_to_hep_0h <- weight_list[[1]][,6]
w_to_hep_6h <- weight_list[[2]][,6]
w_to_hep_12h <- weight_list[[3]][,6]
w_to_hep_24h <- weight_list[[4]][,6]
w_to_hep_48h <- weight_list[[5]][,6]
w_to_hep_96h <- weight_list[[6]][,6]
w_to_hep_l <- rbind(w_to_hep_0h,w_to_hep_6h,w_to_hep_12h,w_to_hep_24h,w_to_hep_48h,w_to_hep_96h)
rownames(w_to_hep_l) <- c("D0","8h", "D1", "D2", "D3", "D7")
write.table(w_to_hep_l, "PHx_celltype_to_hepatocyte_interaction_weight_table.txt")
#pheatmap
pheatmap::pheatmap(w_to_hep_l, cluster_rows = F, cluster_cols = F, border_color = "white", fontsize = 14)

# EGFR signaling
regeneration_proliferation <- c("Cd81","Egfr","Lrp1","Ptprf")
regeneration_proliferation_df <- findEnrichedSignaling(list_cellchat_object[[2]], features = regeneration_proliferation, pattern ="incoming", thresh = 0.05)
regeneration_proliferation_df_endo_to_hep <- subset(regeneration_proliferation_df, target == "Hepatocyte" & source == "Endothelial")
Egfr_endo_to_hep_signalling <- regeneration_proliferation_df_endo_to_hep[regeneration_proliferation_df_endo_to_hep$receptor == "Egfr",]
Egfr_endo_to_hep_signalling <- Egfr_endo_to_hep_signalling[order(Egfr_endo_to_hep_signalling$prob, decreasing = T),]
library(tidyverse)
library(circlize)
circos.par(start.degree = 90)
chordDiagram(
  Egfr_endo_to_hep_signalling, 
  directional = 1,
  # grid.col = grid.col, 
  annotationTrack = "grid", link.arr.width = 0.1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  # self.link = "1",
  # link.auto = F,
  # group = group, 
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
circos.clear()

# 如果所研究的信号没有被测到，可以采用如下函数进行探查，trim设为0.1或者0.05或者0
computeAveExpr(list_cellchat_object[[2]], features = c("Hbegf","Egfr"),type =  "truncatedMean",trim = 0)
# B  Basophil Cholangio  Dendritic Endothelial Hepatocyte       HSC        ILC Macrophage
# Hbegf 0.01078465 0.2444278 1.5310260 0.01148199  0.15742274 0.01803702 0.2552233 0.02857425 0.06229267
# Egfr  0.02158630 0.0000000 0.1160137 0.08882196  0.01675303 0.62416636 0.2849083 0.01816524 0.03061208
# Neutr        NK         NKT          T
# Hbegf 0.02805878 0.0114758 0.000000000 0.03794876
# Egfr  0.03766230 0.0000000 0.003104015 0.02322422
