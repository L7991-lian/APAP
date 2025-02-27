# 细胞类型随时间变化的差异表达基因分析、样本相关性分析和GO功能富集分析
# 肝脏驻留细胞类型，6小时的样本显著差异
# All为所有细胞类型的scRNA对象

### scRNA-seq
marker_list <- list()
library(dplyr)
for (m in seq_along(celltype_id)) {
  obj <- subset(All, celltype == celltype_id[m])
  Idents(obj) <- obj$Group
  min.pct = 0.1 
  logfc.threshold = 0.1
  diff.wilcox = FindAllMarkers(obj, only.pos = T, min.pct = min.pct, logfc.threshold = logfc.threshold)
  diff.wilcox_sig <- diff.wilcox %>% filter(p_val_adj <= 0.01)
  marker_list[[celltype_id[m]]] <- diff.wilcox_sig
}
marker_df_scRNA <- do.call(rbind, marker_list)
marker_scRNA_avg <- AverageExpression(All, features = unique(marker_df_scRNA$gene), group.by = c("celltype", "Group"), assays = "SCT")$SCT
marker_scRNA_avg <- as.data.frame(marker_scRNA_avg)

### scATAC-seq
snATAC_obj <- subset(coembed_data, Datatype == "snATAC")
snATAC_obj$Group <- droplevels(snATAC_obj$Group)
DefaultAssay(snATAC_obj) <- 'ACTIVITY'
marker_list_snatac <- list()
library(dplyr)
for (m in seq_along(celltype_id)) {
  obj <- subset(snATAC_obj, celltype == celltype_id[m])
  Idents(obj) <- obj$Group
  min.pct = 0.1 
  logfc.threshold = 0.1
  diff.wilcox = FindAllMarkers(obj, only.pos = T, min.pct = min.pct, logfc.threshold = logfc.threshold)
  diff.wilcox_sig <- diff.wilcox %>% filter(p_val_adj <= 0.01)
  marker_list_snatac[[celltype_id[m]]] <- diff.wilcox_sig
}
marker_df_snATAC <- do.call(rbind, marker_list_snatac)
marker_snATAC_avg <- AverageExpression(snATAC_obj, features = unique(marker_df_snATAC$gene), group.by = c("celltype", "Group"), assays = "ACTIVITY")$ACTIVITY
marker_snATAC_avg <- as.data.frame(marker_snATAC_avg)


### 相关性分析
a_scrna <- intersect(unique(c(rownames(marker_snATAC_avg), rownames(marker_scRNA_avg))), rownames(All))
a_scATAC <- intersect(unique(c(rownames(marker_snATAC_avg), rownames(marker_scRNA_avg))), rownames(snATAC_obj))
scrna_scatac_deg_inter <- intersect(a_scrna, a_scATAC)
# scRNA
bind_marker_scRNA_avg <- AverageExpression(All, features = scrna_scatac_deg_inter, group.by = c("celltype", "Group"), assays = "SCT")$SCT
bind_marker_scRNA_avg <- as.data.frame(bind_marker_scRNA_avg)
colnames(bind_marker_scRNA_avg) <- paste0("RNA_", colnames(bind_marker_scRNA_avg))
# scATAC
DefaultAssay(snATAC_obj) <- "RNA"
bind_marker_snATAC_avg <- AverageExpression(snATAC_obj, features = scrna_scatac_deg_inter, group.by = c("celltype", "Group"), assays = "RNA")$RNA
bind_marker_snATAC_avg <- as.data.frame(bind_marker_snATAC_avg)
colnames(bind_marker_snATAC_avg) <- paste0("ATAC_", colnames(bind_marker_snATAC_avg))
bind_marker_snATAC_avg <- bind_marker_snATAC_avg[rownames(bind_marker_scRNA_avg), ]
# CORRELATION
bind_marker_df <- cbind(bind_marker_scRNA_avg, bind_marker_snATAC_avg)
bind_marker_df_order <- bind_marker_df[, c(1:18, 49:54, 43:48, 19:30, 37:42, 31:36, #scRNA
                                           55:63, 79:81, 76:78, 64:69, 73:75, 70:72 #snATAC
)]
bind_marker_df_cor_order <- cor(bind_marker_df_order)
pheatmap::pheatmap(bind_marker_df_cor_order, scale = "none", border_color = "white", fontsize = 12,
                          cluster_rows = F, cluster_cols = F, show_rownames = F, 
                          treeheight_row = 0, main = "Correlation between samples of celltypes",
                          color = colorRampPalette(colors = c("blue","white","red3"))(100))


### 肝脏驻留细胞类型的6h差异基因热图和功能富集分析
# DEG pheatmap
hsc_ec_hep_maf_chol_6h_intersect <- c(subset(marker_list[[2]], cluster == "6h")$gene, subset(marker_list[[9]], cluster == "6h")$gene, subset(marker_list[[1]], cluster == "6h")$gene, subset(marker_list[[3]], cluster == "6h")$gene, subset(marker_list[[8]], cluster == "6h")$gene)
hsc_ec_hep_maf_chol_6h_scRNA_avg <- marker_scRNA_avg[hsc_ec_hep_maf_chol_6h_intersect, c(match(c("Hepatocyte_6h", "Endothelial_6h", "Macrophage_6h", "HSC_6h", "Cholangiocyte_6h"), colnames(marker_scRNA_avg)))]
colnames(hsc_ec_hep_maf_chol_6h_scRNA_avg) <- sub("_6h", "", colnames(hsc_ec_hep_maf_chol_6h_scRNA_avg))
pheatmap::pheatmap(hsc_ec_hep_maf_chol_6h_scRNA_avg, scale = "row", 
                          cluster_rows = T, cluster_cols = T, show_rownames = F, cutree_rows =4,
                          treeheight_row = 0,
                          fontsize = 12,
                          treeheight_col = 0,clustering_method = "average", cellwidth = 30,
                          border_color = "white")
# GO
hsc_ec_hep_maf_chol_6h_nocluter_df <- data.frame(gene = c(subset(marker_list[[2]], cluster == "6h")$gene, subset(marker_list[[9]], cluster == "6h")$gene, subset(marker_list[[1]], cluster == "6h")$gene, subset(marker_list[[3]], cluster == "6h")$gene, subset(marker_list[[8]], cluster == "6h")$gene), cluster = c(rep("Endothelial", length(c(subset(marker_list[[2]], cluster == "6h")$gene))), rep("HSC", length(c(subset(marker_list[[9]], cluster == "6h")$gene))), rep("Hepatocyte", length(c(subset(marker_list[[1]], cluster == "6h")$gene))), rep("Macrophage", length(c(subset(marker_list[[3]], cluster == "6h")$gene))), rep("Cholangiocyte", length(c(subset(marker_list[[8]], cluster == "6h")$gene))))) 
hsc_ec_hep_maf_chol_6h_nocluter_df_unique <- hsc_ec_hep_maf_chol_6h_nocluter_df %>% distinct()
hsc_ec_hep_maf_chol_6h_nocluter_df_unique$cluster <- factor(hsc_ec_hep_maf_chol_6h_nocluter_df_unique$cluster, levels = c("Endothelial","HSC" ,"Hepatocyte" , "Macrophage","Cholangiocyte"))
cluster_markers <- hsc_ec_hep_maf_chol_6h_nocluter_df_unique
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
hsc_hep_endo_maf_6h <- "hsc_ec_hep_maf_chol_6h_deg_notcluster"
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
pheatmap::pheatmap(df_wide, border_color = "white", fontsize = 12, cluster_rows = F, cluster_cols = F, cellwidth = 30)


### 不是肝脏本身驻留细胞类型
# DEG pheatmap
T_B_NK_Neu_6h_intersect <- c(subset(marker_list[[4]], cluster == "6h")$gene, subset(marker_list[[5]], cluster == "6h")$gene, subset(marker_list[[6]], cluster == "6h")$gene, subset(marker_list[[7]], cluster == "6h")$gene)
T_B_NK_Neu_6h_scRNA_avg <- marker_scRNA_avg[T_B_NK_Neu_6h_intersect, c(match(c("T_6h", "B_6h", "NK_6h", "Neutrophil_6h"), colnames(marker_scRNA_avg)))]
colnames(T_B_NK_Neu_6h_scRNA_avg) <- sub("_6h", "", colnames(T_B_NK_Neu_6h_scRNA_avg))
pheatmap::pheatmap(T_B_NK_Neu_6h_scRNA_avg[,c("T", "B","NK","Neutrophil")], scale = "row", 
                          cluster_rows = T, cluster_cols = F, show_rownames = F, 
                          # cutree_rows =4,
                          treeheight_row = 0,
                          fontsize = 12,
                          treeheight_col = 0,clustering_method = "complete", cellwidth = 30,
                          # color = colorRampPalette(colors = c("blue","white","red3"))(100),
                          border_color = "white")
# GO
T_B_NK_Neu_6h_nocluter_df <- data.frame(gene = c(subset(marker_list[[4]], cluster == "6h")$gene, subset(marker_list[[5]], cluster == "6h")$gene, subset(marker_list[[6]], cluster == "6h")$gene, subset(marker_list[[7]], cluster == "6h")$gene), cluster = c(rep("T", length(c(subset(marker_list[[4]], cluster == "6h")$gene))), rep("B", length(c(subset(marker_list[[5]], cluster == "6h")$gene))), rep("Neutrophil", length(c(subset(marker_list[[6]], cluster == "6h")$gene))), rep("NK", length(c(subset(marker_list[[7]], cluster == "6h")$gene))))) 
T_B_NK_Neu_6h_nocluter_df_unique <- T_B_NK_Neu_6h_nocluter_df %>% distinct()
T_B_NK_Neu_6h_nocluter_df_unique$cluster <- factor(T_B_NK_Neu_6h_nocluter_df_unique$cluster, levels = c("T","B" ,"NK" , "Neutrophil"))
write.table(T_B_NK_Neu_6h_nocluter_df_unique, "T_B_NK_Neu_6h_deg_nocluter_table.txt")
cluster_markers <- T_B_NK_Neu_6h_nocluter_df_unique
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
hsc_hep_endo_maf_6h <- "T_B_NK_Neu_6h_deg_notcluster"
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
result_df_2 <- data.frame(cluster = character(),
                          description = character(),
                          pvalue = numeric(),
                          qvalue = numeric(),
                          stringsAsFactors = FALSE)
for (i in 1:length(enrichment_results)) {
  cluster <- clusters[i]
  enr <- enrichment_results[[i]]
  
  cluster_result <- enr@result
  cluster_result$cluster <- cluster
  
  result_df_2 <- rbind(result_df_2, cluster_result)
}
# terms pheatmap
set_terms1 <- c(
  #T
  "GO:0140236", #translation at presynapse
  "GO:0042254", #ribosome biogenesis
  "GO:0000028", #ribosomal small subunit assembly
  "GO:0006914", #autophagy
  "GO:0022900", #electron transport chain
  "GO:0042773", #ATP synthesis coupled electron transport
  "GO:0001906", #cell killing
  "GO:0001913", #T cell mediated cytotoxicity
  "GO:0002475", #antigen processing and presentation via MHC class Ib
  "GO:0032611", #interleukin-1 beta production
  
  #B
  "GO:1902229", #regulation of intrinsic apoptotic signaling pathway in response to DNA damage
  "GO:0048002", #antigen processing and presentation of peptide antigen
  "GO:0019886", #antigen processing and presentation of exogenous peptide antigen via MHC class II
  "GO:0071353", #cellular response to interleukin-4
  "GO:0032609", #type II interferon production
  
  # NK
  "GO:0070231", #T cell apoptotic process
  "GO:0070228", #regulation of lymphocyte apoptotic process
  "GO:0034063", #stress granule assembly
  
  #Neutrophil
  "GO:0031098", #stress-activated protein kinase signaling cascade
  "GO:0071674", #mononuclear cell migration
  "GO:0032635", #interleukin-6 production
  "GO:0051403", #stress-activated MAPK cascade
  "GO:0060326", #cell chemotaxis
  "GO:0071706", #tumor necrosis factor superfamily cytokine production
  "GO:0002755" #MyD88-dependent toll-like receptor signaling pathway
)
term_ids <- c()
for (i in seq_along(set_terms1)) {
  id <- grep(set_terms1[i], result_df_2sig$ID)
  term_ids <- c(term_ids, id)
}
T_B_NK_Neu_6h_set_terms_df <- result_df_2sig[term_ids,]
T_B_NK_Neu_6h_set_terms_df <- T_B_NK_Neu_6h_set_terms_df[,c("Description","p.adjust","cluster")]
T_B_NK_Neu_6h_set_terms_df$log10_p.adjust <- -log10(T_B_NK_Neu_6h_set_terms_df$p.adjust)
T_B_NK_Neu_6h_set_terms_df <- T_B_NK_Neu_6h_set_terms_df[, c("cluster", "log10_p.adjust", "Description")]

df_wide <- T_B_NK_Neu_6h_set_terms_df %>%
  pivot_wider(names_from = cluster, values_from = c(log10_p.adjust))
rownames(df_wide) <- df_wide$Description
rn <- df_wide$Description
df_wide <- df_wide[, c(2,3,5,4)]
rownames(df_wide) <- rn
pheatmap::pheatmap(df_wide, border_color = "white", fontsize = 12, cluster_rows = F, cluster_cols = F, cellwidth = 30)



