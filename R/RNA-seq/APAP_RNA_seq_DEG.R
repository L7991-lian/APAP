# 公共APAP的RNA-seq数据进行差异基因表达分析和功能富集分析。以及细胞周期基因的ssGSEA分析。

GSE195753_genes_fpkm_expression <- read.delim("/data2/lijinlian/APAP_project/spatial/APAP_RNA_seq_public/GSE195753_genes_fpkm_expression.txt")
expMatrix <- GSE195753_genes_fpkm_expression[, c(6, 16:19)]
expMatrix <- expMatrix[!duplicated(expMatrix$gene_name),]
rownames(expMatrix) <- expMatrix$gene_name
expMatrix <- expMatrix[,c(2:5)]
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)
tpms <- tpms[, c(2,4,1,3)]
colnames(tpms) <- c("control1", "control2", "APAP1", "APAP2")
head(tpms)

#### limma 差异表达基因分析
group_list=c(rep('Control',2),rep('APAP',2))
group_list <- factor(group_list,levels = c("Control","APAP"), ordered = F)
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
exprSet <- log2(exprSet+1)
dat <- exprSet
design=model.matrix(~factor(group_list))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
deg=topTable(fit,coef=2,adjust='BH',number = Inf)

if(T){
  logFC_t=1.2
  deg$g=ifelse(deg$P.Value > 0.05,'stable',
               ifelse( deg$logFC > logFC_t,'UP',
                       ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
  )
  table(deg$g)
  head(deg)
  deg$symbol=rownames(deg)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Mm.eg.db)
  head(df)
  DEG=deg
  head(DEG)
  
  DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
  head(DEG)
  
  saveRDS(DEG,file = 'APAP_6h_anno_DEG.Rds')
  gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
  gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
}
table(DEG$g)
# DOWN stable     UP 
# 208  39620    389  

gene_up1= DEG[DEG$g == 'UP',] 
gene_down1=DEG[DEG$g == 'DOWN',] 
select_up_down_mat <- tpms[c(match(gene_up1$symbol, rownames(tpms)), match(gene_down1$symbol, rownames(tpms))),]
pheatmap::pheatmap(select_up_down_mat, cluster_cols = F, cluster_rows = F, scale = "row", show_rownames = F, border_color = "white", cellwidth = 60, main = "Differential expressed genes between APAP and control")

DEG_up <- gene_up
if(T){
  enrichKK <- enrichGO(gene = DEG_up,
                       OrgDb = org.Mm.eg.db, 
                       readable = TRUE,
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
}
up_terms <- c("cell chemotaxis","myeloid leukocyte migration", "neutrophil migration", "reactive oxygen species metabolic process", "regulation of DNA-templated transcription in response to stress", "regulation of epithelial cell proliferation", "blood coagulation", "regeneration", "cell cycle G1/S phase transition", "regulation of DNA-templated transcription initiation")
# down_terms <- c("fatty acid metabolic process")
up_id_vector <- c()
for (i in 1:length(up_terms)) {
  up_id <- grep(paste("^", up_terms[i], "$", sep = ""), aa$Description)
  up_id_vector <- c(up_id_vector, up_id)
}
A = enrichKK@result[up_id_vector,]
A$log_p.adjust <- -log10(A$p.adjust)
A <- arrange(A, log_p.adjust)
A$Description <- factor(A$Description, levels = unique(A$Description))
pgo <- ggplot(A,
              aes(x=Description,y=log_p.adjust, fill= log_p.adjust)) +
  geom_bar(stat="identity", width=0.8) +
  scale_fill_viridis_c(option = "C", direction = -1)+
  coord_flip() +
  xlab("") +
  ylab("-log10(p.adjust)") +
  theme_bw()+
  labs(title = "GO enrichment of up-expression genes in APAP 6h") +
  scale_x_discrete(labels=function(x) stringr::str_wrap(x, width=40))+
  theme(text = element_text(size = 14, face = "bold"), 
        axis.title = element_text(size = 14), axis.text.y = element_text(size = 14))
pgo

### ssGSEA
##基因集需要是list为对象。
library(homologene)
library(GSVA)
gene_list <- list("S_genes" = human2mouse(cc.genes$s.genes)$mouseGene, "G2M_genes" = human2mouse(cc.genes$g2m.genes)$mouseGene)
ssgsea <- gsva(tpms, gene_list, method='ssgsea',
               kcdf='Gaussian',
               abs.ranking=TRUE)
# ssgsea
# control1 control2 APAP1 APAP2
# S_genes      7.659    7.479 7.424 7.318
# G2M_genes    6.659    6.707 6.794 6.841

library(tibble)
library(magrittr)
library(reshape2)
library(ggplot2)
library(ggpubr)

tmp1 <- ssgsea %>% t()%>% as.data.frame() %>%
  rownames_to_column("Sample")
tmp1$Group <- c(rep("Control",2),rep("APAP",2))
tmp1 <- melt(tmp1)
colnames(tmp1) <- c("Sample","Group","Signature","Score")
tmp1$Group <- factor(tmp1$Group, levels = c("Control", "APAP"))
pdf("ssGSEA_Sgenes_G2Mgenes_boxplot.pdf", width = 4, height = 5)
ggplot(tmp1,aes(Signature,Score)) + 
  geom_boxplot(aes(fill = Group),outlier.shape = 21)+
  theme_bw() +
  labs(x = NULL, y = "Score") +
  scale_fill_manual(values = c("blue", "red"))+
  stat_compare_means(aes(group = Group,label = after_stat(p.signif)),
                     method = "wilcox.test",
                     label = "p.format",
                     hide.ns = F)+
  theme(plot.margin=unit(c(1,1,1,1),'cm'),
        plot.title = element_text(size = 12,color="black",hjust = 0.5),
        axis.title = element_text(size = 12,color ="black"), 
        axis.text = element_text(size= 12,color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid=element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12))
dev.off()
