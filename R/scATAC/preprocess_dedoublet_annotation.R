library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(future)
plan("multisession", workers = 4)

counts_0h <- Read10X_h5("/data5/lijinlian22/APAP/cellranger_ATAC/APAP0h_atac.2.1.0/APAP_0h/outs/filtered_peak_bc_matrix.h5")
counts_12h <- Read10X_h5("/data5/lijinlian22/APAP/cellranger_ATAC/APAP12h_atac.2.1.0/outs/filtered_peak_bc_matrix.h5")
counts_96h <- Read10X_h5("/data5/lijinlian22/APAP/cellranger_ATAC/APAP96h_atac.2.1.0/APAP_96h/outs/filtered_peak_bc_matrix.h5")

# metadata
metadata_0h <- read.csv(
  file = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP0h_atac.2.1.0/APAP_0h/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_12h <- read.csv(
  file = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP12h_atac.2.1.0/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
metadata_96h <- read.csv(
  file = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP96h_atac.2.1.0/APAP_96h/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

## Create ChromatinAssay
# ChromatinAssay
assay_0h <- CreateChromatinAssay(
  counts = counts_0h,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP0h_atac.2.1.0/APAP_0h/outs/fragments.tsv.gz",
  min.cells = 10,                                                                                           ##这个阈值参考官网，可以修改
  min.features = 200
)
assay_12h <- CreateChromatinAssay(
  counts = counts_12h,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP12h_atac.2.1.0/outs/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)
assay_96h <- CreateChromatinAssay(
  counts = counts_96h,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP96h_atac.2.1.0/APAP_96h/outs/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

## create seurat object
seurat_0h <- CreateSeuratObject(
  counts = assay_0h,
  assay = 'peaks',
  project = '0h',                  ##指定后面画图的坐标SeuratProject，如下面的QC图
  meta.data = metadata_0h
)
seurat_12h <- CreateSeuratObject(
  counts = assay_12h,
  assay = 'peaks',
  project = '12h',
  meta.data = metadata_12h
)
seurat_96h <- CreateSeuratObject(
  counts = assay_96h,
  assay = 'peaks',
  project = '96h',
  meta.data = metadata_96h
)

## gene annotations
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
head(annotations)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
head(annotations)
# add the gene information to the object
Annotation(seurat_0h) <- annotations
Annotation(seurat_12h) <- annotations
Annotation(seurat_96h) <- annotations
# show annotation
table(seurat_0h@assays[["peaks"]]@annotation$gene_biotype)
table(seurat_0h@assays[["peaks"]]@annotation$type)

## Quality Control
# 0h  add blacklist ratio and fraction of reads in peaks
seurat_0h <- NucleosomeSignal(object = seurat_0h)
seurat_0h <- TSSEnrichment(seurat_0h, fast = FALSE)
seurat_0h$pct_reads_in_peaks <- seurat_0h$peak_region_fragments / seurat_0h$passed_filters * 100
seurat_0h$blacklist_ratio <- seurat_0h$blacklist_region_fragments / seurat_0h$peak_region_fragments
# 12h
seurat_12h <- NucleosomeSignal(object = seurat_12h)
seurat_12h <- TSSEnrichment(seurat_12h, fast = FALSE)
seurat_12h$pct_reads_in_peaks <- seurat_12h$peak_region_fragments / seurat_12h$passed_filters * 100
seurat_12h$blacklist_ratio <- seurat_12h$blacklist_region_fragments / seurat_12h$peak_region_fragments
# 96h
seurat_96h <- NucleosomeSignal(object = seurat_96h)
seurat_96h <- TSSEnrichment(seurat_96h, fast = FALSE)
seurat_96h$pct_reads_in_peaks <- seurat_96h$peak_region_fragments / seurat_96h$passed_filters * 100
seurat_96h$blacklist_ratio <- seurat_96h$blacklist_region_fragments / seurat_96h$peak_region_fragments

## subset high quality cells
sce.list <- list("0h" = seurat_0h, "12h" = seurat_12h, "96h" = seurat_96h)
subset.sce.list <- list()
for(i in 1:length(sce.list)){
  obj <- sce.list[[i]]
  name <- names(sce.list[i])
  obj <- subset(obj,subset = peak_region_fragments > 1000 &
                  peak_region_fragments < 20000 &
                  pct_reads_in_peaks > 40 &
                  blacklist_ratio < 0.05 &
                  nucleosome_signal < 2 &
                  TSS.enrichment > 4)
  subset.sce.list[[name]] <- obj
}

## Normalization and linear dimensional reduction
subset.sce.normalize.list <- list()
for (i in 1:length(subset.sce.list)) {
  obj <- subset.sce.list[[i]]
  name <- names(subset.sce.list[i])
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = 'q0')
  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:30)
  obj <- RunTSNE(obj,reduction = 'lsi',dims = 2:30)
  obj <- FindNeighbors(obj, reduction = 'lsi', dims = 2:30)
  obj <- FindClusters(obj, verbose = FALSE, algorithm = 3)
  subset.sce.normalize.list[[name]] <- obj
}

### merge data
## CreateFragmentObject
frags.APAP0h <- CreateFragmentObject(
  path = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP0h_atac.2.1.0/APAP_0h/outs/fragments.tsv.gz",
  cells = colnames(subset.sce.normalize.list[[1]])
)
frags.APAP12h <- CreateFragmentObject(
  path = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP12h_atac.2.1.0/outs/fragments.tsv.gz",
  cells = colnames(subset.sce.normalize.list[[2]])
)
frags.APAP96h <- CreateFragmentObject(
  path = "/data5/lijinlian22/APAP/cellranger_ATAC/APAP96h_atac.2.1.0/APAP_96h/outs/fragments.tsv.gz",
  cells = colnames(subset.sce.normalize.list[[3]])
)
## Quantify peaks granges in each dataset
# get peaks
gr.APAP0h <- granges(subset.sce.normalize.list[[1]])          ##等同于subset.sce.normalize.list[[1]]@assays$peaks
gr.APAP12h <- granges(subset.sce.normalize.list[[2]])
gr.APAP96h <- granges(subset.sce.normalize.list[[3]])
# combine peaks
combined.peaks <- reduce(x = c(gr.APAP0h, gr.APAP12h, gr.APAP96h))
#Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20] 
# Quantify:this function is parallelized using the future package
APAP0h.counts <- FeatureMatrix(
  fragments = frags.APAP0h,
  features = combined.peaks,
  cells = colnames(subset.sce.normalize.list[[1]])
)
APAP12h.counts <- FeatureMatrix(
  fragments = frags.APAP12h,
  features = combined.peaks,
  cells = colnames(subset.sce.normalize.list[[2]])
)
APAP96.counts <- FeatureMatrix(
  fragments = frags.APAP96h,
  features = combined.peaks,
  cells = colnames(subset.sce.normalize.list[[3]])
)
## create Seurat object
APAP0h_assay <- CreateChromatinAssay(APAP0h.counts, fragments = frags.APAP0h)
APAP0h <- CreateSeuratObject(APAP0h_assay, assay = "ATAC", meta.data = subset.sce.normalize.list[[1]]@meta.data) 
APAP12h_assay <- CreateChromatinAssay(APAP12h.counts, fragments = frags.APAP12h) 
APAP12h <- CreateSeuratObject(APAP12h_assay, assay = "ATAC", meta.data = subset.sce.normalize.list[[2]]@meta.data) 
APAP96h_assay <- CreateChromatinAssay(APAP96.counts, fragments = frags.APAP96h)
APAP96h <- CreateSeuratObject(APAP96h_assay, assay = "ATAC", meta.data = subset.sce.normalize.list[[3]]@meta.data) 
# add group/dataset message
APAP0h$Group <- '0h'
APAP12h$Group <- '12h'
APAP96h$Group <- '96h'

##merge object
combined <- merge(x = APAP0h, y = list(APAP12h, APAP96h), add.cell.ids = c("APAP0h", "APAP12h", "APAP96h")) # 213013*88415
combined$Group <- factor(combined$Group, levels = c('0h', '12h', '96h'))
combined[["ATAC"]]

## gene annotations
# extract gene annotations from EnsDb
library(EnsDb.Mmusculus.v79)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) #human:EnsDb.Hsapiens.v86
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(combined) <- annotations

## combined normalization
DefaultAssay(combined) <- 'ATAC'
combined <- RunTFIDF(combined)

## harmony
## test PC 
ATAC.PC.clustering <- function(Harmony.scATAC, mincutoff = "q25", var_harmony = "Group", PC = 10, resolution = 0.1, group = "Group"){
  library(Signac)
  library(Seurat)
  library(harmony)
  Harmony.scATAC <- FindTopFeatures(Harmony.scATAC, min.cutoff = mincutoff)    #上步用的是q0
  print(paste0("mincutoff:", mincutoff))
  print(paste0("variablefeatures:", length(VariableFeatures(Harmony.scATAC))))
  print(paste0("PC:", PC))
  name <- paste(mincutoff, length(VariableFeatures(Harmony.scATAC)), PC, sep = ".")
  Harmony.scATAC <- RunSVD(Harmony.scATAC, features = VariableFeatures(Harmony.scATAC))
  hm.integrated2 <- RunHarmony(
    object = Harmony.scATAC,
    group.by.vars = var_harmony,
    assay.use = 'ATAC',reduction = "lsi", project.dim = FALSE, verbose = FALSE
  )
  Harmony.scATAC <- RunUMAP(Harmony.scATAC, reduction = "harmony", dims = 2:PC, verbose = FALSE)
  Harmony.scATAC <- FindNeighbors(object = Harmony.scATAC, reduction = 'harmony', dims = 2:PC)
  Harmony.scATAC <- FindClusters(object = Harmony.scATAC, resolution = resolution, verbose = FALSE, algorithm = 3)
  p1 <- DimPlot(Harmony.scATAC, reduction = 'umap')+ggtitle(name)|DimPlot(Harmony.scATAC, group.by = group, reduction = 'umap')
  ggsave(paste(name, "Harmony.UMAP.pdf", sep = "."), p1, width = 14, height = 8)
  p2 <- DimPlot(Harmony.scATAC, reduction = 'tsne')+ggtitle(name)|DimPlot(Harmony.scATAC, group.by = group, reduction = 'tsne')
  ggsave(paste(name, "Harmony.tSNE.pdf", sep = "."), p2, width = 14, height = 8)
  return(Harmony.scATAC)
}
## run
PC.vector <- seq(15)
for (i in 1:length(PC.vector)) {
  hm.integrated <- ATAC.PC.clustering(combined, var_harmony = 'Group', mincutoff = "q25", PC = PC.vector[i], resolution = 0.1, group = "Group")
}

### compute gene activities
gene.activities <- GeneActivity(hm.integrated)
# add the gene activity matrix to the Seurat object as a new assay
hm.integrated[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
hm.integrated <- NormalizeData(
  object = hm.integrated,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(hm.integrated$nCount_ATAC)
)
DefaultAssay(hm.integrated) <- 'ACTIVITY'

# input for scRNA 
scRNA <- readRDS("APAP_scRNA_sce.Rds")  
transfer.anchors <- FindTransferAnchors(reference = scRNA, 
                                        query = hm.integrated, 
                                        features = VariableFeatures(object = scRNA),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY", 
                                        reduction = "cca") 
celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = scRNA@meta.data[,"celltype"], 
                                     weight.reduction = hm.integrated[["lsi"]], 
                                     dims = 2:30)   
hm.integrated <- AddMetaData(object = hm.integrated, metadata = celltype.predictions)
scATAC.data.high <- subset(hm.integrated, subset = prediction.score.max>=0.5)


#### dedoublet
library(ArchR)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(future)
addArchRThreads(threads = 6)
addArchRGenome("mm10")

# input
yyy =c("input_file/APAP_Oh.fragments.tsv.gz","input_file/APAP_12h.fragments.tsv.gz","input_file/APAP_96h.fragments.tsv.gz")
names(yyy) <- c("APAP_Oh","APAP_12h","APAP_96h")

ArrowFiles <- createArrowFiles(
  inputFiles = yyy,
  sampleNames = names(yyy),
  minTSS = 4, #Don't set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "APAP_output",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme2 <- filterDoublets(projHeme1, filterRatio = 1.5)

#Get all the cells in the specified sample 
#96h
idxSample.96h <- BiocGenerics::which(projHeme1$Sample %in% "APAP_96h")           #ArchR 
cellsSample.96h <- projHeme1$cellNames[idxSample.96h]

idxSampleFilter.96h <- BiocGenerics::which(projHeme2$Sample %in% "APAP_96h")
cellsSampleFiler.96h <- projHeme2$cellNames[idxSampleFilter.96h]
cell_diff.96h <- setdiff(cellsSample.96h,cellsSampleFiler.96h)
cell_diff.96h <- gsub("#","_",gsub("_","",cell_diff.96h))

#12h
idxSample.12h <- BiocGenerics::which(projHeme1$Sample %in% "APAP_12h")           #ArchR 
cellsSample.12h <- projHeme1$cellNames[idxSample.12h]

idxSampleFilter.12h <- BiocGenerics::which(projHeme2$Sample %in% "APAP_12h")
cellsSampleFiler.12h <- projHeme2$cellNames[idxSampleFilter.12h]
cell_diff.12h <- setdiff(cellsSample.12h,cellsSampleFiler.12h)
cell_diff.12h <- gsub("#","_",gsub("_","",cell_diff.12h))

#0h
idxSample.0h <- BiocGenerics::which(projHeme1$Sample %in% "APAP_Oh")           #ArchR 
cellsSample.0h <- projHeme1$cellNames[idxSample.0h]

idxSampleFilter.0h <- BiocGenerics::which(projHeme2$Sample %in% "APAP_Oh")
cellsSampleFiler.0h <- projHeme2$cellNames[idxSampleFilter.0h]
cell_diff.0h <- setdiff(cellsSample.0h,cellsSampleFiler.0h)
cell_diff.0h <- gsub("#","_",gsub("_","",cell_diff.0h))

##
# all doublet cells barcode
celldiff_combine <- c(cell_diff.0h,cell_diff.12h,cell_diff.96h)
scATAC.data.high$filter <- ifelse(rownames(scATAC.data.high@meta.data) %in% celldiff_combine, "Doublet","Normal_Cell")
scATAC.data.high.rmDoublets <-subset(x = scATAC.data.high, subset = filter == "Normal_Cell")

Idents(scATAC.data.high.rmDoublets) <- scATAC.data.high.rmDoublets$predicted.id
saveRDS(scATAC.data.high.rmDoublets, "scATAC.data.high.rmDoublets.Rds")


### Co-embedding
DefaultAssay(scATAC.data.high.rmDoublets) ="ACTIVITY"
transfer.anchors2 <- FindTransferAnchors(reference = scRNA, 
                                         query = scATAC.data.high.rmDoublets, 
                                         features = VariableFeatures(object = scRNA),
                                         reference.assay = "RNA",
                                         query.assay = "ACTIVITY", 
                                         reduction = "cca") 
genes.use <- transfer.anchors@anchor.features 
DefaultAssay(scRNA) <- "SCT"
scRNA.object <- NormalizeData(scRNA)
refdata <- GetAssayData(scRNA.object, assay = "SCT", slot = "data")[genes.use, ]
# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors2, refdata = refdata, weight.reduction = scATAC.data.high.rmDoublets[["lsi"]], dims = 2:30)
scATAC.data.high.rmDoublets[["RNA"]] <- imputation
DefaultAssay(scATAC.data.high.rmDoublets) <- "RNA"     

scRNA.object <- AddMetaData(scRNA.object, metadata = rep("scRNA", nrow(scRNA.object@meta.data)), col.name = "Datatype")
scATAC.object <- AddMetaData(scATAC.data.high.rmDoublets, metadata = rep("scATAC", nrow(scATAC.data.high.rmDoublets@meta.data)), col.name = "Datatype")

coembed <- merge(x = scRNA.object, y = scATAC.object)  
coembed$Group <- factor(coembed$Group, levels = c("0h", "6h", "12h", "24h", "48h", "96h"))
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use)
coembed <- RunHarmony(object = coembed, group.by.vars = "Datatype")
coembed <- RunUMAP(coembed, dims = 1:30, reduction = 'harmony')
DimPlot(coembed, group.by = "Datatype", reduction = "umap", label = F) |DimPlot(coembed, group.by = "celltype", reduction = "umap", label = T)
saveRDS(coembed, "coembed.Rds")

# End!
# coembed.Rds