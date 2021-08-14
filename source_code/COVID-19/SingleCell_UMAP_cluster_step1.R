library(Seurat)
library(cowplot)
library(patchwork)
library("hdf5r")
createOb <- function(sample){
  #rawcount <- read.csv(paste0(sample,'_hs_RSEC_MolsPerCell.csv'),comment.char="#",row.names=1)
  samplename <- gsub("_.*","",sample)
  dat <- Read10X_h5(sample)
  nCoV.seurat=CreateSeuratObject(counts = dat, project = samplename, min.cells = 3, min.features = 200)
###########
nCoV.seurat[['percent.mito']] <- PercentageFeatureSet(nCoV.seurat, pattern = "^MT-")
pdf(paste0(samplename,'feature_ncount_percentmito.pdf'),width=15)
print(VlnPlot(object = nCoV.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
nCoV.seurat <- subset(nCoV.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 20)
return(nCoV.seurat)

}
samples=list.files('./','*_filtered_feature_bc_matrix.h5')
sceList=lapply(samples,createOb)
sceList
options(future.globals.maxSize = 4000 * 1024^2)
for (i in 1:length(sceList)) {
sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE)
sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
}
names <- gsub('_.*','',samples)
names(sceList) <- names
sceList.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30)
sceList.integrated <- IntegrateData(anchorset = sceList.anchors, dims = 1:30)
library(data.table)
saveRDS(sceList.integrated,file="sceList.integrated_step1.RDS")

sceList.integrated <- ScaleData(sceList.integrated, verbose = FALSE)
sceList.integrated <- RunPCA(sceList.integrated, npcs = 50, verbose = FALSE)
COVIDRunner.integrated <- sceList.integrated
print(COVIDRunner.integrated[["pca"]], dims = 1:5, nfeatures = 5)
pdf('PCA_VizDimLoading.pdf')
VizDimLoadings(COVIDRunner.integrated, dims = 1:2, reduction = "pca")
dev.off()
pdf('PCA_Dim_plot.pdf')
DimPlot(COVIDRunner.integrated, reduction = "pca")
dev.off()

pdf('PC1_15.pdf')
DimHeatmap(COVIDRunner.integrated, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#COVIDRunner.integrated <- JackStraw(COVIDRunner.integrated, num.replicate = 100)
#COVIDRunner.integrated <- ScoreJackStraw(COVIDRunner.integrated, dims = 1:20)
#pdf('ElbowPlot.pdf')
#ElbowPlot(COVIDRunner.integrated)
#dev.off()

COVIDRunner.integrated <- FindNeighbors(COVIDRunner.integrated, dims = 1:50)
COVIDRunner.integrated <- FindClusters(COVIDRunner.integrated, resolution = 1)
head(Idents(COVIDRunner.integrated), 5)
COVIDRunner.integrated <- RunUMAP(COVIDRunner.integrated, dims = 1:50)
COVIDRunner.integrated <- RunTSNE(object = COVIDRunner.integrated, dims = 1:50, do.fast = TRUE)
pdf('tSNE_plot.pdf')
DimPlot(COVIDRunner.integrated, reduction = "tsne")
dev.off()
pdf('Umap_plot.pdf')
DimPlot(COVIDRunner.integrated, reduction = "umap")
dev.off()

saveRDS(COVIDRunner.integrated,file="sceList.integrated_step2.RDS")

sce.markers <- FindAllMarkers(COVIDRunner.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(sce.markers,file="sce.markers.RDS")
write.table(as.data.frame(sce.markers),file="clusters_markers.xls",sep="\t",quote = F)

