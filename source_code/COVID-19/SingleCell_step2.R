library(Seurat)
library(ggsci)
library(ggplot2)
library(RColorBrewer)
COVIDRunner.integrated <- readRDS('sceList.integrated_step2.RDS')

#######Supplementary Figure 4
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colourCount = length(unique( COVIDRunner.integrated@meta.data$seurat_clusters))
#######Supplementary Figure 4A
pdf('Total_Umap_last_version.pdf')
DimPlot(COVIDRunner.integrated) + scale_color_manual(values = getPalette(colourCount))
dev.off()

#######Supplementary Figure 4B
DefaultAssay(COVIDRunner.integrated) <- "RNA"

markersDotPlot <- c('CD79A','CD79B','MZB1','TNFRSF17','IGHG4','IGHA1','FCGR3B','CSF3R',
	'LILRA4','BCL11A','SERPINF1','TPPP3','KRT18','KRT19','ELF3','MUC4','IFITM3','S100A9',
	'LYZ','CST3','S100A8','FCN1','CD163','MARCO','FCGR1A',"CD1C",'LAMP3','CD3D',
	'CD3G','CD3E','CD2','CD96','CD6','CD8A','CD4','PDCD1','CTLA4','LAG3','TIGIT','HAVCR2',
	'GNLY','KLRD1','NKG7','PRF1','CCL5','ACE2','TMPRSS2')
pdf('All_clusters_DotPlot.pdf',width=12)
DotPlot(COVIDRunner.integrated, features = markersDotPlot)+theme(axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 90))+ 
scale_color_gradient(low = "white",high = "red")
dev.off()

###########Figure 4

celltypesAnno <- read.table('cells_type_v2.xls',header=T,sep="\t")
sampleinfo <- read.table('samples_group_info1.xls',sep="\t",header=T)
loc <- match(COVIDRunner.integrated@meta.data$orig.ident,sampleinfo$GEO_ACCESSION)
COVIDRunner.integrated@meta.data[,'Group'] = sampleinfo$patientGroup[loc]

loc <- match(COVIDRunner.integrated@meta.data$seurat_clusters,celltypesAnno$ClusterID)
COVIDRunner.integrated@meta.data[,'Celltypes'] = celltypesAnno$celltype_abbre[loc]
Idents(COVIDRunner.integrated) = 'Celltypes'

#######Figure 4A
pdf('Umap_groups.pdf',width=20)
DimPlot(COVIDRunner.integrated,split.by="Group") + scale_color_d3("category20")
dev.off()

##########Figure 4B
pdf('All_cell_types_markers3.pdf',width=11,height=9)
VlnPlot(COVIDRunner.integrated,features=c('MUC4','KRT18','IGHG4','CD79A','BCL11A','LILRA4',
	'CD1C','FCGR3B','CD163','MARCO','FCGR1A','FCGR3A','PDCD1','CTLA4','TIGIT','HAVCR2','LAG3','CD2','CD3E','CD3D','KLRD1'),stack=T,flip=T)+
scale_x_discrete(limits=c('NK','T','T_CD4','T_CD8','Exh T_CD8','Macro/T','Macro','Mono','Neu','DC','pDC','B','Plasma','Epi'))+
scale_fill_d3(palette='category20')
dev.off()

##########write data and metadat of severe COVID-19 for cellphnodb analysis
sce_group.list = SplitObject(COVIDRunner.integrated,split.by="Group")
Severe = sce_group.list$` Severe`
write.table(as.matrix(Severe@assays$RNA@data), 'BALF_Severe_cellphonedb_count.txt', sep='\t', quote=F)
meta_data <- cbind(rownames(Severe@meta.data), Severe@meta.data[,'Celltypes', drop=F])
write.table(meta_data, 'BALF_Severe_cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


##############sigmoid function
sigmoid = function(x,a=1){1/(1+exp(-a*x))}
##############Figure 5A
##############cytokine storm scores
cytokinestormMarkers = c('IL2', 'IL7', 'CSF3', 'CXCL10', 'CCL2', 'CCL3', 'TNF', 'IL6')
cytokinestormMarkersList = list(cytokinestormMarkers)
COVIDRunner.integrated <- AddModuleScore(COVIDRunner.integrated, features = cytokinestormMarkersList,ctrl = 5, name = 'Cytokine_Storm_score')
Cytokine_Storm_score_scale = scale(COVIDRunner.integrated@meta.data$Cytokine_Storm_score1)
COVIDRunner.integrated@meta.data[,'Cytokine_Storm_score_scale_sigmod'] = sigmoid(Cytokine_Storm_score_scale,4)
pdf('BALF_Cytokine_Storm_score.pdf',width=20)
FeaturePlot(COVIDRunner.integrated,features="Cytokine_Storm_score_scale_sigmod",split.by='Group',cols=c('white','red'))
dev.off()

#############immune exhaustion scores
exhaustionMarkers = c('CTLA4','LAG3','BTLA','HAVCR2','PDCD1','TIGIT')
exhaustionMarkersList = list(exhaustionMarkers)
COVIDRunner.integrated <- AddModuleScore(COVIDRunner.integrated, features = exhaustionMarkersList,ctrl = 5, name = 'Immune_Exhaustion_score')
Immune_Exhaustion_score_scale = scale(COVIDRunner.integrated@meta.data$Immune_Exhaustion_score1)
COVIDRunner.integrated@meta.data[,'Immune_Exhaustion_score_scale_sigmod'] = sigmoid(Immune_Exhaustion_score_scale,4)
pdf('BALF_Immune_Exhaustion_score.pdf',width=20)
FeaturePlot(COVIDRunner.integrated,features="Immune_Exhaustion_score_scale_sigmod",split.by='Group',cols=c('white','red'))
dev.off()

#################Figure4D
celltypes.list = SplitObject(COVIDRunner.integrated,split.by="Celltypes")
MacroT.sce = celltypes.list$MacroT
CD3D= colnames(MacroT@assays$RNA@counts)[MacroT@assays$RNA@counts['CD3D',] >0]
CD163= colnames(MacroT@assays$RNA@counts)[MacroT@assays$RNA@counts['CD163',] >0]
LAG3= colnames(MacroT@assays$RNA@counts)[MacroT@assays$RNA@counts['LAG3',] >0]
MacroTx=list(CD3D=CD3D,CD163=CD163,LAG3=LAG3)
saveRDS(MacroTx,file="MacroT_venn_dat.RDS")


MacroCD3D = colnames(Macro.sce@assays$RNA@counts)[Macro.sce@assays$RNA@counts['CD3D',] >0]
MacroLAG3 = colnames(Macro.sce@assays$RNA@counts)[Macro.sce@assays$RNA@counts['LAG3',] >0]
MacroCD163 = colnames(Macro.sce@assays$RNA@counts)[Macro.sce@assays$RNA@counts['CD163',] >0]
MacroX = list(CD3D=MacroCD3D,LAG3=MacroLAG3,CD163=MacroCD163)
saveRDS(MacroX,file="Macro_venn_dat.RDS")

############Macrophage re-clustering
############Figure 5C
library(Seurat)

Macro.sce = celltypes.list$Macro
DefaultAssay(Macro.sce) <- "integrated"
Macro.sce <- FindVariableFeatures(Macro.sce)
Macro.sce <- ScaleData(Macro.sce)
Macro.sce <- FindNeighbors(Macro.sce, dims = 1:50)
Macro.sce <- FindClusters(Macro.sce, resolution = 0.2)

Macro.sce <- RunTSNE(Macro.sce, dims = 1:50)

pdf('Macro_tsne.pdf',width=20)
DimPlot(Macro.sce,label=T,split.by="Group",reduction="tsne")+scale_color_d3()
dev.off()

DefaultAssay(Macro.sce) = "RNA"

cytokines = c('CXCL10','CXCL11','CCL2','CCL3','CCL4','TNF','TGFB1','IDO1','FCN1')


Macro.sce@meta.data[,'Macro_subclasses'] = paste0('C',Macro.sce@meta.data$seurat_clusters)
pdf('Macro_dotplot.pdf',height=4)
DotPlot(Macro.sce,group.by="Macro_subclasses", features = cytokines)+theme(axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 90))+ 
scale_color_gradient(low = "white",high = "red")
dev.off()


