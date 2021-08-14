library(NMF)
#argv <- commandArgs(T)
argvs <- commandArgs(T)
load('NMFbrunet.RData')
nk = as.numeric(argvs[1])
print(nk)
estim.r <- nmf(gene_exp_matrix_filter, nk, nrun=40, seed=123456)
w <- basis(estim.r)
h <- coef(estim.r)


library(GSVA)
gmt <- read.table('/home/myang/ALL/old/RNA/TCGA/estimate_genesignatures.txt',header=T)
gmt1 <- gmt[,-1]
gmt2 <- as.data.frame(t(gmt1))
gmt.list=as.list(gmt2)
es.dif <- gsva(gene_exp_matrix_filter, gmt.list, mx.diff=TRUE, verbose=FALSE, parallel.sz=1,method="ssgsea")
NMFconsensushc <- consensushc(estim.r,dendrogram=FALSE)
NMFconsensushc
pdf('coefmap.pdf')
coefmap(estim.r)
dev.off()
pdf('consensusmap.pdf')
consensusmap(estim.r)
dev.off()
NMFconsushctree=as.data.frame(cutree(NMFconsensushc,k=nk))
colnames(NMFconsushctree) <- 'NMFconsushctree'
write.table(data.frame(samples=rownames(NMFconsushctree),NMFconsushctree),file="NMFconsushctree_cluster.xls",quote=F,sep="\t")
nmfcluster <- read.table('NMFconsushctree_cluster.xls',header=T)
#nmfcluster <- read.table('NMFconsushctree_cluster.xls',header=T)
nmfcluster <- nmfcluster[order(nmfcluster[,2]),]

loc <- match(nmfcluster$samples,colnames(h))
h1 <- h[,loc]
write.table(as.data.frame(h1),file="samples_nmf_loading_score.csv",sep="\t",quote=F)
library(ComplexHeatmap)
ha=HeatmapAnnotation(NMF=as.character(nmfcluster[,2]),annotation_name_side = "left")
pdf('NMF_samples_loading_score.pdf')
Heatmap(h1,border=T,height = unit(10, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")
dev.off()
s <- extractFeatures(estim.r)
for(i in 1:length(s)){
genes=rownames(gene_exp_matrix_filter)[s[[i]]]
df = as.data.frame(genes)
write.table(df,file=paste0("NMF_group",as.character(i),'_metagenes.xls'),sep="\t",quote=F)
}


order <- match(nmfcluster$samples,colnames(es.dif))
es.diforder <- es.dif[,order]
write.table(es.diforder,file="immune_score_order.xls",sep="\t",quote=F)
library(ComplexHeatmap)
#Heatmap(es.diforder,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")+Heatmap(es.dif,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue1")
#pdf('immune_score.pdf')
#Heatmap(es.diforder,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")+Heatmap(es.dif,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue1")
#dev.off()
ha=HeatmapAnnotation(NMF=as.character(nmfcluster[,2]),annotation_name_side = "left")
pdf('immue_score_NMF.pdf')
Heatmap(es.diforder,border=T,height = unit(10, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")
dev.off()
nmfcluster_Clindat <- fread('nmfcluster_Clindat_filters.xls')
library(data.table)
nmfcluster_Clindat <- fread('nmfcluster_Clindat_filters.xls')
nmfcluster_Clindat <- fread('nmfcluster_Clindat.xls')
nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,'class'] = 'Immune'
nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,'class']
nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,]
dim(nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,])
library(survival)
my.surv <- Surv(as.numeric(nmfcluster_Clindat$survi_all_time),nmfcluster_Clindat$vital_status=='dead')
kmfit3 <- survfit(my.surv~nmfclusterClindat11$class,data=nmfclusterClindat11)
kmfit3 <- survfit(my.surv~nmfclusterClindat11$class,data=nmfcluster_Clindat)
kmfit3 <- survfit(my.surv~nmfcluster_Clindat$class,data=nmfcluster_Clindat)
cox.res <-coxph(Surv(as.numeric(nmfcluster_Clindat$survi_all_time)/30,nmfcluster_Clindat$status)~nmfcluster_Clindat$class)
cox.res <-coxph(Surv(as.numeric(nmfcluster_Clindat$survi_all_time)/30,nmfcluster_Clindat$vital_status=='dead')~nmfcluster_Clindat$class)
summary(cox.re)
summary(cox.res)
library(survminer)
ggsurvplot(kmfit3)
ggsurvplot(kmfit3,pval=T)
write.table(nmfcluster_Clindat,file="nmf_2_7_immune_cluster_clindat.xls",sep="\t",quote=F)
getwd()
write.table(nmfcluster_Clindat,file="nmf_2_7_immune_cluster_clindat.xls",sep="\t",quote=F,row.names=F)
