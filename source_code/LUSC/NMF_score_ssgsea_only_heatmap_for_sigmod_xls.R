###USAGE:Rscript NMF_score_ssgsea_only_heatmap_for_sigmod_xls.R NMFconsushctree_cluster.xls CD8_Tcell_exhausted.xls LGG_CGGA_expr.csv
library(NMF)
library(data.table)
#argv <- commandArgs(T)
argvs <- commandArgs(T)
nmfcluster <- argvs[1]
gmtf <- argvs[2]


dat<- fread(argvs[3])
dat <- as.data.frame(dat)
#gene_exp_matrix <- dat[,-1]
#rownames(gene_exp_matrix) <- dat[,1]
#gene_exp_matrix <- as.matrix(gene_exp_matrix)
#print(head(gene_exp_matrix))
geneidfactor<-factor(dat[,1])
gene_exp_matrix<-apply(dat[,-1],2,function(x) tapply(x,geneidfactor,mean))
fwrite(cbind(genes=rownames(gene_exp_matrix),gene_exp_matrix),file="Uniq_genes_exp_matrix.csv")

library(GSVA)
gmt <- fread(gmtf,header=T)
print(gmt)
gmt.list=as.list(gmt)
for(i in 1:length(gmt.list)){
	gmt.list[[i]]=gmt.list[[i]][!gmt.list[[i]]==""]
}
print(gmt.list)
#library(GSA)
#gmts=GSA.read.gmt('/home/myang/soft/software/GSEA/c2.cp.kegg.v6.2.symbols.gmt')
#gmt.list=gmts$genesets
# gmt1 <- gmt[,-1]
# gmt2 <- as.data.frame(t(gmt1))
# gmt.list=as.list(gmt2)
es.dif <- gsva(gene_exp_matrix, gmt.list, mx.diff=TRUE, verbose=FALSE, parallel.sz=1,method="ssgsea")
print(es.dif)
sigmoid = function(x,a=1){1/(1+exp(-a*x))}
es.dif=t(scale(t(es.dif)))
es.dif=sigmoid(es.dif)
es.dif=as.data.frame(es.dif)
print(es.dif)
# x=es.dif
# center <- sweep(x, 1, apply(x, 1, min),'-')
# R <- apply(x, 1, max) - apply(x,1,min)
# es.dif<- sweep(center, 1, R, "/")

nmfcluster <- read.table(nmfcluster,header=T)
#nmfcluster <- read.table('NMFconsushctree_cluster.xls',header=T)
nmfcluster <- nmfcluster[order(nmfcluster[,2]),]
order <- match(nmfcluster$samples,colnames(es.dif))
es.diforder <- es.dif[,order]
write.table(es.diforder,file=paste0(argvs[2],'immue_score_NMF.xls'),quote=F,sep="\t")
library(ComplexHeatmap)
library(circlize)
f1 = colorRamp2(seq(min(es.diforder), max(es.diforder), length = 3), c("#dce6f1", "white","#ff0000"))
#Heatmap(es.diforder,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")+Heatmap(es.dif,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue1")
#pdf('immune_score.pdf')
#Heatmap(es.diforder,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")+Heatmap(es.dif,border=T,height = unit(5, "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue1")
#dev.off()
ha=HeatmapAnnotation(NMF=as.character(nmfcluster[,2]),annotation_name_side = "left")
pdf(paste0(argvs[2],'immue_score_NMF.pdf'),width=9)
Heatmap(es.diforder,border=T,col =f1,height = unit(4*nrow(es.diforder), "mm"),show_column_names=F,top_annotation =ha,cluster_rows = FALSE,cluster_columns=FALSE,row_names_side = "left",name="immue")
dev.off()
# nmfcluster_Clindat <- fread('nmfcluster_Clindat_filters.xls')
# library(data.table)
# nmfcluster_Clindat <- fread('nmfcluster_Clindat_filters.xls')
# nmfcluster_Clindat <- fread('nmfcluster_Clindat.xls')
# nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,'class'] = 'Immune'
# nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,'class']
# nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,]
# dim(nmfcluster_Clindat[nmfcluster_Clindat$NMFconsushctree==7,])
# library(survival)
# my.surv <- Surv(as.numeric(nmfcluster_Clindat$survi_all_time),nmfcluster_Clindat$vital_status=='dead')
# kmfit3 <- survfit(my.surv~nmfclusterClindat11$class,data=nmfclusterClindat11)
# kmfit3 <- survfit(my.surv~nmfclusterClindat11$class,data=nmfcluster_Clindat)
# kmfit3 <- survfit(my.surv~nmfcluster_Clindat$class,data=nmfcluster_Clindat)
# cox.res <-coxph(Surv(as.numeric(nmfcluster_Clindat$survi_all_time)/30,nmfcluster_Clindat$status)~nmfcluster_Clindat$class)
# cox.res <-coxph(Surv(as.numeric(nmfcluster_Clindat$survi_all_time)/30,nmfcluster_Clindat$vital_status=='dead')~nmfcluster_Clindat$class)
# summary(cox.re)
# summary(cox.res)
# library(survminer)
# ggsurvplot(kmfit3)
# ggsurvplot(kmfit3,pval=T)
# write.table(nmfcluster_Clindat,file="nmf_2_7_immune_cluster_clindat.xls",sep="\t",quote=F)
# getwd()
# write.table(nmfcluster_Clindat,file="nmf_2_7_immune_cluster_clindat.xls",sep="\t",quote=F,row.names=F)
