library(NMF)
###brunet lee ns
argv <- commandArgs(T)
library(data.table)
dat<- fread(argv[1])
dat <- as.data.frame(dat)
geneidfactor<-factor(dat[,1])
gene_exp_matrix<-apply(dat[,-1],2,function(x) tapply(x,geneidfactor,mean))
gene_exp_matrix_filter <- gene_exp_matrix[rowMeans(gene_exp_matrix)>1,]
#load('.RData')
method = 'brunet'
estim.r <- nmf(gene_exp_matrix_filter, 2:10,method, nrun=40, seed=123456)
save.image(paste0('NMF',method,'.RData'))
pdf(paste0('NMF_',method,'_measure.pdf'))
plot(estim.r)
dev.off()
pdf(paste0('NMF_',method,'_heatmap.pdf'))
consensusmap(estim.r)
dev.off()
