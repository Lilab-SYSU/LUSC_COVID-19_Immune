library(data.table)
args <- commandArgs(T)
TCGARNAfiread <- function(dir){
fis <- dir(dir)
fis <- grep('FPKM',fis,value = T)
data <- fread(fis[1],sep="\t",header=F)
print(fis[1])
colnames(data) <-c('Genes',fis[1])
for(i in fis[2:length(fis)]){
fir <- fread(i,sep="\t",header=F)
colnames(fir) <- c('Genes',i)
data = merge(data,fir,by="Genes")
}
return(data)
}
datout <- TCGARNAfiread(args[1])###give a dir

names <- args[2]### project name
colnames(datout) <- gsub('\\.FPKM.*gz','',colnames(datout))

fwrite(datout,file=paste0(names,"_TCGA_ALL_Samples.txt"))
Genes <- gsub('\\..*','',datout$Genes)

library(rtracklayer)
library(SummarizedExperiment)
#setwd('~/DB/dog/gtf')
gtf1 <- rtracklayer::import('/home/myang/DB/TCGA/gencode.v22.annotation.gtf.gz')
gtf2 <- as.data.frame(gtf1[gtf1$type=="gene",])

loc <- match(Genes,gsub('\\..*','',gtf2$gene_id))
annotInfo <- gtf2[loc,c('gene_name','gene_id','seqnames','start','end','strand','gene_type')]
genesName <- annotInfo$gene_name
datout1 <- data.frame(Genes=genesName,datout[,-1])
colnames(datout1) <- gsub('\\.FPKM.*gz','',colnames(datout))
datlincRNA = datout1[annotInfo$gene_type=="lincRNA",]
datprotein = datout1[annotInfo$gene_type=="protein_coding",]
datother = datout1[annotInfo$gene_type !="protein_coding" & annotInfo$gene_type !="lincRNA",]
#colnames(datout1) <- gsub('\\.FPKM.*gz','',colnames(datout))
fwrite(datout1,file=paste0(names,"_Symbols_TCGA_ALL_Samples.txt"))
fwrite(datlincRNA,file=paste0(names,"_lincRNA_Symbols_TCGA_ALL_Samples.txt"))
fwrite(datprotein,file=paste0(names,"_protein_Symbols_TCGA_ALL_Samples.txt"))
fwrite(datother,file=paste0(names,"_other_Symbols_TCGA_ALL_Samples.txt"))




