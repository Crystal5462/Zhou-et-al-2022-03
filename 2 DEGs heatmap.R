######

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

library(limma)
library(pheatmap)

setwd("D:\\bio\\diff")       #set working direcory
inputFile="symbol.TMB.txt"                                   #input file
fdrFilter=0.05                                               #fdr cut-off 
logFCfilter=1                                                #logFC cut-off
conNum=184                                                   #the number of low TMB samples
treatNum=183                                                 #the number of high TMBsamples

#read input file
outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.2,]

#differential analysis
for(i in row.names(data)){
	  geneName=unlist(strsplit(i,"\\|",))[1]
	  geneName=gsub("\\/", "_", geneName)
	  rt=rbind(expression=data[i,],grade=grade)
	  rt=as.matrix(t(rt))
	  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
	  conGeneMeans=mean(data[i,1:conNum])
	  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	  logFC=log2(treatGeneMeans)-log2(conGeneMeans)  
	  pvalue=wilcoxTest$p.value
	  conMed=median(data[i,1:conNum])
	  treatMed=median(data[i,(conNum+1):ncol(data)])
	  diffMed=treatMed-conMed
	  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
			outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

#output the diffenrences in all genes 
write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)
write.table(outDiff,file="diff.txt",sep="\t",row.names=F,quote=F)

#output the expression of differential genes
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)

#heatmap of differential genes
geneNum=20
diffSig=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(diffSig[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(geneNum*2) ){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
hmExp=log2(hmExp+0.001)
Type=c(rep("low",conNum),rep("high",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=4,width=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8)
dev.off()


######