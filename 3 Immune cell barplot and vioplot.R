######

setwd("D:\\bio\\barplot")    #set working directory
input="CIBERSORT-Results.txt"       #input file
outpdf="barplot.pdf"                #output file name
pFilter=0.05                        #CIBERSORT result filter condition

#read the input file and organize the data
immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
data=t(immune)
col=rainbow(nrow(data),s=0.7,v=0.7)

#draw the barplot
pdf(outpdf,height=10,width=22)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()


######
######

#install.packages("vioplot")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#quote the package
library(vioplot)
library(limma)

pFilter=0.05               #CIBERSORT results filter condition
setwd("D:\\bio\\vioplot")       #set working directory

#read the result of CIBERSORT and organize the data
immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(immune))
immune=avereps(immune)

#read the TMB value and organize the data
tmb=read.table("TMB.txt",sep="\t",header=T,check.names=F,row.names=1)
tmb=as.matrix(tmb)
row.names(tmb)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(tmb))
tmb=avereps(tmb)
lowTmb=tmb[tmb[,"TMB"]<=median(tmb[,"TMB"]),]
highTmb=tmb[tmb[,"TMB"]>median(tmb[,"TMB"]),]
lowTmbName=names(lowTmb)
highTmbName=names(highTmb)

#extract the immune cell content of high and ow TMB groups
lowTmbImm=intersect(row.names(immune),lowTmbName)
highTmbImm=intersect(row.names(immune),highTmbName)
rt=rbind(immune[lowTmbImm,],immune[highTmbImm,])
lowTmbNum=length(lowTmbImm)
highTmbNum=length(highTmbImm)

#output the vioplot
pdf("vioplot.pdf",height=8,width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#vioplot is plotted for the cycle of each immune cells, with low TMB in green and high TMB in red
for(i in 1:ncol(rt)){
	  if(sd(rt[1:lowTmbNum,i])==0){
	    rt[1,i]=0.001
	  }
	  if(sd(rt[(lowTmbNum+1):(lowTmbNum+highTmbNum),i])==0){
	    rt[(lowTmbNum+1),i]=0.001
	  }
	  lowTmbData=rt[1:lowTmbNum,i]
	  highTmbData=rt[(lowTmbNum+1):(lowTmbNum+highTmbNum),i]
	  vioplot(lowTmbData,at=3*(i-1),lty=1,add = T,col = 'green')
	  vioplot(highTmbData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
	  wilcoxTest=wilcox.test(lowTmbData,highTmbData)
	  p=wilcoxTest$p.value
	  mx=max(c(lowTmbData,highTmbData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()


######
