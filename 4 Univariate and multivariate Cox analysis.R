######

#install.packages("survival")

library(survival)

setwd("D:\\bio\\uniCox")                    #set working directory
pFilter=0.05                                                                #filter condition
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)    #read the input file
rt$futime=rt$futime/365
#rt[,3:ncol(rt)]=log2(rt[,3:ncol(rt)]+1)

outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
	   if(sd(rt[,gene])<0.01){next}
	   if(grepl("-", gene)){next}
	   cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
	   coxSummary = summary(cox)
	   coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	   if(coxP<pFilter){
		   	   group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
		       diff=survdiff(Surv(futime, fustat) ~group,data = rt)
		       pValue=1-pchisq(diff$chisq,df=1)
		       if(pValue<pFilter){
			       sigGenes=c(sigGenes,gene)
			       outTab=rbind(outTab,
			                    cbind(gene=gene,
			                         #KM=pValue,
			                         HR=coxSummary$conf.int[,"exp(coef)"],
			                         HR.95L=coxSummary$conf.int[,"lower .95"],
			                         HR.95H=coxSummary$conf.int[,"upper .95"],
					                 coxPvalue=coxP) )
				}
		}
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)    #output genes and p-value
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


######
#install.packages('survival')

library(survival)                                         #packages
setwd("D:\\bio\\multiCox")       #set working directory
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)    #read input file

#COX regression model
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#output model parameters
outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=outTab[,c(1,2)]
outTab=gsub("`","",outTab)
write.table(outTab,file="coef.txt",sep="\t",row.names=F,quote=F)

#output the risk scores of each patient
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
    row.names=F)

for(gene in row.names(outTab)){
	group=ifelse(rt[,gene]>median(rt[,gene]),"high","low")
	diff=survdiff(Surv(futime, fustat) ~group,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	#plot the survival curve
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
	pdf(file=paste0("survival.",gene,".pdf"), width=5, height=5)
	plot(fit, 
		lwd=2,
		col=c("red","blue"),
		xlab="Time (year)",
		mark.time=T,
		ylab="Survival rate",
		main=paste0(gene,"(",pValue,")") )
	legend("topright", 
			c("High expression","Low expression"), 
			lwd=2, 
			col=c("red","blue"))
	dev.off()
}


######
