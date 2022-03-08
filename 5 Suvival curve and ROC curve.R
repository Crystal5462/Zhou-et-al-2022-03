######

#install.packages("survival")
#install.packages("survivalROC")

library(survival)
library(survivalROC)
setwd("D:\\biowolf\\TMBimmune\\26.survival")       #set working directory

rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
#plot the survival curve
pdf(file="survival.pdf", width=5, height=5)
plot(fit, lwd=2, col=c("red","blue"),mark.time=T,
     xlab="Time (year)", ylab="Survival rate",
     main=paste0("Survival curve(",pValue,")") )
legend("topright", lwd=2, 
     c("High risk","Low risk"), col=c("red","blue"))
dev.off()
summary(fit)    #five-year survival rates


#plot the ROC curve
pdf(file="ROC.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
      predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
    xlab="False positive rate", ylab="True positive rate",
    main=paste("ROC curve (", "AUC = ",sprintf("%.3f",roc$AUC),")"),
    lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()


######