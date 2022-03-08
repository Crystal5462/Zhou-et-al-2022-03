######

#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)           #package
setwd("D:\\maftools")      #set working directory
maf = read.maf(maf = 'input.maf')          #read inpt file

#summary figure
pdf(file="summary.pdf",width=7,height=6)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#waterfall plot
pdf(file="waterfall.pdf",width=7,height=6)
oncoplot(maf = maf, top = 30, fontSize = 0.8 ,showTumorSampleBarcodes = F )
dev.off()

#corplot
pdf(file="interaction.pdf",width=7,height=6)
somaticInteractions(maf = maf, top = 20, pvalue = c(0.05, 0.001),showCounts = FALSE, fontSize = 0.6)
dev.off()

#gene cloud figure
pdf(file="Genecloud.pdf",width=7,height=6)
geneCloud(input = maf, minMut = 30)
dev.off()


######
######

