source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("cummeRbund")


library(cummeRbund)
setwd("/Volumes/DOWNLOAD/CMD/CMD_cuffdiff")

cuff<-readCufflinks()
cuff

sample.names<-samples(genes(cuff))
head(sample.names, n=20L)

gene.diff<-diffData(genes(cuff))
head(gene.diff,n=30L)

s<-csScatter(genes(cuff),"A","E",smooth=T)
s<-csScatter(isoforms(cuff),"A","E",smooth=T)
s

dend<-csDendro(isoforms(cuff))
dend

dend.rep<-csDendro(isoforms(cuff),replicates=T)
dend.rep

b<-csBoxplot(genes(cuff))
b

brep<-csBoxplot(genes(cuff),replicates=T)
brep

myProfile<-c(500,0,500,0,0,500,0)
mySimilar2<-findSimilar(cuff,myProfile,n=20)
mySimilar2.expression<-expressionPlot(mySimilar2,logMode=T,showErrorbars=F)
mySimilar2.expression

mySimilar<-findSimilar(cuff,"Manes.17G079800",n=20)
mySimilar.expression<-expressionPlot(mySimilar,logMode=T,showErrorbars=F)
mySimilar.expression

A_vs_E.sigGenesIds<-getSig(cuff,x='A',y='E',alpha=0.05,level='isoforms')
A_vs_E_SigGenes<-getGenes(cuff,A_vs_E.sigGenesIds)

ic<-csCluster(A_vs_E_SigGenes,k=12)
icp<-csClusterPlot(ic)
icp