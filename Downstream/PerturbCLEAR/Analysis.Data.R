library(qs)
library(openxlsx)
library(dplyr)
library(tidyr)

print("Load data")
obj=qread("Obj.Scala.qs") ##data from Scala et al
dat=read.csv("251208_all_metrics_filledBlanks.csv")
meta=read.csv("251121_metadata.csv")
met=read.csv("251208_MetricsSet1.csv")



tab=read.csv("251210_Scala_new_filledBlanks.csv")
rownames(tab)=tab[,1]
tab=tab[,3:ncol(tab)]
tab <- t(apply(tab, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))
obj[[2]]=tab #[,inter]

dat_all=dat

dat=dat[,met[,1]]

print("Cluster")
set.seed(1)
source('Cluster.R')
dat <- t(apply(dat, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))
res=ProcessData(meta=meta,dat=dat,exclude=c(),k=5)
qsave(res,"Cluster.results.comb.qs")
qsave(meta,"meta.comb.qs")
qsave(dat,"dat.comb.qs")
qsave(dat_all,"dat_all.comb.qs")

print("Plot")
library(ggplot2);library(cowplot);theme_set(theme_cowplot())
p=ggplot(res,aes(x=PC1,y=PC2,color=ClusterRaw))+geom_point()
ggsave("Cluster.orig.pdf",p)

p=ggplot(res,aes(x=PC1,y=PC2,color=ClusterPC))+geom_point()
ggsave("ClusterPC.orig.pdf",p)

print("Compare to Scala")
#mn=apply(dat1,2,function(x) sum(is.na(x)))
#dat1=dat1[,mn==0]
dat1=dat_all[,met[,1]]
dat1 <- t(apply(dat1, 1, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}))
COR=cor(t(scale(obj[[2]][,colnames(dat1)])),t(scale(dat1))) #correlation between Scala and our data
mn=apply(COR,2,which.max)
CT=obj[[3]][mn,"RNA.type"]
print(table(CT))
res["Scala"]=CT
p=ggplot(res,aes(x=UMAP1,y=UMAP2,color=Scala))+geom_point()+theme(legend.position="none")
ggsave("Cluster.Scala.pdf",p)

UseSVM<-function(dat,obj)
{
    library(e1071)
    inter=intersect(colnames(dat),colnames(obj[[2]]))
    dat=dat[,inter]
    dat <- t(apply(dat, 1, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    }))
    
    obj[[3]]["Clust"]=map_chr(obj[[3]]$RNA.type,function(x) strsplit(x,"_")[[1]][1])
    mod=svm(x=scale(obj[[2]][,inter]),y=factor(obj[[3]]$Clust),probability=T)
    pred=predict(mod,scale(dat[,inter]),probability=T)
    prob=attr(pred,"prob")
    print(head(prob))
    return(as.numeric(prob[,1]))
}




res["predSVM"]=UseSVM(dat,obj)

res["Label"]="L2/3 IT"
res[res$predSVM>.5,"Label"]="L4/5 IT"
p=ggplot(res,aes(x=PC1,y=PC2,color=Label))+geom_point()
ggsave("Cluster.Label.pdf",p)

qsave(res,"Cluster.results.comb.qs")
qsave(meta,"meta.comb.qs")
qsave(dat,"dat.comb.qs")
qsave(dat_all,"dat_all.comb.qs")

print("Test for perts")
res=res[meta$Type=="Pert",]
dat=dat[meta$Type=="Pert",]
dat_all=dat_all[meta$Type=="Pert",]
meta=meta[meta$Type=="Pert",]

print("Test for differences")
library(purrr)
out=map(c("L2/3 IT","L4/5 IT"),function(x){
    TestPert(scale(dat_all[res$Label==x,met[,1]]),res[res$Label==x,],form=Met~Perturbation+(1|Brain))
})
names(out)=c("L2/3 IT","L4/5 IT")
qsave(out,"Test.results.qs")
out2=out
names(out2)=c("L23_IT","L45_IT")
write.xlsx(out2,"Test.results.xlsx")

library(TRADE)
library(ashr)
out1=split(out[[1]],out[[1]]$Test)
val=map_dbl(out1,function(x){colnames(x)[2]="SE";colnames(x)[5]="P";TRADE(mode="univariate",results1=x,log2FoldChange="Estimate",lfcSE="SE",pvalue="P")$distribution_summary$transcriptome_wide_impact})
val1=val
dat1=data.frame("Perturbation"=names(val1),"TWI"=val1)
dat1["CT"]="L2/3"
out2=split(out[[2]],out[[2]]$Test)
val=map_dbl(out2,function(x){colnames(x)[2]="SE";colnames(x)[5]="P";TRADE(mode="univariate",results1=x,log2FoldChange="Estimate",lfcSE="SE",pvalue="P")$distribution_summary$transcriptome_wide_impact})
val2=val
dat2=data.frame("Perturbation"=names(val2),"TWI"=val2)
dat2["CT"]="L4/5"
dat=rbind(dat1,dat2)
dat["Gene"]=sub("Perturbation","",dat$Perturbation)
qsave(dat,"TRADE.TWI.byType.qs")
write.table(dat,"TRADE.TWI.byType.txt",sep="\t",row.names=F,quote=F)




print("Test for perts")
dat=qread("dat.comb.qs")
meta=qread("meta.comb.qs")
dat=dat[meta$Type=="Pert",]
meta=meta[meta$Type=="Pert",]

print("Test for differences")
library(purrr)
out2=map(unique(res$ClusterRaw),function(x){
    TestPert(scale(dat_all[res$ClusterRaw==x,met[,1]]),res[res$ClusterRaw==x,],form=Met~Perturbation+(1|Brain))
})
names(out2)=unique(res$ClusterRaw)

qsave(out2,"Test.results.cluster.qs")
write.xlsx(out2,"Test.results.cluster.xlsx")
print("Done!")

print("Get NMF")
set.seed(1234)
nmfRet=GetNMF(dat_all[res$Label=="L4/5 IT",met[,1]],5)
nmf=nmfRet[[1]]
nmfH=nmfRet[[2]]
nmfMrk=TestPert(nmf,res[res$Label=="L4/5 IT",],form=Met~Perturbation+Tracing+(1|Brain))
qsave(nmf,"NMF.L4_5IT.qs")
qsave(nmfRet,"NMFRet.L4_5IT.qs")
qsave(nmfMrk,"NMF.L4_5IT.Test.qs")
print(head(nmfMrk))
COR=cor(nmf,dat[res$Label=="L4/5 IT",met[,1]])
library(ComplexHeatmap)
pdf("Heatmap.NMF.L4_5IT.pdf",width=7,height=10)
Heatmap(t(COR),name="r")
dev.off()

nmfRet2=GetNMF(dat_all[res$Label=="L2/3 IT",met[,1]],5)
nmf2=nmfRet2[[1]]
nmf2H=nmfRet2[[2]]
nmfMrk2=TestPert(nmf2,res[res$Label=="L2/3 IT",],form=Met~Perturbation+(1|Brain))
qsave(nmf2,"NMF.L2_3IT.qs")
qsave(nmfMrk2,"NMF.L2_3IT.Test.qs")
qsave(nmfRet2,"NMFRet.L2_3IT.qs")

COR=cor(nmf2,dat[res$Label=="L2/3 IT",met[,1]])
library(ComplexHeatmap)
pdf("Heatmap.NMF.L2_3IT.pdf",width=7,height=10)
Heatmap(t(COR),name="r")
dev.off()


nmfRet=GetNMF(dat_all[,met[,1]],10)
nmf=nmfRet[[1]]
nmfH=nmfRet[[2]]
nmfMrk3=TestPert(nmf,res,form=Met~Perturbation+Tracing+(1|Brain))
qsave(nmf,"NMF.all.qs")
qsave(nmfMrk3,"NMF.all.Test.qs")
qsave(nmfRet,"NMFRet.all.qs")


COR=cor(nmf,dat[,met[,1]])
library(ComplexHeatmap)
pdf("Heatmap.NMF.pdf",width=7,height=10)
Heatmap(t(COR),name="r")
dev.off()

ret=list("L23_IT"=nmfMrk2,"L45_IT"=nmfMrk,"All"=nmfMrk3)
write.xlsx(ret,"NMF.Test.results.xlsx")

out=map(1:5,function(x){
    tab<-res %>% group_by(Perturbation,Brain) %>% summarise(Num=mean(ClusterRaw==x)) %>% as.data.frame()
    tab[1]=factor(tab[,1])
    tab[1]=relevel(tab[,1],ref="ctrl")
    fit=lm(Num~Perturbation,tab)
    coef=summary(fit)$coefficients
    print(coef)
    coef=data.frame(coef)
    coef["Pert"]=rownames(coef)
    coef=coef[grep("Perturbation",rownames(coef)),]
    coef["Cluster"]=x
    return(coef)
})
out=do.call(rbind,out)
colnames(out)[1:4]=c("Estimate","StdError","tvalue","pval")
out=out[order(out$pval),]
out["padj"]=p.adjust(out$pval,method="fdr")
qsave(out,"ChangeCluster.qs")

print("subcluster L2/3 IT")
res_sub=ProcessData(meta=meta[res$Label=="L2/3 IT",],dat=dat[res$Label=="L2/3 IT",],exclude=c(),k=2)
dat_sub=dat_all[res$Label=="L2/3 IT",]
qsave(res_sub,"Cluster.results.L2_3_IT.sub.qs")
print("plot it")
p=ggplot(res_sub,aes(x=PC1,y=PC2,color=ClusterRaw))+geom_point()
ggsave("Cluster.L2_3_IT.sub.pdf",p)
print("Test for differences between Perts in each cluster")
out=map(unique(res_sub$ClusterRaw),function(x){
    TestPert(scale(dat_sub[res_sub$ClusterRaw==x,met[,1]]),res_sub[res_sub$ClusterRaw==x,],form=Met~Perturbation+(1|Brain))
})
names(out)=unique(res_sub$ClusterRaw)
qsave(out,"Test.results.L2_3_IT.sub.qs")
write.xlsx(out,"Test.results.L2_3_IT.sub.xlsx")

res=qread("Cluster.results.comb.qs")
meta=qread("meta.comb.qs")
dat=qread("dat.comb.qs")
dat_all=qread("dat_all.comb.qs")

res=res[meta$Type!="Pert",]
dat=dat[meta$Type!="Pert",]
dat_all=dat_all[meta$Type!="Pert",]
meta=meta[meta$Type!="Pert",]



library(purrr)
out=map(unique(res$ClusterRaw),function(x){
    tryCatch({TestPert(scale(dat_all[res$ClusterRaw==x,met[,1]]),res[res$ClusterRaw==x,],form=Met~Age+(1|Brain),ctrl="P10",pertCol="Age")},error=function(e){NULL})
})
names(out)=unique(res$ClusterRaw)
out <- out[!sapply(out, is.null)]
for(i in names(out)){
    colnames(out[[i]])[5]="pval"
}

qsave(out,"Test.results.by.Cluster.ageData.qs")
write.xlsx(out,"Test.results.by.Cluster.ageData.xlsx")


out=map(unique(res$Label),function(x){
    tryCatch({TestPert(scale(dat_all[res$Label==x,met[,1]]),res[res$Label==x,],form=Met~Age+(1|Brain),ctrl="P10",pertCol="Age")},error=function(e){NULL})
})
names(out)=sub("/","_",unique(res$Label),fixed=T)
out <- out[!sapply(out, is.null)]
for(i in names(out)){
    colnames(out[[i]])[5]="pval"
}

qsave(out,"Test.results.by.Type.ageData.qs")
write.xlsx(out,"Test.results.by.Type.ageData.xlsx")

print("NMF")
nmf=GetNMF(dat_all[,met[,1]],5)[[1]]
print(head(nmf))
qsave(nmf,"NMF.age.qs")
nmfMrk=TestPert(nmf,res,form=Met~Age+(1|Brain),ctrl="P10",pertCol="Age")
qsave(nmf,"NMF.age.qs")
qsave(nmfMrk,"NMF.age.qs")


res=qread("Cluster.results.comb.qs")
meta=qread("meta.comb.qs")
dat=qread("dat.comb.qs")



print("Map NMF to RNA")
dat=obj[[1]]
meta=obj[[3]]
tab=obj[[2]]
tab=tab[meta$RNA.type!="L2/3 IT_3",]
dat=dat[,meta$RNA.type!="L2/3 IT_3"]
dat=data.frame(as.matrix(dat))
for(i in colnames(dat)){
    dat[i]=1000000*dat[,i]/sum(dat[,i])
}
nmfRet=qread("NMFRet.L4_5IT.qs")
w=nmfRet[[2]]
tab=tab[,colnames(w)]
h=RcppML::project(data=t(tab),w=w)
mn=apply(dat,1,mean)
dat=dat[mn>1,]
inter=intersect(rownames(dat),res$Gene)
w2=RcppML::project(data=scale(t(log(dat[inter,]+1)),center=F),w=h);
w2=t(w2);
w2=data.frame(w2);
w2["Gene"]=rownames(w2);
w2=w2[order(w2$nmf3,decreasing=T),]



