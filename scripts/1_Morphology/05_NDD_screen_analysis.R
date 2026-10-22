# Perturb-CLEAR NDD screen analysis.

suppressPackageStartupMessages({
  library(qs)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(purrr)
  library(e1071)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_directory <- if (length(script_argument) == 1L) {
  dirname(normalizePath(sub("^--file=", "", script_argument), winslash = "/"))
} else {
  getwd()
}
cluster_script <- list.files(
  script_directory,
  pattern = "^04_Morphology_clustering[.]R$",
  full.names = TRUE
)
if (length(cluster_script) != 1L) {
  stop("Expected 04_Morphology_clustering.R in the script directory.")
}

# Command-line arguments, in order:
# 1. Scala reference object in qs format.
# 2. Perturb-CLEAR morphometric measurements CSV.
# 3. Perturb-CLEAR cell-metadata CSV.
# 4. Prespecified 36-feature list CSV.
# 5. Scala morphometric measurements CSV.
# 6. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Expected 6 command-line arguments in the order documented above.")
}
scala_object_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
morphometrics_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
metadata_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
feature_file <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
scala_morphometrics_file <- normalizePath(args[[5]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[6]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
starting_directory <- getwd()
setwd(output_dir)
on.exit(setwd(starting_directory), add = TRUE)

print("Load data")
obj=qread(scala_object_file) ##data from Scala et al.
dat=read.csv(morphometrics_file)
meta=read.csv(metadata_file)
met=read.csv(feature_file)



tab=read.csv(scala_morphometrics_file)
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
source(cluster_script)
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
theme_set(theme_cowplot())
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
out=map(c("L2/3 IT","L4/5 IT"),function(x){
    TestPert(scale(dat_all[res$Label==x,met[,1]]),res[res$Label==x,],form=Met~Perturbation+(1|Brain))
})
names(out)=c("L2/3 IT","L4/5 IT")
qsave(out,"Test.results.qs")
out2=out
names(out2)=c("L23_IT","L45_IT")
write.xlsx(out2,"Test.results.xlsx")





print("Test for perts")
dat=qread("dat.comb.qs")
meta=qread("meta.comb.qs")
dat=dat[meta$Type=="Pert",]
meta=meta[meta$Type=="Pert",]

print("Test for differences")
out2=map(unique(res$ClusterRaw),function(x){
    TestPert(scale(dat_all[res$ClusterRaw==x,met[,1]]),res[res$ClusterRaw==x,],form=Met~Perturbation+(1|Brain))
})
names(out2)=unique(res$ClusterRaw)

qsave(out2,"Test.results.cluster.qs")
write.xlsx(out2,"Test.results.cluster.xlsx")
print("Done!")


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

