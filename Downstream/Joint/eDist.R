library(scperturbR)
library(sva)
library(Seurat)
library(purrr)
library(ComplexHeatmap)
library(tidyr)
library(dplyr)



get_edist_morpho=function(dat,meta,batch="Brain",pert="Perturbation",minCells=5,verbose=T,sample_correction=F)
{
    vals=table(meta[,pert])
    nams=names(vals)[vals>minCells]
    dat=dat[meta[,pert] %in% nams,]
    meta=meta[meta[,pert] %in% nams,]
    dat=scale(dat)
    mod=model.matrix(~as.factor(meta[,pert]))
    
    tab=ComBat(t(dat),as.factor(meta[,batch]),mod=mod,par.prior=TRUE,prior.plots=FALSE)
    emb=t(tab)
    

    ##From edist command in scperturbR
    print("edist")
    labels <- meta[,pert]
    groups <- unique(labels)
    df <- setNames(data.frame(matrix(ncol = length(groups), nrow = length(groups)), 
        row.names = groups), groups)
    if (verbose) {
        print("Computing E-test statistics for each group.")
    }
    completed_groups <- c()
    for (groupx in groups) {
        for (groupy in groups) {
            if (groupy %in% completed_groups) {
                next
            }
            x <- as.matrix(emb)[labels == groupx, ]
            y <- as.matrix(emb)[labels == groupy, ]
            N <- nrow(x)
            M <- nrow(y)
            dist_xy <- rdist::cdist(x, y)
            dist_x <- rdist::pdist(x)
            dist_y <- rdist::pdist(y)
            if (sample_correction) {
                ed <- 2 * (sum(dist_xy)/(N * M)) - (sum(dist_x)/(N * 
                  (N - 1))) - (sum(dist_y)/(M * (M - 1)))
            }
            else {
                ed <- 2 * (sum(dist_xy)/(N * M)) - (sum(dist_x)/(N * 
                  N)) - (sum(dist_y)/(M * M))
            }
            df[groupx, groupy] <- df[groupy, groupx] <- ed
        }
        completed_groups <- c(completed_groups, groupx)
    }
    return(df)
    ##End of code from edist command in scperturbR

    
}

get_edist_rna<-function(seur,batch="Sample",pert="Assign",dims=1:25,minCells=90)
{
    #seur=RunHarmony(seur,group.by.vars = batch,reduction = "pca",dims = dims)
    vals=table(seur@meta.data[,pert])
    nams=names(vals)[vals>minCells]
    seur=subset(seur,cells=names(seur@active.ident)[seur@meta.data[,pert] %in% nams])
    df=edist(seur,pert,reduction = "harmony")
    #vals=table(seur@meta.data[,pert])
    #nams=names(vals)[vals>minCells]
    #df=df[nams,nams]
    nams=colnames(df)
    nams[nams=="ST_2"]="ctrl"
    nams=map_chr(nams,function(x) strsplit(x,"_")[[1]][1])
    colnames(df)=nams
    rownames(df)=nams
    nams=nams[grep("ST",nams,invert=T)]
    nams=nams[grep("NT",nams,invert=T)]
    return(df[nams,nams])
}


RunPip<-function(type="L2/3")
{
    print("Run for RNA")
    seurfil=paste0("/stanley/levin_xinjin/JinLab/ssimmons/BoLi/CASTCombined/seur.",sub("/","",type,fixed=T),".harmony.RDS")
    lab=paste0(type," IT")

    seur=readRDS(seurfil)
    df1=get_edist_rna(seur,batch="Sample",pert="Assign",dims=1:30,minCells=90)

    print("Run for morpho")
    res=readRDS("/stanley/levin_xinjin/JinLab/ssimmons/BoLi/CASTCombined/Morpho/Data_Jan14_26/Cluster.results.comb.RDS")
    res=res[res$Perturbation!="None" & res$Label==lab,]
    dat=res[,13:48]
    dat=scale(dat)
    meta=res
    df2=get_edist_morpho(dat,meta,batch="Brain",pert="Perturbation",minCells=5,verbose=T,sample_correction=F)

    print("Combine")
    
    ret=list("RNA"=df1,"Morpho"=df2)
    inter=intersect(colnames(df1),colnames(df2))
    df1=df1[inter,inter]
    df2=df2[inter,inter]
    df1=df1/max(df1)
    df2=df2/max(df2)
    ret[["Combined_Ave"]] = (df1[inter,inter]+df2[inter,inter])/2
    df=(df1[inter,inter]+df2[inter,inter])/2
    clust=hclust(as.dist(df))
    inter=inter[clust$order]

    #inter=intersect(colnames(df1),colnames(df2))
    
    df=df1[inter,inter]

  
    for(i in 1:(length(inter)-1))
    {
        for(j in (i+1):length(inter))
        {
            print(i)
            print(j)
            df[i,j]=df2[inter[i],inter[j]]
        }
    }

    ret[["Combined"]] = df
    ret[["inter"]]=inter
    return(ret)

}

RunL23<-function()
{
    return(RunPip())
    
}

RunL45<-function()
{
    return(RunPip("L4/5"))
}

PlotHeatmap<-function(ret)
{
    df=ret$Combined
    #clust=hclust(as.dist(ret$Combined_Ave))
    Heatmap(df,cluster_rows = F,cluster_columns = F)
}


metaAnalysis_edist<-function(seur,batch="Sample",pert="Assign",minCells=90,ctrl="ST_2")
{
    vals=table(seur@meta.data$Assign)
    seur=subset(seur,cells=names(seur@active.ident)[seur@meta.data$Assign %in% names(vals)[vals>minCells]])
    res=map(unique(seur@meta.data$Sample),function(x){tryCatch({out=etest(subset(seur,Sample==x),"Assign",ctrl,"harmony");out=out[order(out[,1]),];out["FDR"]=p.adjust(out[,1],"fdr");out["Pert"]=rownames(out);out},error=function(y){return(NULL)})})
    res=do.call(rbind,res)

    GetPval=function(ps){qs=log(ps);chi=-2*sum(qs);1-pchisq(chi,2*length(ps),F)}

    tab=res %>% group_by(Pert) %>% summarise(pval=GetPval(pval),edist=mean(edist),numTest=length(edist)) %>% data.frame();tab=tab[order(tab$pval),];tab["fdr"]=p.adjust(tab$pval,"fdr")
    return(tab)

}

etest_morpho_easy<-function(type="L2/3")
{
    res=readRDS("/stanley/levin_xinjin/JinLab/ssimmons/BoLi/CASTCombined/Morpho/Data_Jan14_26/Cluster.results.comb.RDS")
    lab=paste0(type," IT")
    res=res[res$Perturbation!="None" & res$Label==lab,]
    dat=res[,13:48]
    dat=scale(dat)
    meta=res
    res2=etest_morpho(dat,meta,batch="Brain",pert="Perturbation",ctrl="ctrl",minCells=5)
    return(res2)
}


etest_morpho<-function(dat,meta,batch="Brain",pert="Perturbation",ctrl="ctrl",minCells=5)
{
    vals=table(meta[,pert])
    nams=names(vals)[vals>minCells]
    dat=dat[meta[,pert] %in% nams,]
    meta=meta[meta[,pert] %in% nams,]
    dat=scale(dat)
    mod=model.matrix(~as.factor(meta[,pert]))
    
    tab=ComBat(t(dat),as.factor(meta[,batch]),mod=mod,par.prior=TRUE,prior.plots=FALSE)
    emb=t(tab)

    seur <- CreateSeuratObject(counts = t(emb), min.cells = 0, min.features =0)
    for(i in colnames(meta))
    {
        seur@meta.data[i]=meta[,i]
    }
    emb=scale(emb)
    colnames(emb)=paste0("emb_",1:ncol(emb))
    seur[["emb"]]=CreateDimReducObject(embeddings = emb, key = "emb_", assay = DefaultAssay(seur))
    res=etest(seur,pert,ctrl,reduction="emb")
    res=res[order(res[,1]),]
    res["padj"]=p.adjust(res$pval,"fdr")
    return(res)
}