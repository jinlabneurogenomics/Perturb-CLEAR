library(umap)
library(purrr)
library(lmerTest)
library(RcppML)

ProcessData<-function(meta=NULL,dat=NULL,comb=NULL,exclude=c("y-extent.(um)","Trunk.length.(um)"),k=5,nPCs=5,CVcutoff=-1,SkipUMAP=F)
{

    if(!is.null(comb))
    {
        meta=comb[,1:8]
        dat=comb[,9:dim(comb)[2]]
        dat=dat[,colnames(dat)!="cell_id"]
    }
    if(CVcutoff>0)
    {
        mn=apply(dat,2,function(x) sd(x)/mean(x))
        dat=dat[,mn>CVcutoff]
    }
    #dat[is.na(dat)]=0
    dat=scale(dat)
    dat[is.na(dat)]=0
    print("Cluster no exclude")
    out3=kmeans(scale(dat),k)

    #dat=dat[,!(colnames(dat) %in% exclude)]

    print("Run pca")
    pc=prcomp(t(scale(dat)))# ,center=F,scale=F)
    print((pc$sdev)^2/sum((pc$sdev)^2))
    rot=pc$rotation
    print(dim(rot))

    print("Cluster--PC")
    out1=kmeans(rot[,1:nPCs],k)

    print("Cluster--not PC")
    out2=kmeans(scale(dat),k)

    
    if(SkipUMAP)
    {
        comb=cbind(meta,dat,rot)
        comb["ClusterPC"]=as.character(out1$cluster)
        comb["ClusterRaw"]=as.character(out2$cluster)
    

        return(comb)
    }
    print("UMAP")
    umap <- umap(scale(dat))$layout
    umap=data.frame(umap)
    colnames(umap)=c("UMAP1","UMAP2")
    print(dim(umap))

    umap2 <- umap(rot[,1:nPCs])$layout
    umap2=data.frame(umap2)
    colnames(umap2)=c("UMAPPC1","UMAPPC2")
    print(dim(umap2))

    print("Combine and return")
    comb=cbind(meta,dat,umap,umap2,rot)
    comb["ClusterPC"]=as.character(out1$cluster)
    comb["ClusterRaw"]=as.character(out2$cluster)
    #comb["ClusterRawNoExclude"]=as.character(out3$cluster)

    return(comb)


}


TestComposition<-function(comb,clusterCol="ClusterRaw",pertCol="Perturbation",brainCol="Brain",ref="ctrl",useBrain=T)
{
    print("Prep data")
    tab=comb[,c(clusterCol,pertCol,brainCol)]
    colnames(tab)=c("Cluster","Pert","Brain")
    tab["Pert"]=factor(tab[,"Pert"])
    tab["Pert"]=relevel(tab[,"Pert"],ref=ref)
    ##List perts
    perts=table(tab[,"Pert"])
    perts=names(perts)[perts>5]
    #List clusts
    clusts=table(tab[,"Cluster"])
    clusts=names(clusts)[clusts>10]

    out=map(clusts,function(clust){
        tab["Lab"]=tab[,"Cluster"]==clust
        form=Lab~Pert+Brain
        if(!useBrain){form=Lab~Pert}
        fit=glm(form,data=tab,family="binomial")
        coef=data.frame(summary(fit)$coef)
        coef["Test"]=rownames(coef)
        colnames(coef)=c("Estimate","SE","Z","pval","Test")
        coef=coef[grep("^Pert",rownames(coef)),]
        rownames(coef)=NULL
        coef["Cluster"]=clust
        return(coef)
    })

    mrk=do.call(rbind,out)
    mrk=mrk[order(mrk$pval),]
    mrk["padj"]=p.adjust(mrk$pval,"fdr")
    return(mrk)

}



TestChanges<-function(comb,covCols,clusterCol="ClusterPC",pertCol="Perturbation",brainCol="X...Brain",ref="ctrl",useBrain=T)
{
    print("Prep data")
    tab=comb[,c(clusterCol,pertCol,brainCol)]
    colnames(tab)=c("Cluster","Pert","Brain")
    tab["Pert"]=factor(tab[,"Pert"])
    tab["Pert"]=relevel(tab[,"Pert"],ref=ref)
    ##List perts
    perts=table(tab[,"Pert"])
    perts=names(perts)[perts>5]
    #List clusts
    clusts=table(tab[,"Cluster"])
    clusts=names(clusts)[clusts>10]

    out=map(covCols,function(clust){
        tab["Lab"]=comb[,clust]
        form=Lab~Pert+Brain
        if(!useBrain){form=Lab~Pert}
        fit=lm(form,data=tab)
        coef=data.frame(summary(fit)$coef)
        coef["Test"]=rownames(coef)
        colnames(coef)=c("Estimate","SE","Z","pval","Test")
        coef=coef[grep("^Pert",rownames(coef)),]
        rownames(coef)=NULL
        coef["Cluster"]=clust
        return(coef)
    })

    mrk=do.call(rbind,out)
    mrk=mrk[order(mrk$pval),]
    mrk["padj"]=p.adjust(mrk$pval,"fdr")
    return(mrk)

}


TestPert<-function(dat,meta,form=Met~Perturbation+Tracing+(1|Brain),ctrl="ctrl",pertCol="Perturbation")
{
    meta[pertCol]=factor(meta[,pertCol])
    meta[pertCol]=relevel(meta[,pertCol],ref=ctrl)

    mets=colnames(dat)
    out=map(mets,function(met){
        tryCatch({
            print(paste0("Testing ",met))
            tab=meta
            tab["Met"]=dat[,met]
            #tab[is.na(tab$Met),"Met"]=mean(tab$Met,na.rm=T) ##Remove
            fit=lmer(form,data=tab[!is.na(tab$Met),])
            coef=data.frame(summary(fit)$coefficients)
            coef["Test"]=rownames(coef)
            coef["Metric"]=met
            coef=coef[grep(pertCol,rownames(coef)),]
            rownames(coef)=NULL
            return(coef)
        },error=function(e){NULL})
    })

    mrk=do.call(rbind,out)
    print(head(mrk))
    mrk=mrk[order(mrk[,5]),]
    mrk["padj"]=p.adjust(mrk[,5],"fdr")
    pert=unique(mrk$Test)
    for(i in pert)
    {
        mrk[mrk$Test==i,"padj_byTest"]=p.adjust(mrk[mrk$Test==i,5],"fdr")
    }
    return(mrk)

}


TestPertLM<-function(dat,meta,form=Met~Perturbation+Tracing+Brain,ctrl="ctrl",pertCol="Perturbation")
{
    meta[pertCol]=factor(meta[,pertCol])
    meta[pertCol]=relevel(meta[,pertCol],ref=ctrl)

    mets=colnames(dat)
    out=map(mets,function(met){
        tryCatch({
            print(paste0("Testing ",met))
            tab=meta
            tab["Met"]=dat[,met]
            #tab[is.na(tab$Met),"Met"]=0
            fit=lm(form,data=tab[!is.na(tab$Met),])
            coef=data.frame(summary(fit)$coefficients)
            coef["Test"]=rownames(coef)
            coef["Metric"]=met
            coef=coef[grep(pertCol,rownames(coef)),]
            rownames(coef)=NULL
            return(coef)
        },error=function(e){NULL})
    })

    mrk=do.call(rbind,out)
    print(head(mrk))
    mrk=mrk[order(mrk[,4]),]
    mrk["padj"]=p.adjust(mrk[,4],"fdr")
    pert=unique(mrk$Test)
    for(i in pert)
    {
        mrk[mrk$Test==i,"padj_byTest"]=p.adjust(mrk[mrk$Test==i,4],"fdr")
    }
    return(mrk)

}




TestPertRNA<-function(dat,meta,form=Met~Assign+(1|orig.ident),ctrl="NT_2",pertCol="Assign")
{
    meta[pertCol]=factor(meta[,pertCol])
    meta[pertCol]=relevel(meta[,pertCol],ref=ctrl)

    mets=colnames(dat)
    out=map(mets,function(met){
        print(paste0("Testing ",met))
        tab=meta
        tab["Met"]=dat[,met]
        fit=lmer(form,data=tab)
        coef=data.frame(summary(fit)$coefficients)
        coef["Test"]=rownames(coef)
        coef["Metric"]=met
        coef=coef[grep(pertCol,rownames(coef)),]
        rownames(coef)=NULL
        return(coef)
    })

    mrk=do.call(rbind,out)
    print(head(mrk))
    mrk=mrk[order(mrk[,5]),]
    mrk["padj"]=p.adjust(mrk[,5],"fdr")
    pert=unique(mrk$Test)
    for(i in pert)
    {
        mrk[mrk$Test==i,"padj_byTest"]=p.adjust(mrk[mrk$Test==i,5],"fdr")
    }
    return(mrk)

}




GetNMF<-function(tab,k=8)
{
    #print(head(scale(tab,center=F)))
    #tab[is.na(tab)]=0
    tab=scale(tab,center=F)
    print(head(tab))
    for(i in colnames(tab)){
        tab[is.na(tab[,i]),i]=mean(tab[,i],na.rm=T)
    }
    
    model <- RcppML::nmf(tab, k = k,seed=1) #,mask = "NA") ##k chosen with crossValidate command
    W=model@w
    H=model@h
    return(list(W=W,H=H))
}

CombNMF<-function(tab1,tab2,k=8)
{
    inter=intersect(colnames(tab1),colnames(tab2))
    tab1=tab1[,inter]
    tab2=tab2[,inter]
    model <- RcppML::nmf(scale(tab1,center=F), k = k,seed=1)
    h=model@h
    w1=model@w

    w2=RcppML::project(data=t(scale(tab2,center=F)),w=h)

    return(list(w1=w1,w2=t(w2),h=t(h)))


}



RNATest<-function(mrk,dat,cutoff=.05,sign=1)
{
    genes=mrk$Gene[mrk$adj_pval<cutoff & mrk$lfc*sign>0]
    #genes=mrk$Gene[mrk$lfc*sign>0]
    print(length(genes))
    #genes=head(genes,100)
    dat=data.frame(as.matrix(dat))
    for(i in colnames(dat)){
        dat[i]=10000*dat[i]/sum(dat[i])
    }
    dat=t(scale(t(dat)))
    print(length(genes))
    print(genes)
    print(sum(genes %in% rownames(dat)))
    print(head(rownames(dat)))
    dat=dat[rownames(dat) %in% genes,]
    mn1=apply(dat,1,mean)
    dat=dat[!is.nan(mn1),]
    ##print(dat)
    print(dim(dat))
    mn=apply(dat,2,mean)
    return(mn)

}


CompareRNAMorpho<-function(mrk,w,dat,cutoff=.05,method="spearman",sign=1)
{
    mn=RNATest(mrk,dat,cutoff=cutoff,sign=sign)
    Tests=map(1:dim(w)[2],function(i){
        cor.test(mn,w[,i],method=method)
    })
    pval=map_dbl(Tests,function(x) x$p.value)
    cors=map_dbl(Tests,function(x) x$estimate)
    ret=data.frame("pvalue"=pval,"rho"=cors,"Factor"=colnames(w))
    ret["FDR"]=p.adjust(ret$pvalue,"fdr")
    ret=ret[order(ret$pvalue),]
    lst=list("Cor"=ret,"MeanRNA"=mn,"W"=w)
    return(lst)
}


CompareRNAMorpho_byGene<-function(genes,w,dat,method="spearman")
{
    print(length(intersect(genes,rownames(dat))))
    res=map(genes,function(x){
        tryCatch({mn=as.numeric(dat[x,])
            print(head(mn))
            Tests=map(1:dim(w)[2],function(i){
                cor.test(mn,w[,i],method=method)
            })
            pval=map_dbl(Tests,function(x) x$p.value)
            cors=map_dbl(Tests,function(x) x$estimate)
            ret=data.frame("pvalue"=pval,"rho"=cors,"Factor"=colnames(w))
            ret["FDR"]=p.adjust(ret$pvalue,"fdr")
            ret=ret[order(ret$pvalue),]
            ret["Gene"]=x
            lst=ret
            print("Yay!")
            return(lst)
        }, error=function(e) {
            message(paste("Error in gene", x, ":", e$message))
            return(NULL)
        })
    })
    res=do.call(rbind,res)
    res=res[order(res$pvalue),]
    res["FDR"]=p.adjust(res$pvalue,"fdr")
    return(res)
}


GetPower=function(out,CT)
{
    power=map_dbl(out,function(x){
        mean(map_dbl(x[[CT]],function(y) sum(y$padj<.05 & !is.na(y$padj)) ))
    })
    return(power)
}

GetPower_CT=function(out,tab)
{
    CT=names(out[[1]])
    power=map(CT,function(x){
        tab1=tab
        tab1["CT"]=x
        tab1["Power"]=GetPower(out,x)
        tab1
    })
    dat=do.call(rbind,power)
    return(dat)
}



TestCluster<-function(dat,meta,form=Met~Clust+Age+(1|Brain),pertCol="ClusterRaw")
{
    meta[pertCol]=factor(meta[,pertCol])
    clusts=levels(meta[,pertCol])
   
    out=map(clusts,function(clust){
        meta["Clust"]=meta[,pertCol]==clust
        meta$Clust=factor(meta$Clust,levels=c("FALSE","TRUE"))
        mets=colnames(dat)
        out=map(mets,function(met){
                tryCatch({print(paste0("Testing ",met))
                tab=meta
                tab["Met"]=dat[,met]
                fit=lmer(form,data=tab)
                coef=data.frame(summary(fit)$coefficients)
                coef["Test"]=rownames(coef)
                coef["Cluster"]=clust
                coef["Met"]=met
                coef=coef[grep("Clust",rownames(coef)),]
                rownames(coef)=NULL
                return(coef)
            },error=function(e){NULL})
        })
    
        mrk=do.call(rbind,out)
        mrk=mrk[order(mrk[,5]),]
        mrk["padj"]=p.adjust(mrk[,5],"fdr")
        pert=unique(mrk$Test)
        #for(i in pert)
        #{
        #    mrk[mrk$Test==i,"padj_byTest"]=p.adjust(mrk[mrk$Test==i,5],"fdr")
        #}
        return(mrk)
    })

    return(out)

}
