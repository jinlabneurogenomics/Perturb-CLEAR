library(glmGamPoi)
library(Matrix)
library(HDF5Array)
library(purrr)
#library(Seurat)
library(tictoc)


runGlmGamPoi<-function(dat=NULL,meta=NULL,seur=NULL,form=~Assign+BioSamp,filter=.01,cpc=.01,filenamePref="/broad/hptmp/ssimmons/tmp/tmp",ref="NT_0",cond="Assign",compute_lfc_se=T)
{
    if(is.null(dat) && is.null(meta) && is.null(seur))
    {
        stop("Please provide either dat and meta or seur")
    }
    if(!is.null(seur))
    {
        dat=seur@assays$RNA@counts
        meta=seur@meta.data
    }
    print("Make reference")
    meta[cond]=factor(meta[,cond],levels=unique(meta[,cond]))
    meta[cond]=relevel(meta[,cond],ref)

    print("make model matrix")

    mod_mat=model.matrix(form,meta)
    #print(head(mod_mat))

    ##Filter dat
    mn=rowMeans(dat>0)
    dat=dat[mn>filter,]
    genes=rownames(dat)

    ##Save dat as an hdf5 file using HDF5Array
    #if(!file.exists(filename))
    #{
    print("Saving data to HDF5 file")
    rand=0;
    filename=paste0(filenamePref,rand,".h5")
    while(file.exists(filename))
    {
        rand=rand+1
        filename=paste0(filenamePref,rand,".h5")
    }
    print(filename)
    writeHDF5Array(dat, filepath = filename, name = "temp",with.dimnames=F)
    #}
    #else
    #{
    #    print("HDF5 file already exists, skipping save")
    #}
    ##Make into HDF5Array object
    dat_h5 <- HDF5Array::HDF5Array(filename,"temp")
    ##Test with glmgampoi
    print("Running glmGamPoi")
    tic()
    fit=glmGamPoi::glm_gp(dat_h5, design = mod_mat, verbose=TRUE,on_disk=T,subsample=T)
    toc()

    print("Fitting done, estimate DE next")
    res=map(grep(cond,colnames(mod_mat),value=T),function(coef){
        print(coef);tic();mrk=test_de(fit, coef, pval_adjust_method="BH",compute_lfc_se=compute_lfc_se);toc();mrk["Gene"]=genes;mrk=mrk[order(mrk$pval),]
    }) #,sort_by="pval",compute_lfc_se=T
    names(res)=colnames(mod_mat)[grep(cond,colnames(mod_mat))]
    map(names(res),function(x){print(x);print(head(res[[x]]))})
    return(res)

}


runGlmGamPoi_CT<-function(dat=NULL,meta=NULL,seur=NULL,form=~Assign+BioSamp,filter=.01,cpc=.01,filenamePref="/broad/hptmp/ssimmons/tmp/tmp",ref="NT_0",CT="CT",cond="Assign",minNumCT=500,minNumAssign=10,collapsePerCT=F,compute_lfc_se=T)
{
    if(is.null(meta)){meta=seur@meta.data}
    if(cond!="Assign")
    {
    meta["Assign"]=meta[,cond]
    }
    types=table(meta[,CT])
    types=names(types)[types>minNumCT]
    out=map(1:length(types),function(x){
        print(x)
        print(types[x])
        #seur2=subset(seur,CT==types[x])
        dat2=dat[,meta[,CT]==types[x]]
	meta2=meta[meta[,CT]==types[x],]
	assign=table(meta2[,"Assign"])
        assign=names(assign)[assign>minNumAssign]
        if(!(ref %in% assign))
        {
            print(paste("Reference",ref,"not in cell type",types[x]))
            return(NULL)
        }
        
        #seur2=subset(seur2,Assign %in% assign)
        dat2=dat2[,meta2$Assign %in% assign]
	meta2=meta2[meta2$Assign %in% assign,]
	filenamePref=paste0(filenamePref,"_",x,"_")
        res=runGlmGamPoi(dat=dat2,meta=meta2,form=form,filter=filter,filenamePref=filenamePref,ref=ref,cond=cond,compute_lfc_se=compute_lfc_se)
        if(collapsePerCT)
        {
            res=map(names(res),function(y){
                res2=res[[y]]
                res2["CT"]=types[x]
                res2["Assign"]=y
                return(res2)
            })
            res=do.call(rbind,res)
            res=res[order(res$pval),]
        }
	print(head(res))
        return(res)
    })
    names(out)=types
    out=out[!sapply(out,is.null)]
    return(out)
}
