library(purrr)

##
##Takes input from CellRanger and relabels, adding new column Assign
##
RelabelPerts<-function(tab)
{
    ##Gets number of guides
    numGuide=map_chr(tab$feature_call,function(x){
        s=strsplit(x,"|",fixed=T)[[1]];
        return(length(s))
    })

    ##combines paired guides for those with 2 guides
    Pert2guides=map_chr(tab$feature_call,function(x){
        s=strsplit(x,"|",fixed=T)[[1]];
        if(length(s)!=2){
            return("Unsure")
        };
        gene1=strsplit(s[1],"_",fixed=T)[[1]][1];
        gene2=strsplit(s[2],"_",fixed=T)[[1]][1];
        num1=as.numeric(sub("g","",strsplit(s[1],"_",fixed=T)[[1]][2]));
        num2=as.numeric(sub("g","",strsplit(s[2],"_",fixed=T)[[1]][2]));
        num1=floor((num1-1)/2);
        num2=floor((num2-1)/2);
        if(gene1==gene2 && num1==num2){
            return(paste(gene1,num1,sep="_"))
        }else{
            return("Doublet")
        }
    })

    ##combines paired guides for those with 1 guides
    Pert1guides=map_chr(tab$feature_call,function(x){
        s=strsplit(x,"|",fixed=T)[[1]];
        if(length(s)!=1){
            return("Unsure")
        };
        gene1=strsplit(s[1],"_",fixed=T)[[1]][1];
        num1=as.numeric(sub("g","",strsplit(s[1],"_",fixed=T)[[1]][2]));
        num1=floor((num1-1)/2);
        return(paste(gene1,num1,sep="_"))
    })

    dat=tab
    ##Add to dataframe
    dat["Assign"]="Doublet"
    dat[numGuide==1,"Assign"]=Pert1guides[numGuide==1]
    dat[numGuide==2,"Assign"]=Pert2guides[numGuide==2]

    ##Return
    return(dat)
}


##
##Test to see expected number combined vs actual number
##
TestRelabelPerts<-function(tab)
{
    tab=RelabelPerts(tab)
    assign=tab$Assign
    numGuide=map_chr(tab$feature_call,function(x){
        s=strsplit(x,"|",fixed=T)[[1]];
        return(length(s))
    })

    tab=tab[numGuide==2,]
    print(head(tab))
    guide1=map_chr(tab$feature_call,function(x){
        s=strsplit(x,"|",fixed=T)[[1]][1];
        return(s)
    })
    guide2=map_chr(tab$feature_call,function(x){
        s=strsplit(x,"|",fixed=T)[[1]][1];
        return(s)
    })

    numAgree=sum(tab$Assign!="Doublet")

    guides=table(c(guide1,guide2))
    print(guides)
    tab=data.frame(Guide=names(guides),Num=as.numeric(guides))
    tab["AfterCorrect"]=map_chr(names(guides),function(x){
        gene1=strsplit(x,"_",fixed=T)[[1]][1];
        num1=as.numeric(sub("g","",strsplit(x,"_",fixed=T)[[1]][2]));
        num1=floor((num1-1)/2);
        return(paste(gene1,num1,sep="_"))
    })
    print(head(tab))

    Prob=sum(map_dbl(1:dim(tab)[1],function(x){
        guide=tab[x,"Guide"]
        after=tab[x,"AfterCorrect"]
        prob1=sum(tab[tab$Guide==guide,"Num"])/sum(tab[,"Num"])
        prob2=sum(tab[tab$Guide!=guide & tab$AfterCorrect==after,"Num"])/sum(tab[tab$Guide!=guide,"Num"])
        return(prob1*prob2)
    }))
    numExpected=Prob*sum(tab[,"Num"])
    return(c(numAgree,numExpected,numAgree/numExpected))


}