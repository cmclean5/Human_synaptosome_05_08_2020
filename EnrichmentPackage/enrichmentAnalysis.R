source('../setUp.R')

#require("eva");      #GPD fit
#require("ismev");    #GPD fit
#require("fExtremes");#GPD fit
#require("MASS");     #Poission fit

# Gene Set Analysis (GSA) object
setClass(Class="STAT",representation(
                          SUM="matrix",
                          SUM2="matrix",
                          SUM3="matrix",
                          SUM4="matrix",
                          CAN="data.frame"
                      )
         )			   

fisher <- function( mu, P, F, N ){

    mm  <- matrix(c(mu,(P-mu),(F-mu),(N-P-F-mu)),2,2)
    res <- fisher.test(mm)
    return(res$p.value)
}

plotRatio <- function(xx, plotDIR, subDIR, desc="", anno="", LEGtextSize=1.5, LEGlineSize=4){
    
    #---For p.values
    #xx  = statsR1@SUM3 #Fe
    Nxx = length(colnames(xx))
    df  = data.frame()

    #--- labels
    xlab = colnames(xx)[3:Nxx]
    xval = seq(0,(length(xlab)-1),1)
    xlim = c(0,25,50,75,100)

    indx = match(xlim,xval)
    xval = xval[indx]
    xlab = xlab[indx]
    #---
    
    #--- test intervals
    X1 = 7
    X2 = 54
    X3 = 90

    rank <- matrix("",ncol=3,nrow=length(xx[,1]))

    for( i in 1:length(xx[,1]) ){

        size  = rep(1.8,length(xx[i,1]))
        alpha = rep(0.7,length(xx[i,1]))
        col   = rep("grey", length(xx[i,1]))
        xlabs = colnames(xx)[3:Nxx]
        
        tmp = cbind(xx[i,1],seq(1,(Nxx-2),1),as.numeric(xx[i,3:Nxx])/as.numeric(xx[i,2]), size, alpha, col, xlabs)  

        df = rbind(df,tmp)

        zz0 = as.numeric(as.vector(xx[i,2]))
        zz  = as.numeric(as.vector(xx[i,3:Nxx]))
    
        rank[i,1] = xx[i,1]
        rank[i,2] = ifelse( zz0 == 0, 0, (zz[7]  - zz[54])/zz0)
        rank[i,3] = ifelse( zz0 == 0, 0, (zz[54] - zz[90])/zz0)
    
    }

    #---
    rank = rank[order(as.numeric(rank[,2]),decreasing=T),]
    colnames(rank) <- c("Alg","FacComsEnriched_log2(FE)>0.5_log2(FE)<4.8","FacComsEnriched_log2(FE)>4.8_log2(FE)<8.0")    
    write.table(rank, sprintf("ranking_%s.csv",desc), sep="\t", row.names=F, col.names=T, quote=F)

    
    colnames(df) <- c("ALG","X","Y","LSIZE","ALPHA","COL","XLAB")
    df           <- df[order(match(df[,1],rank[,1])),]
    df$ALG       <- factor(df$ALG, levels=rank[,1])
    df           <- as.data.frame(df)        

    if( length(which(df$ALG == "Spectral")) != 0 ){
        df$COL[df$ALG   == "Spectral"]="royalblue"#lawngreen"
        df$ALPHA[df$ALG == "Spectral"]= 1.0
    
     }

    if( length(which(df$ALG == "SVI")) != 0 ){
        df$COL[df$ALG   == "SVI"]="green2"
        df$ALPHA[df$ALG == "SVI"]= 1.0
    
    }

    
    if( length(which(df$ALG == "louvain2")) != 0 ){
        df$COL[df$ALG   == "louvain2"]="magenta"    
        df$ALPHA[df$ALG == "louvain2"]= 1.0
    }
    #---
    
    #---Generate plot
    SIZEa=2
    SIZEb=2

    #---legend
    #LEGtextSize=0.75 #ALL Algs
    #LEGtextSize=1.5   #Selected ALgs

    #LEGlineSize=2    #ALL Algs
    #LEGlineSize=4    #Selected ALgs
    
    colours = df$COL[match(levels(factor(df$ALG)),df$ALG)]
    
    gplot <- ggplot(df,aes(x=(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
    geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
    labs(x="log2(Fe)",y="Fraction of Enriched Communities",title=anno)+
    theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
          axis.title.y=element_text(face="bold",size=rel(1.5)),
          legend.text=element_text(face="bold",size=rel(LEGtextSize)),
          plot.title=element_text(face="bold",size=rel(1.5)),
          legend.position="bottom")+
        scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = FALSE,
               size  = FALSE)

png(sprintf("%s/%sFe_%s.png",plotDIR,subDIR,desc),width=WIDTH,height=HEIGHT,units="px");
print(gplot)
dev.off()

  
    
    gplot2 <- ggplot(df,aes(x=log(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
        geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
        labs(x="log(log2(Fe))",y="Fraction of Enriched Communities",title=anno)+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(LEGtextSize)),
              plot.title=element_text(face="bold",size=rel(1.5)),
               legend.position="bottom")+
        scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=log(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = FALSE,
               size  = FALSE)

png(sprintf("%s/%slogFe_%s.png",plotDIR,subDIR,desc),width=WIDTH,height=HEIGHT,units="px");
print(gplot2)
dev.off()


library(viridis)
    gplot3 <- ggplot(df,aes(x=(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
    geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
    labs(x="log2(Fe)",y="Fraction of Enriched Communities",title=anno)+
    theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
          axis.title.y=element_text(face="bold",size=rel(1.5)),
          legend.text=element_text(face="bold",size=rel(LEGtextSize)),
          plot.title=element_text(face="bold",size=rel(1.5)),
          legend.position="bottom")+
        scale_color_viridis("",discrete = TRUE, option = "D")+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = FALSE,
               size  = FALSE)

png(sprintf("%s/%sFe_col_%s.png",plotDIR,subDIR,desc),width=WIDTH,height=HEIGHT,units="px");
print(gplot3)
dev.off()

 gplot4 <- ggplot(df,aes(x=log(as.numeric(df$X)),y=as.numeric(as.vector(df$Y)),colour=df$ALG))+
        geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
        labs(x="log(log2(Fe))",y="Fraction of Enriched Communities",title=anno)+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(LEGtextSize)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        scale_color_viridis("",discrete = TRUE, option = "D")+
        scale_y_continuous(expand=c(0,0),limits=c(0,1))+
        scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_vline(xintercept=log(as.numeric(X1)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X2)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        geom_vline(xintercept=log(as.numeric(X3)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
        guides(color = guide_legend(override.aes = list(size=LEGlineSize)),
               alpha = FALSE,
               size  = FALSE)

png(sprintf("%s/%slogFe_col_%s.png",plotDIR,subDIR,desc),width=WIDTH,height=HEIGHT,units="px");
print(gplot4)
dev.off()
    
}

addDescription <- function(AFILE, CAND, ID, IDr){

    #ID="GO:", IDr="GO."
    indx = gsub(ID,IDr,AFILE[,1])

    desc = AFILE[match(CAND[,2],indx),2]

    CN <- colnames(CAND)
    oo <- cbind(CAND[,1],CAND[,2],desc,CAND[,3],CAND[,4])
    colnames(oo) <- c(CN[1],CN[2],"DES",CN[3],CN[4])

    return(oo)
}

summaryStats <- function( RES, ALPHA, usePadj=FALSE, FeMAX=0, FcMAX=0 ){

    cand <- data.frame(a=as.character(),b=as.character(),c=as.character(),d=as.character())
    
    CN <- colnames(RES[[1]])

    ALGi   = which(CN=="Alg")[1]
    Pvi    = which(CN=="Pv")[1]
    PvALTi = which(CN=="PvALT")[1]
    if( usePadj ){
        Pvi    = which(CN=="Ap")[1]
        PvALTi = which(CN=="ApALT")[1]
    }
    ORi    = which(CN=="OR")[1]
    CIli   = which(CN=="CIl")[1]
    Fei    = which(CN=="Fe")[1]
    Fci    = which(CN=="Fc")[1]

    Ci     = which(CN=="C")[1]
    Cni    = which(CN=="Cn")[1]
    Mui    = which(CN=="Mu")[1]

    Ni     = which(CN=="N")[1]
    N      = as.numeric(RES[[1]][1,Ni])

    hh0  = c("alg","FN","CN","FNxCN","Psig","PALTsig","OR>1","ORsig","Psig&ORsig","PALTsig&ORsig","FEsig","Psig&ORsig&FEsig","EnrichedComs(%)","p.value")
    sum1 = matrix("",ncol=length(hh0),nrow=length(names(RES)))
    colnames(sum1) <- hh0

    cmin = seq(0,10,0.1)
    hh   = sprintf("Cmin_%.1f&Cmax_%.1f",cmin,10)
    sum2 = matrix("",nrow=length(names(RES)), ncol=(2+length(hh)) )
    colnames(sum2) = c("Alg","Psig&ORsig",hh)

    steps = 100
    hh2i  = round(FeMAX/steps,3)
    femin = seq(0,ceiling(FeMAX), hh2i) 
    hh2   = sprintf("%.1f",femin)
    sum3  = matrix("",nrow=length(names(RES)), ncol=(2+length(hh2)) )
    colnames(sum3) = c("Alg","Psig&ORsig",hh2)

    hh3i  = round(FcMAX/steps,3)
    fcmin = seq(0,ceiling(FcMAX), hh3i) 
    hh3   = sprintf("%.1f",fcmin)
    sum4  = matrix("",nrow=length(names(RES)), ncol=(2+length(hh3)) )
    colnames(sum4) = c("Alg","Psig&ORsig",hh3)
    
    for( i in 1:length(RES) ){

    Ncn = length(RES[[i]][,1])
    P   = as.numeric(RES[[i]][,Pvi])
    Palt= as.numeric(RES[[i]][,PvALTi])
    OR  = as.numeric(RES[[i]][,ORi])          
    CI  = as.numeric(RES[[i]][,CIli])          
    FE  = as.numeric(RES[[i]][,Fei])
    FC  = as.numeric(RES[[i]][,Fci])

    CNo = as.numeric(RES[[i]][,Cni])
        
    Cmax = max(as.numeric(RES[[i]][,Ci]))    
    
    sum1[i,1] = as.character(RES[[i]][1,ALGi])

    sum1[i,2] = Ncn / Cmax

    sum1[i,3] = Cmax
        
    sum1[i,4] = Ncn ##

    sum1[i,5] = sum(P <= ALPHA)

    sum1[i,6] = sum(Palt <= ALPHA)
        
    sum1[i,7] =  sum(OR > 1)
    
    sum1[i,8] =  sum(OR > 1 & CI > 1)

    sum1[i,9] =  sum(OR > 1 & CI > 1 & P <= ALPHA)

    sum1[i,10] =  sum(OR > 1 & CI > 1 & Palt <= ALPHA)
        
    sum1[i,11] = sum( log2(FE) > 0.5 & log2(FE) < 4.8 )
        
    sum1[i,12] = sum(OR > 1 & CI > 1 & P <= ALPHA & log2(FE) > 0.5 & log2(FE) < 4.8 )

   sum2[i,1] =  as.character(RES[[i]][1,ALGi])
   sum2[i,2] =  sum(OR > 1 & CI > 1 & P <= ALPHA)
    for( j in 1:length(hh) ){
        sum2[i,(j+2)] = sum(OR > 1 & CI > 1 & P <= ALPHA & CNo > ((cmin[j]*N)/100) & CNo < ((10*N)/100) )
    }

   sum3[i,1] =  as.character(RES[[i]][1,ALGi])
   sum3[i,2] =  sum(OR > 1 & CI > 1 & P <= ALPHA)
    for( j in 1:length(hh2) ){
        sum3[i,(j+2)] = sum(OR > 1 & CI > 1 & P <= ALPHA & log2(FE) > femin[j] )
    }
      
   sum4[i,1] =  as.character(RES[[i]][1,ALGi])
   sum4[i,2] =  sum(OR > 1 & CI > 1 & P <= ALPHA)
    for( j in 1:length(hh3) ){
        sum4[i,(j+2)] = sum(OR > 1 & CI > 1 & P <= ALPHA & log2(FC) > fcmin[j] )
    }
      
    sum1[i,13] = (as.numeric(sum1[i,9]) / as.numeric(sum1[i,4])) * 100

    sum1[i,14] = fisher(as.numeric(sum1[i,12]), as.numeric(sum1[i,9]),
                        as.numeric(sum1[i,11]), as.numeric(sum1[i,4]))
        
   #store candidate enriched functional communities
   indx <- OR > 1 & CI > 1 & P <= ALPHA & log2(FE) > 0
   n    <- sum(indx)     
   A    <- rep(sum1[i,1],n)
   B    <- RES[[i]][indx,1]
   C    <- RES[[i]][indx,2]
   D    <- RES[[i]][indx,Mui]
   cand <- rbind(cand,data.frame(A,B,C,D))
    }

  colnames(cand) <- c("ALG","Fn","C","Mu")
        
    return(new("STAT",
               SUM=sum1,
               SUM2=sum2,
               SUM3=sum3,
               SUM4=sum4,
               CAN=cand))
}

addColour <- function( CVEC, ALGN, COL ){

    indx = which(CVEC==ALGN)
    if( length(indx) != 0 ){
        CVEC[indx[1]] = COL
    }

    return(CVEC)
}

Scale <- function(X){

    X = as.numeric(as.vector(X))
    X = (X-min(X))/(max(X)-min(X))
    return(X)
    
}

getCI <- function(x, indx=1){
    strsplit(x,",")[[1]][indx]
}

#---Get all edges internal to a community
intraEdges <- function(GG, ALG, CC, INTRA=NULL, INTER=NULL){

    intra = NULL #edges in the community CC
    inter = NULL #edges going out from community CC   
    
    if( !is.null(igraph::get.vertex.attribute(GG,ALG)) ){

        coms <- get.vertex.attribute(GG,ALG)

        if( length(which(coms == CC)) != 0 ){
        
            ed_cc = E(GG)[inc(coms == CC)]

            all_edges_m <- get.edges(GG, ed_cc) #matrix representation

            inter = (ed_cc[!(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])

            intra = (ed_cc[(all_edges_m[, 1] %in% V(GG)[coms == CC] & all_edges_m[, 2] %in% V(GG)[coms == CC])])
            
        }
    }

    if( INTRA==TRUE && !is.null(intra) ){
        intra_m = get.edges(GG,intra)
        intra   = cbind(V(GG)$name[intra_m[,1]],V(GG)$name[intra_m[,2]])
        return(intra)
    }

    if( INTER==TRUE && !is.null(inter) ){
        inter_m = get.edges(GG,inter)
        inter   = cbind(V(GG)$name[inter_m[,1]],V(GG)$name[inter_m[,2]])
        return(inter)
    }

    return(NULL)
        
}



#---OUT Dir
OUT    <- vector(length=4)
OUT[1] <- DIRS[grepl("EnrichmentPackage",DIRS)]
OUT[2] <- DIRS[grepl("Graphs",DIRS)]
OUT[3] <- DIRS[grepl("parameterFiles",DIRS)]
OUT[4] <- DIRS[grepl("Annotations",DIRS)]

#---Check or create output dir
enrdir <- sprintf("%s/%s/%s",OUT[1],"RESULTS",subDIR[S])
if( !file_test("-d",enrdir) ){
    dir.create(enrdir)
}

grdir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}

plotDIR <- sprintf("%s/PLOTS",OUT[1])
if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}
#---

#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
ids <- V(gg)$name;
N   <- length(V(gg))

#---READ IN ALGORITHM 
algs <- read.delim(sprintf("%s/clusteringAlg.csv",OUT[3]),header=F,sep="\t",quote="")
ALGS <- as.vector(algs[which(as.vector(algs[,1]) == 1),3])

#---READ IN ANNOTATION
anno  <- read.delim(sprintf("%s/annotations.csv",OUT[3]),header=F,sep="\t",quote="")
ANNO  <- as.vector(anno[which(as.vector(anno[,1]) == 1),3])
Afile <- as.vector(anno[which(as.vector(anno[,1]) == 1),2])
afile <- read.delim(sprintf("%s/%s",OUT[4],Afile),sep="\t", skip=1,header=F)

RES <- list()

Anno=1
k=1

FcMAX=0
FeMAX=0

searchTERM="GO"
if(grepl("GO",ANNO)){
    searchTERM="GO"
}

if(grepl("topOnto_ovg",ANNO)){
    searchTERM="DOID"
}

for( a in 1:length(ALGS) ){

    str <- sprintf("%s/%s/permute_p_values_%s.csv",enrdir,ALGS[a],ANNO[Anno])

     if( file.exists(str) && file.info(str)$size!=0 ){
    
         #--- load functional enrichment file
         tt <- read.delim(str,sep="\t", header=T)

         indx <- grep(searchTERM,colnames(tt))         
         fn   <- colnames(tt)[indx]
         FN   <- length(fn)

         CN   <- length(tt[,1])

         DF <- matrix("",nrow=(FN*CN), ncol=17)
         colnames(DF) <- c("Fn","C","Cn","N","Mu","OR","CIl","CIu","Pv","Ap","PvALT","ApALT","Fe","E(Mu)","F","Alg","Fc")

         DF[,1]  = rep(fn,CN)
         DF[,4]  = rep(N,(FN*CN))
         DF[,15] = rep(as.vector(unlist(tt[1,indx])),CN)
         DF[,16] = rep(ALGS[a],(FN*CN))

         temp1  <- c()
         temp2  <- c()
         temp3  <- c()
         temp4  <- c()
         temp5  <- c()
         temp6  <- c()
         temp7  <- c()
         temp8  <- c()
         temp9  <- c()
         temp10 <- c()
         temp11 <- c()
         temp12 <- c()
         
         for( i in 1:length(tt[,1]) ){

             cc    <- rep(tt[i,1],FN)
             temp1 <- c(temp1, cc)

             cn    <- rep(tt[i,2],FN)
             temp2 <- c(temp2, cn)

             ov = as.vector(unlist(tt[i,grepl("actual",colnames(tt))]))
             temp3 <- c(temp3, ov)

             ep = as.vector(unlist(tt[i,grepl("expected",colnames(tt))]))
             temp12 <- c(temp12, ep)
             
             or = as.vector(unlist(tt[i,grepl("OR",colnames(tt))]))
             temp4 <- c(temp4, or)

             ci = as.vector(unlist(tt[i,grepl("CI",colnames(tt))]))
             ci = gsub("\\[","",ci)
             ci = gsub("\\]","",ci)
             temp5 <- c(temp5,as.vector(sapply(ci, getCI, indx=1)))
             temp6 <- c(temp6,as.vector(sapply(ci, getCI, indx=2)))
                          
             pv = as.vector(unlist(tt[i,grepl("p.value", colnames(tt)) &
                                        !grepl("X.p_value",colnames(tt)) &
                                        !grepl("ALT",colnames(tt))]))
             temp7 <- c(temp7, pv)
             
             pa = as.vector(unlist(tt[i,grepl("adjusted",colnames(tt)) &
                                      !grepl("ALT",colnames(tt))]))
             temp8 <- c(temp8, pa)

             palt = as.vector(unlist(tt[i,grepl("p.value.ALT.",colnames(tt))]))
             temp9 <- c(temp9, palt)
             
             alt = as.vector(unlist(tt[i,grepl("adjusted.ALT.",colnames(tt))]))
             temp10 <- c(temp10, alt)

             fc = (as.numeric(ov)/as.numeric(cn)) / (as.numeric(cn)/as.numeric(N))
             temp11 <- c(temp11, fc)

             if( max(fc) > FcMAX ){
                 FcMAX = max(fc)
             }
             
             #rm(cc,cn,ov,ep,or,ci,pv,pa,palt,alt,fc)
             
         }
    
         DF[,2]  = temp1
         DF[,3]  = temp2
         DF[,5]  = temp3
         DF[,6]  = temp4
         DF[,7]  = temp5
         DF[,8]  = temp6
         DF[,9]  = temp7
         DF[,10] = temp8
         DF[,11] = temp9
         DF[,12] = temp10
         DF[,13] = temp11
         DF[,14] = temp12

         DF[,17] = (as.numeric(DF[,5])/as.numeric(DF[,15])) / (as.numeric(DF[,3])/as.numeric(DF[,4]) )
         if( max(as.numeric(DF[,17])) > FeMAX ){
             FeMAX = max(as.numeric(DF[,17]))
         }
         
         rm(tt)
         rm(cc,cn,ov,or,ci,pv,pa,palt,alt,fc,ep)
         rm(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12)

         DF = as.data.frame(DF)

         RES[[k]] = DF
         names(RES)[k] = sprintf("%s_%s",ALGS[a],ANNO[Anno])
         k=k+1
         rm(DF)
     }

}

FeMAX=log2(FeMAX)
FcMAX=log2(FcMAX)

#p.value
statsR1 <- summaryStats( RES, alpha[1], FeMAX=FeMAX, FcMAX=FcMAX )
stats1  <- statsR1@SUM
ranki   <- which(colnames(stats1)=="EnrichedComs(%)")
stats1  <- stats1[order(as.numeric(stats1[,ranki]),decreasing=T),]

write.table(stats1,      "pvalue_stats.csv", sep="\t", row.names=F, col.names=T, quote=F)

cands <- addDescription(afile, statsR1@CAN, "GO:", "GO.")
write.table(cands, "pvalues_candidate_C.csv", sep="\t", row.names=F, col.names=T, quote=F)

#adjusted p.value
statsR2 <- summaryStats(RES, alpha[1], usePadj=TRUE, FeMAX=FeMAX, FcMAX=FcMAX)
stats2  <- statsR2@SUM
ranki   <- which(colnames(stats2)=="EnrichedComs(%)")
stats2  <- stats2[order(as.numeric(stats2[,ranki]),decreasing=T),]

write.table(stats2,      "adjusted_stats.csv", sep="\t", row.names=F, col.names=T, quote=F)

cands <- addDescription(afile, statsR2@CAN, "GO:", "GO.")
write.table(cands, "adjusted_candidate_C.csv", sep="\t", row.names=F, col.names=T, quote=F)

#safe
write.table(statsR1@SUM3,"test.data.csv", sep="\t", col.names=T, row.names=F, quote=F)


#run=0
run=1

#plot Ratio of Enriched Communities Verses Fold-Enrichment
if( run ){

    #---For p.values
    plotRatio(xx=statsR1@SUM3, plotDIR=plotDIR, subDIR=subDIR[S],
              desc="p.values", anno=ANNO,
              LEGtextSize=0.75, LEGlineSize=2)

    #---For adjusted p.values
    plotRatio(xx=statsR2@SUM3, plotDIR=plotDIR, subDIR=subDIR[S],
              desc="adjusted", anno=ANNO,
              LEGtextSize=0.75, LEGlineSize=2)
    
}
