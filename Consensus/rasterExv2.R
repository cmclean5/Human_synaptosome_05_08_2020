source('../setUp.R')

require(raster)
require(spatstat)
require(colorspace)
require(rasterVis)
library(gridExtra)
require(scatterplot3d)

require(plot3D)
require(rgl)

saveDatainCell <- function(DF=NULL, GG=NULL, GRID=NULL, NBIN=NULL, CELLNo=NULL, ALGname=NULL, ORGX=NULL){

    if( !is.null(GRID) ){

        xy  <- cbind(DF[,1],DF[,2])
        r1  <- rasterize(cbind(DF[,1],DF[,2]), GRID, DF[,3])
        res <- cbind(raster::extract(r1,xy, df=T, cellnumbers=T),DF)

        if( !is.null(ORGX) ){ res <- cbind(res,ORGX) }
        
        indx <- match(res[,6],V(GG)$name)
        if( length(indx) != 0 ){

            gn   = get.vertex.attribute(GG,"GeneName",  indx)
            alg  = get.vertex.attribute(GG,ALGname,     indx)
            dis  = get.vertex.attribute(GG,"TopOntoOVG",indx)
            gobp = get.vertex.attribute(GG,"GOBP",      indx)
            gocc = get.vertex.attribute(GG,"GOCC",      indx)
            gomf = get.vertex.attribute(GG,"GOMF",      indx)

            res2 <- cbind(res,gn,alg,dis,gobp,gocc,gomf)

            if( !is.null(CELLNo) ){
                res3 <- res2[res2[,2]==CELLNo,]
                colnames(res3) <- c("ID","cell","layer","x","y","EntrezID","ORGX","GeneName",ALGname,"topOntoOVG","GOBP","GOCC","GOMF")
                res3 <- as.data.frame(res3)

                if( !is.null(ORGX) ){
                    res3   = res3[c(1,2,3,7,5,6,8,9,10,11,12,13,4)]
                    res3   = res3[-c(13)]
                    colnames(res3) <- c("ID","cell","layer","x","y","EntrezID","GeneName",ALGname,"topOntoOVG","GOBP","GOCC","GOMF")
                }
                
                write.table(res3,sprintf("dataInGrid_%dx%d_CellNo_%d.csv",NBIN,NBIN,CELLNo),sep="\t",row.names=F,col.names=T,quote=F)
            }
        }
        
    }
}

scale <- function(X, VALUE=NULL){

    X = as.numeric(as.vector(X))

    xmin <- min(X,na.rm=T)
    xmax <- max(X,na.rm=T)
   
    if( is.null(VALUE) ){

        X  <- X-xmin
        X  <- ifelse(!is.na(X), X/(xmax-xmin), NA) 
        #X    <- (X-min(X,na.rm=T))/(max(X,na.rm=T)-min(X,na.rm=T))
        return(X)
    }

    VALUE = as.numeric(as.vector(VALUE)[1])
    VALUE = VALUE-xmin
    VALUE = ifelse(!is.na(VALUE), VALUE/(xmax-xmin), NA) 
    return(VALUE)
}

transformQ <- function( cell, q ){

    N   = length(cell[,1])
    arr = rep(0,N)
    nrow=length(q[,1])
    ncol=length(q[1,])
    
    for( k in 0:(N-1) ){
        i = floor(k/ncol)
        j = k %% ncol
        arr[(k+1)] = q[(i+1),(j+1)]    
    }

    return(arr)
}

getCI <- function(x, indx=1){
    strsplit(x,",")[[1]][indx]
}


addNoise <- function( Y, MN=0, SD=0.05 ){
    return( Y+rnorm(length(Y),mean=MN, sd=SD) )
}

nearestSqrt <- function ( n ){

    ss = sqrt(n)
    ff = floor(ss)
    cc = ceiling(ss)

    ff2 = ff*ff
    cc2 = cc*cc

    fdis = n - ff2
    cdis = cc2 - n

    ret=NA
    
    if( fdis < cdis ){
        ret=ff2
    } else {
        ret=cc2
    }

    return(ret)
    
}

findcellNo <- function(X,Y, xcell, ycell, Xmax=1, Ymax=1, COUNT=TRUE, ROUND=4){

    xcell = round(as.vector(xcell),ROUND)
    ycell = round(as.vector(ycell),ROUND)
    
    names(xcell) = seq(1,length(xcell),1)
    names(ycell) = seq(1,length(ycell),1)

    Xres = seq(0,Xmax,length.out=(1+length(unique(xcell))))
    Yres = seq(0,Ymax,length.out=(1+length(unique(ycell))))

    X = round(as.numeric(X),ROUND)
    Y = round(as.numeric(Y),ROUND)
    
    N = length(X)
    cell = rep(0,N)

    if( COUNT ){
    
    for( i in 1:N ){
        cell[i] = names(xcell)[round(X[i],ROUND) == xcell & Y[i] == ycell]
    }

    }
        
    bound = matrix(0,ncol=4,nrow=N)
    
    for( i in 2:length(Xres) ){
        for( j in 2:length(Yres) ){
            indx = X > Xres[i-1] & X <= Xres[i] & Y > Yres[j-1] & Y <= Yres[j]
            bound[indx,1] = Xres[i-1]
            bound[indx,2] = Xres[i]
            bound[indx,3] = Yres[j-1]
            bound[indx,4] = Yres[j]
        }
    }
    
    return(cbind(as.numeric(X),
                 as.numeric(Y),
                 as.numeric(cell),
                 bound))
    
}

readEFile <- function( str="", TERM="GO", algName, ANNO, N=1 ){
    
    RES <- list()
    k=1    
    
     if( file.exists(str) && file.info(str)$size!=0 ){
    
         #--- load functional enrichment file
         tt <- read.delim(str,sep="\t", header=T)

         indx <- grep(TERM,colnames(tt))         
         fn   <- colnames(tt)[indx]
         FN   <- length(fn)

         CN   <- length(tt[,1])

         headers <- c("Fn","C","Cn","N","Mu","OR","CIl","CIu","Pv","Ap","PvALT","ApALT","Alg","E(mu)")
         DF <- matrix("",nrow=(FN*CN), ncol=length(headers))
         colnames(DF) <- headers

         DF[,1]  = rep(fn,CN)
         DF[,4]  = rep(N,(FN*CN))
         DF[,13] = rep(algName,(FN*CN))

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
         DF[,14] = temp12

         rm(tt)
         rm(cc,cn,ov,or,ci,pv,pa,palt,alt,ep)
         rm(temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10, temp11, temp12)

         DF = as.data.frame(DF)

         RES[[k]] = DF
         names(RES)[k] = sprintf("%s_%s",algName,ANNO)
         k=k+1
         rm(DF)
     }


    return(RES)
    
}

#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("EntropyRate",DIRS)]
OUT[3] <- DIRS[grepl("EnrichmentPackage",DIRS)]


##---Check or create output plot dir
plotDIR <- sprintf("%s/PLOTS/",subDIR[S])

if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}

##---Check or create output REGIONS dir
regDIR <- sprintf("%s/REGIONS/",subDIR[S])

if( !file_test("-d",regDIR) ){
    dir.create(regDIR)
}

  

#Read-in all measures    
meas <- read.table(sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]),sep="\t",header=T)

#Read-in Entropy values
ent <- read.table(sprintf("%s/%s/SignalEntropyRate.csv",OUT[2],subDIR[S]),sep="\t",header=T)

algName    = vector(length=3)
algName[1] = "Spectral"
#algName[2] = "louvain2"
#algName[3] = "SVI"

aln = 1

Y=meas[,grepl(sprintf("BRIDGE_%s",algName[aln]),colnames(meas))]


Cent <- vector(length=2)
Cent[1] = "Cl"
Cent[2] = "Entropy"

ce=1

#X value - Cl
if( Cent[c] == "Cl" ){
    X=meas[,8]
    Xl=X
}

#X value - Entropy UP reg.
if( Cent[ce] == "Entropy" ){
    X=ent[,4]
    Xl = X
    IntCon=scale(X,VALUE=0.74024)
}

X=scale(X)

xy <- cbind(X,Y)
xy <- data.frame(xy)
names(xy) <- c('x','y')

#define x-min/max and y-min/max
Xmin=0
Xmax=1

Ymin=0
Ymax=1

#---WIDTH and HEIGHT for plots
WIDTH=750
HEIGHT=750

#create grid
#nbins=c(2,3,4,5,6,7,8,9,10,15,20,25)
nbins=c(25,30,35,40)

#xy points in grid
PP = ppp(x=xy[,1],y=xy[,2], xrange=c(Xmin,Xmax), yrange=c(Ymin,Ymax))

Gr = list()

for( i in 1:length(nbins) ){

    Xbins=nbins[i]
    Ybins=nbins[i]

    rr      = raster(crs=NA, ncol=Xbins, nrow=Ybins, xmn=Xmin, xmx=Xmax, ymn=Ymin, ymx=Ymax)
    Gr[[i]] = rr

    rm(rr)
    
}

Rs    = list()
Qtr   = list()
Xcell = list()
Ycell = list()

for( i in 1:length(nbins) ){

    rs   = rasterize(cbind(xy[,1],xy[,2]),Gr[[i]], fun='count')    
    Q    = quadratcount(PP, nx=Gr[[i]]@ncols, ny=Gr[[i]]@nrows)
    CELL = coordinates(Gr[[i]])
    qtr  = transformQ( CELL, Q )

    Rs[[i]]    = rs
    Qtr[[i]]   = qtr
    Xcell[[i]] = CELL[,1]
    Ycell[[i]] = CELL[,2]

    rm(rs, qtr, CELL)
    
}


plotRaster = 1

if( plotRaster ){

    #See http://hclwizard.org/hclwizard/
    myTheme <- rasterTheme(sequential_hcl(n = 7, h = c(-80, 78), c = c(60, 75, 55), l = c(40, 91), power = c(0.8, 1), rev = TRUE), 
                           layout.heights =
                               list(top.padding       = 0,
                                    main.key.padding  = 0,
                                    key.axis.padding  = 0,
                                    axis.xlab.padding = 0,
                                    xlab.key.padding  = 0,
                                    key.sub.padding   = 0,
                                    bottom.padding    = 0),
                           clip = list(panel="off"),
                           layout.widths =
                               list(left.padding      = 0,
                                    key.ylab.padding  = 0,
                                    ylab.axis.padding = 0,
                                    axis.key.padding  = 0,
                                    right.padding     = 0))

    #sequence for axis
    my.at <- seq(0, 1, by = 0.2)

    if( Cent[ce] == "Cl" ){
        Xtck = seq(min(X), max(X), length.out=length(my.at))
        Xlab = round(seq(min(Xl),max(Xl),length.out=length(Xtck)),3)
    }
    
    if( Cent[ce] == "Entropy" ){
        Xtck = seq(min(X), max(X), length.out=3)#length(my.at))
        Xlab = round(seq(min(Xl),max(Xl),length.out=length(Xtck)),3)
    }

       
    #list of raster objects of different grids
    k=1
    p=1
    while( k <= length(nbins) ){

        m = matrix(1:4, nrow=2)
        
        png(filename=sprintf("%s/%s_%s_rasterPlots_page%d.png",plotDIR,subDIR[S],algName[aln],p), width=WIDTH, height=HEIGHT, units="px")

        for( i in 1:4 ){

            valSize=c(0.5)
            if( nbins[k] > (15) ){ valSize=c(0.25) }
        
            pt = levelplot(Rs[[k]], layers=1, aspect.ratio="iso", margin=list(FUN='sum'), par.settings = myTheme, scales=list(x=list(at=Xtck, labels=Xlab, cex=0.8), y=list(at=my.at, labels=my.at, cex=0.8)), xlab=list(label=Cent[ce], cex=0.8), ylab=list(label="Bridgeness", cex=0.8), main=list(label=sprintf("%s Cells(%dx%d)",algName[aln], nbins[k], nbins[k]), cex=0.8))+
                layer(panel.text(x=Xcell[[k]], y=Ycell[[k]], as.vector(Qtr[[k]]), font=2, cex=valSize, alpha=0.7))+
                layer(panel.abline(h=0.5, col.line="grey40", lty="longdash", cex=1, alpha=c(0.8)))+
                layer(panel.abline(v=0.5, col.line="grey40", lty="longdash", cex=1, alpha=c(0.8))) 

            if( Cent[ce] == "Entropy" ){
                pt = levelplot(Rs[[k]], layers=1, aspect.ratio="iso", margin=list(FUN='sum'), par.settings = myTheme, scales=list(x=list(at=Xtck, labels=Xlab, cex=0.8), y=list(at=my.at, labels=my.at, cex=0.8)), xlab=list(label=Cent[ce], cex=0.8), ylab=list(label="Bridgeness", cex=0.8), main=list(label=sprintf("%s Cells(%dx%d)",algName[aln], nbins[k], nbins[k]), cex=0.8))+
                    layer(panel.text(x=Xcell[[k]], y=Ycell[[k]], as.vector(Qtr[[k]]), font=2, cex=valSize, alpha=0.7))+
                    layer(panel.abline(h=0.5, col.line="grey40", lty="longdash", cex=1, alpha=c(0.8)))+
                    layer(panel.abline(v=0.5, col.line="grey40", lty="longdash", cex=1, alpha=c(0.8)))+
                    layer(panel.abline(v=IntCon, col.line="red", alpha=c(0.8), lty="longdash"))
            }

            
    #Plts[[i]] = pt

            print(pt, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<=4), axes=FALSE, legend=FALSE, asp=0.2)
    
            rm(pt)
            k=k+1
        }
        
    #x11()
        dev.off()
        p=p+1
        
    }

}


genesInCell = 1

if( genesInCell ){

#store gene Entrez IDs in each cell

# cbind(x,y,value)
df = cbind(xy[,1], xy[,2], meas[,1])
colnames(df) <- c("x","y","EntrezID")

for( i in 1:length(nbins) ){

    #rasterize(cbind(x,y), X, value)
    r1   <- rasterize(cbind(df[,1],df[,2]), Gr[[i]], df[,3])
    #r1[] <- 1:ncell(r1)

    res = cbind(raster::extract(r1,xy, df=T, cellnumbers=T),df, meas[,2])

    cc  = cbind(as.character(res[,6]),as.character(res[,2]))

    outfile <- file(sprintf("%s/%s_cells_%d.csv", regDIR, algName[aln], nbins[i]),"w")
    cat("#regions",file=outfile,"\n")
    write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
    close(outfile);

    rm(r1, res, cc)

}

}


cellData=0

if( cellData ){

# cbind(x,y,value)
df = cbind(xy[,1], xy[,2], meas[,1])
colnames(df) <- c("x","y","EntrezID")

#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

#cno=30  #louvain2, Gr[[8]], nbins[8], Cl
#cno=304 #louvain2, Gr[[12]], nbins[12], entropy
#cno=23  #SVI,      Gr[[6]], nbins[6], Cl
    
saveDatainCell( DF=df, GG=gg, GRID=Gr[[8]], NBIN=nbins[8], CELLNo=30, ALGname=algName[aln], ORGX=Xl )

    
}

BINS = nbins#[1:9]

subdir        = vector(length=2)
subdir[1]     = "topOnto_ovg"
subdir[2]     = "GOBP"

searchTERM    = vector(length=2)
searchTERM[1] = "DOID"
searchTERM[2] = "GO"

ANNO          = vector(length=2)
ANNO[1]       = "topOntoOVG"
ANNO[2]       = "GOBP"

TARGET        = vector(length=2)
TARGET[1]     = "AD"
TARGET[2]     = "GO:0000165"

N=4495

s=1

FILE=sprintf("%s/RESULTS/PPI_PSP/CELLS/%s/%s/%s",OUT[3],algName[aln],Cent[ce],subdir[s])


RES <- list()

for( b in 1:length(BINS) ){

    str = sprintf("%s/permute_p_values_%d.csv",FILE,BINS[b])
    
    df <- readEFile( str=str, TERM=searchTERM[s], algName=algName[aln], ANNO=ANNO[s], N=N )
    df <- df[[1]]


    indx  = which(nbins==BINS[b])
    cells = seq(1,length(Xcell[[indx]]),1)
    ic    = match(df[,2],cells)

    DF <- cbind(df[,1], df[,2], df[,7], df[,9], df[,10], Xcell[[indx]][ic], Ycell[[indx]][ic])
    DF <- as.data.frame(DF)

    #set-up XY grid position for plot 
    XYbin = sqrt( nearestSqrt( length(DF[,1]) ))

    XYgrid = raster(crs=NA, ncol=XYbin, nrow=XYbin, xmn=Xmin, xmx=Xmax, ymn=Ymin, ymx=Ymax)
    XYgcor = coordinates(XYgrid)

    XYobs   = findcellNo( DF[,6], DF[,7], Xcell[[indx]], Ycell[[indx]])
    XYobs2  = XYobs[!duplicated(XYobs),]

    XYnew   = findcellNo( XYgcor[,1], XYgcor[,2], Xcell[[indx]], Ycell[[indx]], COUNT=FALSE)

    for( i in 1:length(XYnew[,1]) ){
        INDX = XYnew[i,4] == XYobs2[,4] & XYnew[i,5] == XYobs2[,5] & XYnew[i,6] == XYobs2[,6] & XYnew[i,7] == XYobs2[,7]
        if( sum(INDX) > 0 ){
            XYnew[i,3] = XYobs2[INDX,3]
        }
    }

    XYnewM     = matrix(NA,ncol=8, nrow=length(XYnew[,1]))
    XYnewM[,1] = XYnew[,3]
    XYnewM[,2] = rep("",length(XYnew[,1]))
    XYnewM[,3] = rep(0 ,length(XYnew[,1]))
    XYnewM[,4] = rep(0 ,length(XYnew[,1]))
    XYnewM[,5] = rep(0 ,length(XYnew[,1]))
    XYnewM[,6] = rep(0 ,length(XYnew[,1]))
    XYnewM[,7] = rep("",length(XYnew[,1]))
    XYnewM[,8] = rep(0 ,length(XYnew[,1]))

    for( c in 1:length(cells) ){
        
        tempi = as.numeric(XYobs[,3]) == cells[c]
        temp  = DF[tempi,]
        
        ii = which(XYnewM[,1]==cells[c])

        for( j in 1:length(temp[,1])){
            XYnewM[ii[j],2] = temp[j,1]
            XYnewM[ii[j],3] = temp[j,2]
            XYnewM[ii[j],4] = temp[j,3]
            XYnewM[ii[j],5] = temp[j,4]
            XYnewM[ii[j],6] = temp[j,5]
        }
    
    }
    #-----------

    XYnewM[,2] = gsub(sprintf("%s.",searchTERM[s]), sprintf("%s:",searchTERM[s]),XYnewM[,2])
    XYnewM[,5] = ifelse( (as.numeric(XYnewM[,4]) > 1), as.numeric(XYnewM[,5]),1)
    XYnewM[,6] = ifelse( (as.numeric(XYnewM[,4]) > 1), as.numeric(XYnewM[,6]),1)

    if( searchTERM[s] == "DOID" ){
        XYnewM[,7] = dtype[match(XYnewM[,2],disn)]
    }

    if( searchTERM[s] == "GO" ){
        XYnewM[,7] = XYnewM[,2]
    }
    
    XYnewM[,7] = ifelse(is.na(XYnewM[,7]),"",XYnewM[,7])
    XYnewM[,8] = ifelse(XYnewM[,7]==TARGET[s],1,0)
    XYnewM[,8] = ifelse(is.na(XYnewM[,8]),0,XYnewM[,8])

    XYnewM <- cbind(XYnewM, XYnew[,1], XYnew[,2])   
    
    RES[[b]] = XYnewM
    names(RES)[b] = sprintf("%s cells (%dx%d)", algName[aln], BINS[b], BINS[b])

    rm(cells, df, DF, XYbin, XYgrid, XYgcor, XYobs, XYnew, XYnewM)
}




#Create 3d plot
THETA=40#80
PHI=20
ZLIM=20

LABELalpha <- vector(length=2)
LABELalpha[1] = 0.2
LABELalpha[2] = 1


k=1
p=1

while( k <= length(names(RES)) ){

png(filename=sprintf("%s/%s_%s_%s_rasterPlots_page%d.png",plotDIR,subDIR[S],algName[aln],ANNO[s],p), width=WIDTH, height=HEIGHT, units="px")

par(mfrow=c(2,2))

for( r in 1:4 ){

    X = as.numeric(RES[[k]][,9])
    Y = as.numeric(RES[[k]][,10])

    #cut off
    if( ANNO[s] == "topOntoOVG" ){
        LOGalpha=-log10(alpha[1])
    } 

    if( ANNO[s] == "GOBP" ){
        LOGalpha=-log10(alpha[1])
        LABELalpha[1] = 0.3
    } 

    #p.adusted values
    Z    = -log10(as.numeric(RES[[k]][,6]))
    Zmax = ( (ZLIM/100) * max(Z) ) + max(Z)
    if( Zmax < 6 ){ Zmax = 6 }
    

    
    #which disease to colour
    COLS = as.numeric(RES[[k]][,8])
    
    #labels to plot
    LABELS = ifelse( (Z>=LOGalpha), RES[[k]][,7],"")
    #LABELS = ifelse( (Z>=LOGalpha) | (RES[[k]][,7] == TARGET), RES[[k]][,7],"")

    panelfirst <- function(pmat,NBINS=2){

        x <- y <- seq(0, 1, length=(NBINS+1))
        z <- outer(x,y, FUN = function(x,y) 0*x + 0*y + 0 )
        persp3D(x,y,z, col='#7570B3', add=TRUE, shade=0.75, facets=NA, border="grey10")
 
    
}

scatter3D(x=X,y=Y,z=Z, pch="", type="h", lty=1, lwd=1, xlim=c(0,1), ylim=c(0,1), zlim=c(0,Zmax), xlab=Cent[ce], ylab="Bridgeness", zlab="-log10(p.value)", colkey=FALSE, panel.first=panelfirst(NBINS=BINS[k]), theta=THETA, phi=PHI, colvar=COLS, col=ramp.col(col=c("grey50","darkorange2"),n=2,alpha=c(LABELalpha[1],LABELalpha[2])), ticktype = "detailed", bty="u", col.panel=c("white"), main=names(RES)[k], plot=TRUE)

text3D(x=X, y=Y, z=Z, labels=LABELS, adj=0.5, font=2, add=TRUE, colvar=COLS, col=ramp.col(col=c("grey30","darkorange2"),n=2,alpha=c(LABELalpha[1],LABELalpha[2])), colkey = FALSE, cex = 1.0, plot=TRUE)

x <- y <- seq(0, 1, length.out=2)#length=(1+BINS[r]))#15)
z <- outer(x,y, FUN = function(x,y,VAL=LOGalpha) 0*x + 0*y + VAL )
persp3D(x,y,z, col='#7570B3', add=TRUE, shade=0.75, facets=NA, border="grey10", plot=TRUE)

    k=k+1
    
}

    dev.off()
    p=p+1

}


#---------------------


