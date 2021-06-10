source('../setUp.R')

setClass(Class="DD",representation(
                        II="matrix",
                        EE="matrix",
                        Tr="vector"
                    )
         )

inAnnoSet   <- function( ANNO, targetSET ){

    if( length(targetSET) == 0 ){
        return( (ANNO != "") )
    } else {

        indx = grepl(targetSET[1], ANNO)
        if( length(targetSET) > 1 ){
            for( i in 2:length(targetSET) ){
                indx = indx & grepl(targetSET[i], ANNO)
            }
        }
        return(indx)
    }

    return(-1)
    
}

AnnoDensity <- function( GR, AG, EDag, ANNO, EDanno, NORM=FALSE, targetSET=NULL, MAXnorm=FALSE ){

    N  <- length(GR)
    ii <- matrix(0, ncol=N, nrow=N)
    ee <- matrix(0, ncol=N, nrow=N)

    for( i in 1:N ){

        Ni = sum( (AG==GR[i]) & (inAnnoSet(ANNO, targetSET)))

        if( MAXnorm ){
            Ni = sum( (AG==GR[i]) )
        }
        
        for( j in i:N ){

            Nj = sum( (AG==GR[j]) & (inAnnoSet(ANNO, targetSET)))

            if( MAXnorm ){
                Nj = sum( (AG==GR[j]) )
            }   
            
            E = Ni * Nj

            if( i == j ){
                E = (Ni*(Nj-1))/2
            }

            if( E == 0 ){ E = 1 }
            
            indx = (EDag[,1] == GR[i] & EDag[,2] == GR[j]) | (EDag[,2] == GR[i] & EDag[,1] == GR[j])

            xx      = EDanno[indx,]

            if( !is.null(targetSET) ){
                
                tit = paste(targetSET,collapse=";")
                
                indxA <- inAnnoSet(xx[,1], targetSET)
                indxB <- inAnnoSet(xx[,2], targetSET)

                xx[,1] <- ifelse(indxA, tit, "")
                xx[,2] <- ifelse(indxB, tit, "")
                
            }

            ee[i,j] = E
            
            if( NORM ){            
                ii[i,j] = jacardIndex(xx) / E
            } else {
                ii[i,j] = jacardIndex(xx)
            }

        }
        
    }

    if( is.null(targetSET) ){
        targetSET <- c("")
    }
    
    return(new("DD",II=ii,EE=ee, Tr=targetSET))
    
}



jacardIndex <- function( ED ){

    ji = 0
    
    if( length(ED) != 0 ){
    
        if( !is.vector(ED) ){
    
            N = length(ED[,1])
        
            for( i in 1:N ){

                if( ED[i,1] != "" && ED[i,2] != "" ){ 
                
                    num = length(intersect(strsplit(ED[i,1],";")[[1]],strsplit(ED[i,2],";")[[1]]))

                    dem = length(union(strsplit(ED[i,1],";")[[1]],strsplit(ED[i,2],";")[[1]]))

                    if( dem != 0 ){
                        ji = ji + num/dem
                    }
                }

            }
            
        } else {

            if( ED[1] != "" && ED[2] != "" ){ 
            
                num = length(intersect(strsplit(ED[1],";")[[1]],strsplit(ED[2],";")[[1]]))

                dem = length(union(strsplit(ED[1],";")[[1]],strsplit(ED[2],";")[[1]]))
                
                if( dem != 0 ){
                    ji = ji + num/dem
                }
            }
        }
    }
    
    return(ji)

}

#---OUT Dir
OUT    <- vector(length=2)
OUT[1] <- DIRS[grepl("Phyogeny",DIRS)]
OUT[2] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
grdir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---


#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")

M   <- length(E(gg))

ids <- V(gg)$name;
ag  <- V(gg)$AgeGroup;
bp  <- V(gg)$GOBPID;
mf  <- V(gg)$GOMFID;
cc  <- V(gg)$GOCCID;
ds  <- V(gg)$TopOntoOVG;
dg  <- V(gg)$Dg;
dg  <- ifelse(dg=="NO","",dg)

Gr = unique(ag)
Gr = Gr[Gr != ""]
N  = length(Gr)

ed  <- get.edgelist(gg)
age <- cbind(ag[match(ed[,1],ids)],ag[match(ed[,2],ids)])

bpe <- cbind(bp[match(ed[,1],ids)],bp[match(ed[,2],ids)])
mfe <- cbind(mf[match(ed[,1],ids)],mf[match(ed[,2],ids)])
cce <- cbind(cc[match(ed[,1],ids)],cc[match(ed[,2],ids)])
dse <- cbind(ds[match(ed[,1],ids)],ds[match(ed[,2],ids)])
dge <- cbind(dg[match(ed[,1],ids)],dg[match(ed[,2],ids)])

indx <- !(age[,1] == "" | age[,2] == "")

age <- age[indx,]
bpe <- bpe[indx,]
mfe <- mfe[indx,]
cce <- cce[indx,]
dse <- dse[indx,]
dge <- dge[indx,]

oo <- matrix(0, ncol=N, nrow=N)

ii <- matrix(0,ncol=N, nrow=N)
ee <- matrix(0,ncol=N, nrow=N)

for( i in 1:N ){

    Ni = sum(ag == Gr[i])

    for( j in i:N ){

        Nj = sum(ag == Gr[j])

        E = Ni * Nj
        
        if( i == j ){
            E = Ni*(Nj-1)/2
        }
        
        I=sum((grepl(Gr[i],age[,1]) & grepl(Gr[j],age[,2])) |
              (grepl(Gr[j],age[,2]) & grepl(Gr[i],age[,1])))

        ii[i,j] = I
        ee[i,j] = E                
        oo[i,j] = (I/E)*100
        
    }
    
}

age2 <- cbind(as.numeric(gsub("G","",age[,1])), as.numeric(gsub("G","",age[,2])))

pairs = sum(abs(age2[,1] - age2[,2]) == 0)

sum(pairs == 0)/length(age2[,1])
sum(pairs == 1)/length(age2[,1])
sum(pairs >= 2)/length(age2[,1])


ANNOres <- list()

ANNOres[[1]]       <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=bp, EDanno=bpe )
names(ANNOres)[1] <- "bp"

ANNOres[[2]]       <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=mf, EDanno=mfe )
names(ANNOres)[2] <- "mf"

ANNOres[[3]]       <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=cc, EDanno=cce )
names(ANNOres)[3] <- "cc"

ANNOres[[4]]       <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=ds, EDanno=dse )
names(ANNOres)[4] <- "ds"

ANNOres[[5]]       <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=dg, EDanno=dge )
names(ANNOres)[5] <- "dg"

DS <- list()
ND <- length(dtype)

   for( d in 1:ND ){
       DS[[d]] <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=ds, EDanno=dse, NORM=FALSE, targetSET=dtype[d], MAXnorm=FALSE )
       names(DS)[d] <- dtype[d]

   }
       
       
#--set disease pair of interest
disA = vector(length=5)
disA[1] = "AD"
disA[2] = "AD"
disA[3] = "AD"
disA[4] = "PD"
disA[5] = "MS"

disB = vector(length=5)
disB[1] = "HTN"
disB[2] = "PD"
disB[3] = "MS"
disB[4] = "HTN"
disB[5] = "HTN"


TT <- list()

for( t in 1:length(disA) ){

    tit = sprintf("%s;%s",disA[t],disB[t])

    TT[[t]] <- AnnoDensity( GR=Gr, AG=ag, EDag=age, ANNO=ds, EDanno=dse, NORM=FALSE, targetSET=c(disA[t], disB[t]), MAXnorm=FALSE )
    names(TT)[t] <- tit

}
    
  
