source('../setUp.R')

#Check file headers for additional formatted columns to remove
formatFILES <- function(FILE){

    #Check file columns for following
    del <- vector(length=2)
    del[1] <- FALSE
    del[2] <- ""

    for( d in 1:length(del) ){
        keep           <- colnames(FILE) != del[d]     
        CN             <- colnames(FILE)[keep]    
        FILE           <- FILE[,keep]
        colnames(FILE) <- CN   
    }

    return(FILE)

}


 createSheet <- function( anno, study, root ){

     rm(tt1,tb)
     tt1 = NULL
     st1 = sprintf("%s/%s/permute_p_values_%s.csv",root,study[1],anno)

     if( file.exists(st1) && file.info(st1)$size!=0 ){
     
     tt1 <- read.table(st1,sep="\t",header=T,quote="",check.names=F,stringsAsFactors=F,na.strings=c("NA"," ",""))

     #Check file header
     tt1 <- formatFILES(tt1)        

     rn <- rep("",length=length(tt1[1,]))
     rn[1] = study[1]
     tt1 <- rbind(rn,tt1[1:length(tt1[,1]),])
  
     blank <- rep("",length=length(tt1[1,]))

     if( length(study) > 1 ){
     
         for( i in 2:length(study) ){

                          
             tb <- read.table(sprintf("%s/%s/permute_p_values_%s.csv",root,study[i],anno),sep="\t",header=T,quote="",check.names=F,stringsAsFactors=F,na.strings=c("NA"," ",""))

             #Check file header
             tb <- formatFILES(tb)        
             
             rn <- rep("",length=length(tb[1,]))
             rn[1] = study[i]

             if( !identical(names(tt1),names(tb)) ){
                 names(tb) <- names(tt1)
             }
             
             tb <- rbind(blank,rn,tb[1:length(tb[,1]),])
             
             tt1 <- rbind(tt1,tb)

         }
     }

     }
     
  return(tt1)
  
}


addClusterEntropy <- function( TT, FILE ){


    ent <- read.table(FILE,sep="\t",header=T,quote="",stringsAsFactors=F)

    MM <- vector(length=1)
    MM[1] <- "Robustness"
    
    CN  <- colnames(TT)[colnames(TT) != ""]

    tmp <- matrix("",nrow=length(TT[,1]),ncol=(length(CN)+1))
    colnames(tmp) <- c(CN[1],MM[1],CN[2:length(CN)])
    
    TT[is.na(TT[,1]),1] <- ""
    tmp[,1] <- as.character(TT[,1])
    for( i in 2:length(CN) ){
        TT[is.na(TT[,i]),i] <- ""
        tmp[,(i+length(MM))] <- as.character(TT[,i])
    }


    alg <- unique(ent[,1])

    for( a in 1:length(alg) ){
    
        indxA <- which(tmp[,1]==alg[a])

        if( length(indxA) != 0 ){
        
            indx  <- which(ent[,1]==alg[a])    

            startI <- (indxA+1)
            endI   <- indxA+length(indx)
        
            tmp[startI:endI,2] <- as.numeric(ent$Crob[indx])#Cluster Robustness
            

        }
            
    }

    tmp <- as.data.frame(tmp)
    
    return(tmp)
    
}

#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("Consensus",DIRS)]
OUT[3] <- DIRS[grepl("EnrichmentPackage",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

#---declare clustering algorithms in graph, and with a corresponding consensus matrix
#alg    <- ALGS[c(1:10)]

calgs <- read.table(sprintf("%s/clusteringAlg.csv",pramFILES),header=F,sep="\t",quote="")

alg   <- calgs[calgs[,1] == 1 & match(calgs[,3],ALGS[c(1:10)]),3]

#aa    <- match(as.vector(calgs[,3]),alg)
#alg   <- alg[aa]


#Don't use for interPro annotation, too many singlet communities
alg2 <- alg[!grepl("wt",alg,fixed=T)]

#--- set of annotation files
Anno  <- Anno[Anno!="Cplx"]
Anno2 <- c("Family","Domain","GOBP")

ROOT <- sprintf("%s/RESULTS/%s",OUT[3],subDIR[S])

N <- length(Anno)

tt <- list()
k=1
for( a in 1:N ){

    if(!is.na(match(Anno[a],Anno2))){
        temp <- createSheet(Anno[a],alg2,ROOT)
    } else{
        temp <- createSheet(Anno[a],alg,ROOT)
    }

    if( !is.null(temp) ){
        tt[[k]]      <- temp
        names(tt)[k] <- Anno[a]
        k=k+1
    }
}
    
#--- test, add Cluster Entropy info
#for( a in 1:N ){
#    str <- sprintf("%s/%s/%s_ClustersEntropy.csv",OUT[2],subDIR[S],subDIR[S])    
#    tt[[a]] <- addClusterEntropy( tt[[a]], str )
#}
#---

#---Annotation and Clustering results
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

colNames <- names(igraph::vertex.attributes(gg))

cc <- matrix(NA,ncol=length(colNames),nrow=length(V(gg)))

for( i in 1:length(colNames) ){
    cc[,i] <- as.character(igraph::get.vertex.attribute(gg,names(igraph::vertex.attributes(gg))[i],V(gg)))
}

colNames[which(colNames=="name")]="Entrez.ID"
colnames(cc) <- colNames

cc <- as.data.frame(cc)
#---

#---ADD CENTRALITY MEASURES
str <- sprintf("%s/%s/%s_Measures.csv",OUT[2],subDIR[S],subDIR[S])

meas <- read.table(str,sep="\t",header=T)

indX <- match(cc[,2],meas[,1])

#Don't need to add these columns from meas
rem <- c(1,2,match(alg,colnames(meas)))

cc <- cbind(cc,meas[indX,-rem])
#---

#Make sure GeneName comes after Entrez.ID
str   <- "GeneName"
Cindx <- grep(str, names(cc))
tmp   <- as.data.frame(cc[,Cindx])
colnames(tmp) <- str
cc    <- cc[,-Cindx]
cc    <- cbind(cc[,c(1,2)],tmp,cc[,c(3:ncol(cc))])
#---

#create a final spread-sheet list
ff <- list()

ColMAX <- 350

ff[[1]] <- as.data.frame(cc)
names(ff)[1] <- "Communities"
k=2
for( f in 1:length(names(tt)) ){

    if( length(colnames(tt[[f]])) < ColMAX ){
        ff[[k]] = tt[[f]]
        names(ff)[k] <- Anno[f]
    }

    k=k+1
}
#---

#Write cluster enrichment to spread-sheet
WriteXLS(x=ff,SheetNames=names(ff),row.names=F,col.names=T,ExcelFileName=sprintf("%s/%s.xls",subDIR[S],subDIR[S]))

#Write proDomains file separately (too big) and add later to spread-sheet
for( f in 1:length(names(tt)) ){

    if( length(colnames(tt[[f]])) > ColMAX ){
        write.table(tt[[f]],sprintf("%s/%s_%s.csv",subDIR[S],subDIR[S],names(tt)[f]),row.names=F,col.names=T,sep="\t",quote=F)
    }
}
