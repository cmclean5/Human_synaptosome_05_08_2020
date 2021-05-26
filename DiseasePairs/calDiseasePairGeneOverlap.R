source('../setUp.R')
require(cowplot)

#Overlap of Disease A and B in the interactome
# GG   => igraph network
# GDA  => gda data for this graph
# disA => name of disease A
# disA => name of disease B
diseaseGeneOverlap <- function(GG, GDA, disA, disB){

    #network genes
    N = length(V(GG))
    
    #disease A genes 
    IDS1  <- V(GG)$name[grepl(disA,GDA,fixed=T)]
    NIDS1 <- length(IDS1)

    #disease B genes 
    IDS2  <- V(GG)$name[grepl(disB,GDA,fixed=T)]
    NIDS2 <- length(IDS2)

    Cobs  = sum(IDS1 %in% IDS2)
    Cmin  = ifelse( (NIDS1 > NIDS2), NIDS2, NIDS1)
    Cuin  = length(union(IDS1,IDS2))

    Crand  = ifelse(N    > 0, ((NIDS1 * NIDS2) / N), NA)
    C      = ifelse(Cmin > 0, (Cobs / Cmin), NA)
    J      = ifelse(Cuin > 0, (Cobs / Cuin), NA)

    #stat    = dhyper( (0:Cobs), NIDS1, (N-NIDS1), NIDS2) 
    #pvalueL = sum(stat[stat <= Cobs])
    #pvalueU = sum(stat[stat > Cobs]) 
    pvalueL = phyper( Cobs, NIDS1, (N-NIDS1), NIDS2, lower.tail = TRUE  ) 
    pvalueU = phyper( Cobs, NIDS1, (N-NIDS1), NIDS2, lower.tail = FALSE ) 
    
    res <- vector(length=6)
    res[1] = Cobs
    res[2] = C
    res[3] = J
    res[4] = Crand
    res[5] = pvalueL
    res[6] = pvalueU

    names(res) = c("Cobs","C","J","Crand","pvalueLOW","pvalueUPPER")
    
    return(res)
    
}


#---Directories
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

dadir <- sprintf("%s/%s",subDIR[S],gdaDIR[gdas])
if( !file_test("-d",dadir) ){
    dir.create(dadir)
}
#---

#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

plotdir <- sprintf("PLOTS/%s",subDIR[S])

if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}

FILES    <- vector(length=2)
FILES[1] <- "Disease_overlap_sig"
FILES[2] <- "random_zscores"

#--- load results
ff <- read.table(sprintf("RESULTS/%s/ovg/%s.csv",subDIR[S],FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

#--- load randomised zscores
tests <- read.table(sprintf("RESULTS/%s/ovg/%s.csv",subDIR[S],FILES[2]),sep="\t",header=F,stringsAsFactors=F,quote="", check.names=F)


#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

#--- Find all Gene Disease Associations
GDA <- V(gg)$TopOntoOVG

if( gdas == 2 ){
    GDA <- V(gg)$TopOntoOVPAPERS #See '../setUp.R' for details
}

#--- The number of GDA's in graph
NN  <- length(which(GDA!=""))


#---Check
#---make sure we remove any parent HDO terms first
indx  <- match(pHDO,dtype)
if( length(indx) > 0 ){
    if( !is.na(indx) ){
        disn  <- disn[-indx]
        dtype <- dtype[-indx]
        disl  <- disl[-indx]
    }
}

#---Check
#---Remove Diseases with zero GDA's
remove <- c()
cat("Following Diseases have zero GDA's.\n")
for( d in 1:length(dtype) ){
    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    if( length(IDS) == 0 ){
        cat(dtype[d], " => ", length(IDS),"\n")
	remove <- c(remove,d)
     }
}

if( length(remove) == 0 ){
    cat(" => None.\n")
}

if( length(remove) > 0 ){
   disn  <- disn [-remove]	
   dtype <- dtype[-remove]
   disl  <- disl[-remove]
}    
#---     

OO = matrix(NA, nrow=length(ff[,1]),ncol=8)
OO[,1] = ff[,3]
OO[,2] = ff[,7]
OO[,3] = ff[,9]

for( i in 1:length(ff[,1]) ){

    if( OO[i,1] != OO[i,2] ){
        RES <- diseaseGeneOverlap(gg,GDA,OO[i,1],OO[i,2])

        OO[i,4] = RES[2]
        OO[i,5] = RES[3]
        OO[i,6] = RES[1]/RES[4]
        OO[i,7] = RES[5]
        OO[i,8] = RES[6]
    } else {
        OO[i,1] = NA
        OO[i,2] = NA
        OO[i,3] = NA
    }

}

#remove self disease pairs
OO = OO[!is.na(OO[,1]),]

#log(Jaccard) V. sAB
plot(xy.coords(x=as.numeric(OO[,3]),y=log(as.numeric(OO[,5]))))



