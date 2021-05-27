library(igraph);
library(data.table);

#--- set random number seed
setSeed <- function(SEED=NULL, useTIME=TRUE){

    # make sure we use sample‘s v3.6.0 default behavior
    RNGkind(kind="default")    
    
    if( (is.null(SEED) && useTIME) || (is.null(SEED) && !useTIME) ){
        SEED = as.integer(Sys.time())
        set.seed(SEED)
    } else if ( !is.null(SEED) && useTIME ){
        SEED = SEED + as.integer(Sys.time())
        set.seed(SEED)    
    } else if ( !is.null(SEED) && !useTIME) {
        set.seed(SEED)
    }


}

#---Find Largest CC
findLCC <- function(GG){

    dec <- igraph::decompose.graph(GG)
    d=1
    CC=length(V(dec[[1]]))
    for( i in 1:length(dec) ){
        if(length(V(dec[[i]])) > CC){
            d=i
            CC=length(V(dec[[i]]))
        }
    }   

    GG  <- igraph::decompose.graph(GG)[[d]]
    return(GG)

}

generateRandomEdges <- function( GG, errorRate=seq(0,1,0.1), SEED=NULL ){

    #setSeed(SEED)
    
    Ntot  <- length(V(GG))
    Mtot  <- length(E(GG))

    Mrand <- floor(Mtot * errorRate)
    
    rndmSET <- list()    

    for( i in 1:length(Mrand) ){

        if( Mrand[i] == 0 ){
            rd <- GG
        #} else if( Mrand[i] == Mtot ){
        #    rd <- sample_gnm(n=Ntot,m=Mtot,directed=F,loops=F)
        } else {

            #--- introduces false negative errors, i.e. 
            #    randomly remove edges from graph
            rd <- delete_edges(GG,sample(E(GG),Mrand[i]))
    
            
            #--- introduces false postive errors, i.e. 
            #    spurious edges (i.e. where none existed before) into the graph
            rd <- rd + igraph::edges(sample(V(rd),2*Mrand[i],replace=T))
                
        }	

	rd <- simplify(rd, remove.multiple=TRUE, remove.loops=TRUE)	

	steps =	0
        while( length(E(rd)) < Mtot && steps < 10 ){

           Diff = Mtot - length(E(rd))                
           rd <- rd + igraph::edges(sample(V(rd),2*Diff,replace=T))           
           rd <- simplify(rd, remove.multiple=TRUE, remove.loops=TRUE)

           steps = steps + 1
         }  

	 
	 #---Find Largest CC
	rndmSET[[i]] <- findLCC(rd)
        ##rndmSET[[i]] <- rd
        names(rndmSET)[i] <- sprintf("errRate_%g", errorRate[i])

    }      
   
    return(rndmSET)
    
}


#Overlap of Disease A and B in the interactome
diseaseOverlap <- function(GG, disA, disB, OO){

#disease A genes 
indA  <- which(colnames(OO)==disA)    
IDS1  <- as.character(OO[OO[,indA[1]]!=".",1])
INDX1 <- match(IDS1,V(GG)$name)
INDX1 <- ifelse(is.na(INDX1),FALSE,TRUE)
IDS1  <- IDS1[INDX1]
NIDS1 <- length(IDS1)

#disease B genes 
indB  <- which(colnames(OO)==disB)    
IDS2  <- as.character(OO[OO[,indB[1]]!=".",1])
INDX2 <- match(IDS2,V(GG)$name)
INDX2 <- ifelse(is.na(INDX2),FALSE,TRUE)
IDS2  <- IDS2[INDX2]
NIDS2 <- length(IDS2)

    #disease A given B
    paths  <- igraph::shortest.paths(GG,IDS1,IDS2,weights=NA)
    dsA    <- as.numeric(as.vector(apply(paths,1,min)))	
    
    #dsA <- rep(0,NIDS1)
    #for( i in 1:NIDS1 ){
    #
    #    paths  <- as.vector(shortest.paths(GG,IDS1[i],IDS2,weights=NA))
    #    dsA[i] <- 0
    #    
    #    if( (0 %in% paths) == FALSE  ){
    #        dsA[i] <- as.numeric(min(paths) )
    #    }
    #    
    #}#end
    
    #disease B given A
    paths <- igraph::shortest.paths(GG,IDS2,IDS1,weights=NA)
    dsB   <- as.numeric(as.vector(apply(paths,1,min)))
 	     
    #dsB <- rep(0,NIDS2)
    #for( i in 1:NIDS2 ){
    #
    #    paths  <- as.vector(shortest.paths(GG,IDS2[i],IDS1,weights=NA))        
    #    dsB[i] <- 0  
    #
    #    if( (0 %in% paths) == FALSE  ){
    #        dsB[i] <- as.numeric(min(paths) )
    #    }                
    #}#end

    #network-based separation between disease A and B 
    dAB <- (sum(dsA)+sum(dsB))/(NIDS1+NIDS2)

    #network-based localisation of disease A
    dA   <- mean(as.numeric(as.vector(OO[OO[,indA[1]]!=".",indA[1]])))

    #network-based localisation of disease B
    dB   <- mean(as.numeric(as.vector(OO[OO[,indB[1]]!=".",indB[1]])))
    
    #overlap between disease A and B
    sAB = as.numeric(dAB) - (as.numeric(dA)+as.numeric(dB))/2

    return(sAB)
        
}

# permute biological condition data, i.e
# the sample (given by the cols) and genes
# (given by the rows).
# This needs to be run on 'eddie' ECDF, since ~10,000 iterations needed.
permute <- function(GNS, N){

	temp <- sample(GNS,N,replace=F)

	return(temp)
    
}

#---Script required inputs
args  <- commandArgs(TRUE);
SEED  <- as.numeric(args[1]) #random seed no. 
ITS   <- as.numeric(args[2]) #no: of iterations

#SEED <- 1095855
#ITS  <- 1

#set the random number seed
#set.seed(as.numeric(SEED))
setSeed(SEED);

cat("Running DiseaseLoc.R: \n")
cat("SEED: ", SEED, "\n")
cat("Permutations: ", ITS, "\n")

#---YOUR DATASET TO READ IN
files   <- list.files()
files   <- files[grepl(".gml" ,files,fixed=T)]

studies <- unlist(strsplit(files,".gml"))

#--PPI networks
gg <- read.graph(files[1],format="gml")

#Get all gene Entrez IDS
#GNS <- V(gg)$name

#---HDO ID DISEASES of INTEREST
#---HDO Disease short names
source('annotationTYPES.R')

#load gda's for the ppi network
GDA <- V(gg)$TopOntoOVG


#---Remove Diseases with zero GDA's
remove <- c()
for( d in 1:length(dtype) ){
    IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    if( length(IDS) == 0 ){
        cat(dtype[d], " => ", length(IDS),"\n")
	remove <- c(remove,d)
     }
}
if( length(remove) > 0 ){
   disn  <- disn [-remove]	
   dtype <- dtype[-remove]
}    
#---  

#oo <- matrix(".",nrow=length(GNS),ncol=(length(dtype)+2))
#colnames(oo) <- c("Gene.ID","Gene.Name",dtype)
#oo[,1] <- V(gg)$name
#oo[,2] <- V(gg)$GeneName

#for( d in 1:length(dtype) ){
#
#    #get gene-disease associations (GDAs)
#    IDS  <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
#    for( i in 1:length(IDS) ){
#        indX <- which(oo[,1]==IDS[i])
#        oo[indX,(2+d)] = 1
#    }
#}

# define error Rate
errorRate=seq(0,1,0.1)
###errorRate=c(1)

N     <- length(dtype)
NELE  <- N*(N+1)/2

DDINTS <- list() 

for( e in 1:length(errorRate) ){

   #temp disease-disease overlap storage
   DDints <- matrix(0, nrow=NELE, ncol=(2+ITS))

    for( k in 0:(NELE-1) ){

        #--- indexing for symmetric matrix
        i = floor( (2*N+1 - sqrt( (2*N+1)*(2*N+1) - 8*k ))/2 );
        j = k - N*i + i*(i-1)/2;

        i = i + 1;
        j = j + i; 

        DDints[(k+1),1] <- dtype[i]
        DDints[(k+1),2] <- dtype[j]

    }

    DDINTS[[e]] = DDints
    names(DDINTS)[e] = sprintf("%g",round(errorRate[e],3))

}


##start timing
ptm <- proc.time()

cat("Running...\n")
for( p in 1:ITS ){

    cat("Permutation No: ", p, "\n")    

    cat("generating random graphs...")
    # generate random edges
    rd = generateRandomEdges(gg,errorRate=errorRate) 
    cat("done!\n")

    cat("calculate Disease pair overlaps...\n")
    
    for( e in 1:length(rd) ){

        cat(" for errorRate: ", errorRate[e])
        
	temp           <- matrix(".",nrow=length(V(rd[[e]])),ncol=(length(dtype)+2))
	colnames(temp) <- c("Gene.ID","Gene.Name",dtype)
        temp[,1]       <- V(rd[[e]])$name
        temp[,2]       <- V(rd[[e]])$GeneName

        for( d in 1:length(dtype) ){

         #get gene-disease associations (GDAs)
         rdGDA <- V(rd[[e]])$TopOntoOVG
    	 IDS   <- V(rd[[e]])$name[grepl(dtype[d],rdGDA,fixed=T)]
         INDX  <- match(IDS,V(rd[[e]])$name)
         INDX  <- ifelse(is.na(INDX),FALSE,TRUE)
         IDS   <- IDS[INDX]
    	 NN    <- length(IDS)         

	 XX=igraph::shortest.paths(rd[[e]],IDS,IDS,weights=NA)
         diag(XX) = NA
         ds   = apply(XX,1,min,na.rm=T)
         indX = match(names(ds),temp[,1])
         temp[indX,(2+d)] = as.vector(ds)

	 rm(XX)

         #ds  <- rep(0,NN)
    	 #for( i in 1:NN ){

	 ###paths <- as.vector(shortest.paths(rd[[e]],IDS[i],IDS[-i]))
	 ###ds[i] <- 0
         #ds[i] <- min(as.vector(shortest.paths(rd[[e]],IDS[i],IDS[-i],weights=NA)))
         #
	 ###if( (0 %in% paths) == FALSE ){
	 ###  ds[i] = as.numeric(min(paths))
         ### }
         #
         #indX <- which(temp[,1]==IDS[i])
         #temp[indX,(2+d)] <- as.character(ds[i])
        
    }

  }
  
        for( k in 1:NELE ){

            overlap = 0;
            if(DDints[k,1] != DDints[k,2]){
                overlap <- as.numeric(diseaseOverlap(rd[[e]],DDints[k,1],DDints[k,2],temp))
            }
            
            DDINTS[[e]][k,(2+p)] <- overlap

        }#k

        cat(" done!\n")       
    
    }#errorRate

    cat("done!\n")
        
    ## save output files 	
    for( e in 1:length(errorRate) ){
       st1 = sprintf("%s_%g.csv.gz","sAB_random_separation",errorRate[e])
       res = as.data.table(DDINTS[[e]])	   
       fwrite(res,st1)   
    }

}#permutations

pet <- proc.time() - ptm
cat("Finished! ", sprintf("time = %.3f \n", pet[[1]]), "\n")




