library(igraph);

#Overlap of Disease A and B in the interactome
diseaseOverlap <- function(GG, disA, disB, OO){

#disease A genes 
indA  <- which(colnames(OO)==disA)    
IDS1  <- as.character(OO[OO[,indA[1]]!=".",1])
#INDX1 <- match(IDS1,V(GG)$name)
#INDX1 <- ifelse(is.na(INDX1),FALSE,TRUE)
#IDS1  <- IDS1[INDX1]
NIDS1 <- length(IDS1)

#disease B genes 
indB  <- which(colnames(OO)==disB)    
IDS2  <- as.character(OO[OO[,indB[1]]!=".",1])
#INDX2 <- match(IDS2,V(GG)$name)
#INDX2 <- ifelse(is.na(INDX2),FALSE,TRUE)
#IDS2  <- IDS2[INDX2]
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
set.seed(as.numeric(SEED))

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
GNS <- V(gg)$name

#---HDO ID DISEASES of INTEREST
#---HDO Disease short names
source('annotationTYPES.R')


#parent terms are 12 and 13
#pHDOID <- disn[c(12,13)]
#pHDO   <- dtype[c(12,13)]

#remove parent terms from original list
#disn  <- disn[-c(12,13)]
#dtype <- dtype[-c(12,13)]

#load gda's for the ppi network
GDA <- V(gg)$TopOntoOVG
##GDA <- V(gg)$TopOntoOVPAPERS #See p140CAP paper for details

#find those gene ids connected to parent terms
#pgns <- list()
#for( i in 1:length(pHDO) ){
#     pgns[[i]]      <- V(gg)$name[grepl(pHDO[i],GDA,fixed=T)]
#     names(pgns)[i] <- pHDO[i] 
#}

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

oo <- matrix(0,nrow=length(GNS),ncol=(length(dtype)+2))
colnames(oo) <- c("Gene.ID","Gene.Name",dtype)
oo[,1] <- V(gg)$name
oo[,2] <- V(gg)$GeneName

N     <- length(dtype)
NELE  <- N*(N+1)/2

#mean of Disease-disease overlap
DABm <- matrix(0,nrow=NELE, ncol=3)

#SD of Disease-disease overlap
DABsd <- matrix(0,nrow=NELE, ncol=3)

#temp disease-disease overlap storage
DDints <- matrix(0, nrow=NELE, ncol=(2+ITS))

for( k in 0:(NELE-1) ){

  #--- indexing for symmetric matrix
  i = floor( (2*N+1 - sqrt( (2*N+1)*(2*N+1) - 8*k ))/2 );
  j = k - N*i + i*(i-1)/2;

  i = i + 1;
  j = j + i;

  DABm[(k+1),1] <- dtype[i]
  DABm[(k+1),2] <- dtype[j]

  DABsd[(k+1),1] <- dtype[i]
  DABsd[(k+1),2] <- dtype[j]

  DDints[(k+1),1] <- dtype[i]
  DDints[(k+1),2] <- dtype[j]

 }


##start timing
ptm <- proc.time()

cat("Running...\n")
for( p in 1:ITS ){

      cat("Permutation No: ", p, "\n")

     temp <- matrix(".",nrow=length(GNS),ncol=(length(dtype)+2))
     colnames(temp) <- c("Gene.ID","Gene.Name",dtype)
     temp[,1] <- V(gg)$name
     temp[,2] <- V(gg)$GeneName

     for( d in 1:length(dtype) ){

         #get gene-disease associations (GDAs)
    	 IDS <- V(gg)$name[grepl(dtype[d],GDA,fixed=T)]
    	 NN  <- length(IDS)

	 ###Central Nervous System Parent Term
	 ##if( dtype[d] == "AD" || dtype[d] == "PD" || dtype[d] == "HD" || dtype[d] == "Epi" ){
	 ##PGNS <- pgns[[2]]
	 ##} else {
	 ###Disease of Mental Health Parent Term
         ##PGNS <- pgns[[1]]
	 ##}

         #permute the N GDA's relative to all gene ids
    	 rm(IDS)
	 ##IDS <- permute(PGNS, NN)#restrict to parent HDO ids
	 IDS <- permute(GNS, NN) #case

	 XX=igraph::shortest.paths(gg,IDS,IDS,weights=NA)
         diag(XX) = NA
         ds   = apply(XX,1,min,na.rm=T)
         indX = match(names(ds),temp[,1])
         temp[indX,(2+d)] = as.vector(ds)
	 oo[indX,(2+d)] = as.numeric(oo[indX,(2+d)]) + as.numeric(ds)

	 rm(XX)

         #ds  <- rep(0,NN)
    	 #for( i in 1:NN ){
         #
         #ds[i] <- min(as.vector(shortest.paths(gg,IDS[i],IDS[-i],weights=NA)))
         #
         #indX             <- which(temp[,1]==IDS[i])
         #temp[indX,(2+d)] <- as.character(ds[i])
         #
         #oo[indX,(2+d)] <- as.numeric(oo[indX,(2+d)]) + as.numeric(ds[i])
    #}

  }
  

  for( k in 1:NELE ){

    overlap = 0;
    if(DDints[k,1] != DDints[k,2]){
       overlap <- as.numeric(diseaseOverlap(gg,DDints[k,1],DDints[k,2],temp))
    }

    DDints[k,(2+p)] <- overlap

  } 

  }#permutations

pet <- proc.time() - ptm
cat("Finished! ", sprintf("time = %.3f \n", pet[[1]]), "\n")

 for( k in 1:NELE ){

     if(DDints[k,1] != DDints[k,2]){  
        DABm[k,3]  <- as.character(mean(as.numeric(DDints[k,3:(2+ITS)])))
        DABsd[k,3] <- as.character(sd(as.numeric(DDints[k,3:(2+ITS)])))
        } else {
	  DABm[k,3]  <- 0
          DABsd[k,3] <- 0
	}

}

outfile <- file("gene_disease_separation.csv","w")
write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("sAB_random_separation.csv","w")
write.table(DDints, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("mean_disease_separation.csv","w")
write.table(DABm, file=outfile, append=T, row.names=T, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("sd_disease_separation.csv","w")
write.table(DABsd, file=outfile, append=T, row.names=T, col.names=T, sep="\t",quote=F);
close(outfile);



