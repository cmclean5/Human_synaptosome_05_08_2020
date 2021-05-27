#clean
rm(list=ls())

library(igraph);
library(poweRlaw);

args <- commandArgs(TRUE);

NSEED   <- as.numeric(args[1])
NSIM    <- as.numeric(args[2])
NTHREAD <- as.numeric(args[3])
NSTUDY  <- as.numeric(args[4])

set.seed(NSEED)

maxLSi <- function( XX, BASE=0 ){

    XX  <- as.vector(as.numeric(XX))
    
    XXo <- XX[XX != 0]
    if( BASE == 2 ){
        return ( -sum( XXo * log2(XXo)) )
    }

    if( BASE == 10 ){
        return ( -sum( XXo * log10(XXo)) )
    }

    if( BASE == 0 ){
        return ( -sum( XXo * log(XXo)) )
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


graphEntropy <- function(GG){

#--- initial entropy rate
V    <- length(V(GG))
E    <- length(E(GG))
ki   <- as.vector(igraph::degree(graph=GG))
Kbar <- mean(ki)

#--- get adjacency matrix for graph
A    <- igraph::get.adjacency(GG)

#--- get leading eigenvalue and vector
R     <- eigen(A)
Rindx <- which.max(R$values)
gamma <- R$values[Rindx]
nu    <- R$vectors[,Rindx]

#--- calculate max entropy rate, maxSr
Pij   <- (A * nu) / (gamma * nu)
Pi    <- ki/(2*E)
maxSr <- sum(as.vector(Pi * apply(Pij,1,maxLSi,BASE=0)))


#--- calculate initial configuration
Norm <- as.numeric(V*Kbar) #as.numeric(2*E)
SRo  <- as.numeric(1/Norm)*sum(ki*log(ki))

#if( maxSr == 0 ){ maxSr = 1 }

return( SRo/maxSr )
    
}


FitDegree <- function(DEG, Nsim, Nthreads, Nseed ){

   ret <- vector(length=2)
   ret[1] = 0
   ret[2] = 0   	  

     DEG <- DEG[DEG > 0]
     
     data <- DEG     
     
     m_pl = displ$new(data)

     est = estimate_xmin(m_pl)

     m_pl$setXmin(est)

     gof <- bootstrap_p(m_pl,no_of_sims=Nsim, threads=Nthreads, seed=Nseed )
     
     ret[1] = est$pars
     ret[2] = gof$p
    
     return(ret)
     
 }


#---Semi-local Centrality (Cl)
#   Identifying influential nodes in complex networks, D. Chen et al., Physica A, 2012
Semilocal <- function(gg){

N    <- length(V(gg)$name)
meas <- matrix(0, nrow=N, ncol=3)

for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- neighbors(gg,v=ids,mode="all")

    if( length(neig) > 0 ){
  
        for( w in 1:length(neig) ){
            neig <- c(neig,neighbors(gg,v=as.character(V(gg)$name[neig[w]]),mode="all"))    
        }

        neig <- unique(neig)

        meas[i,1] <- length(neig)-1

  }
    
}

for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- neighbors(gg,v=ids,mode="all")
  
    meas[i,2] <- sum(as.numeric(meas[neig,1]))
  
}


for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)
  
    neig <- neighbors(gg,v=ids,mode="all")
  
    meas[i,3] <- sum(as.numeric(meas[neig,2]))
  
}

 return(as.numeric(meas[,3]))

}

#calculate the mean and sd of the shortest paths for each gene
calShorestPaths <- function(gg){

 N    <- length(V(gg)$name)
 meas <- matrix(0, nrow=N, ncol=3)
		
 for( i in 1:N ){
      sp <- as.numeric(shortest.paths(gg,i))
      sp <- sp[-i]
      sp <- sp[!sp == Inf]
      meas[i,1] <- min(sp)
      meas[i,2] <- round(mean(sp),3)
      meas[i,3] <- round(sd(sp),3)
 } 		

 return(meas)

}


 alpha <- matrix(0,ncol=2,nrow=1)
 #alpha[1,1] = "PPI_Presynaptic_Published"; alpha[1,2] = "2.545";
 #alpha[2,1] = "PPI_PSD_clean_Published";   alpha[2,2] = "2.511";
 #alpha[3,1] = "PPI_PSD_reduced";           alpha[3,2] = "2.420";
 #alpha[1,1] = "PPI_PSD_reduced";           alpha[1,2] = "2.59";

 #alpha[1,1] = "PPI_PSP_clean_Published";   alpha[1,2] = "2.55";
 #alpha[2,1] = "PPI_PSP_reduced";           alpha[2,2] = "2.53";

 #alpha[1,1] = "PPI_PSP";                     alpha[1,2] = "2.51";
 #alpha[2,1] = "PPI_PSP_consensus";           alpha[2,2] = "2.85";
 #alpha[3,1] = "PPI_PSP_consensus2";          alpha[3,2] = "2.93";

 alpha[1,1] = "PPI_PSP";                     alpha[1,2] = "2.41";
 

#---Directories
DIR  <- "datasets/"

#---Load original networks
studies  <- list.files(DIR,"*.gml")
Nstudies <- unlist(strsplit(studies,".gml"))

##gsO = setNames(vector("list",length(studies)),nm=Nstudies)
gsO = list()

for( i in 1:length(Nstudies) ){

    gg <- read.graph(sprintf("%s/%s",DIR,studies[i]),format="gml")
    gsO[[i]]      <- gg
    names(gsO)[i] <- Nstudies[i]       

}

ID <- V(gsO[[1]])$name

if( length(names(gsO)) > 1 ){

    for( i in 2:length(names(gsO)) ){
    	  ID <- unique(c(ID,V(gsO[[i]])$name))
       }	  
}

N <- length(ID)


#--- randomise original graphs

rnd <- list()

if( NSTUDY == 1 ){

for( i in 1:length(studies) ){
     rnd[[i]] <- sample_gnm( n=length(V(gsO[[i]])), m=length(E(gsO[[i]])),directed=F,loops=F)#, loops=T)
     rnd[[i]] <- set.vertex.attribute(rnd[[i]],"name",V(rnd[[i]]),V(gsO[[i]])$name)
     rnd[[i]] <- set.vertex.attribute(rnd[[i]],"GeneName",V(rnd[[i]]),V(gsO[[i]])$GeneName)

    rnd[[i]]      <- findLCC(rnd[[i]])
    names(rnd)[i] <- names(gsO)[i]
}

}

if( NSTUDY == 2 ){

    for( i in 1:length(studies) ){

    ind <- which(alpha[,1]==names(gsO)[i])	 
    if( length(ind) != 0 ){
        Nalpha <- as.numeric(alpha[ind[1],2])
    }

    rnd[[i]] <- sample_fitness_pl(
    no.of.nodes = length(V(gsO[[i]])),
    no.of.edges = length(E(gsO[[i]])), 
    exponent.out = Nalpha, 
    exponent.in = -1,
    loops = FALSE, 
    multiple = FALSE, 
    finite.size.correction = TRUE)
    
     rnd[[i]] <- set.vertex.attribute(rnd[[i]],"name",V(rnd[[i]]),V(gsO[[i]])$name)
     rnd[[i]] <- set.vertex.attribute(rnd[[i]],"GeneName",V(rnd[[i]]),V(gsO[[i]])$GeneName)

    rnd[[i]]      <- findLCC(rnd[[i]])
    names(rnd)[i] <- names(gsO)[i]

}

}

if( NSTUDY == 3 ){

    for( i in 1:length(studies) ){

    Nalpha = 2.0

    rnd[[i]] <- sample_fitness_pl(
    no.of.nodes = length(V(gsO[[i]])),
    no.of.edges = length(E(gsO[[i]])), 
    exponent.out = Nalpha, 
    exponent.in = -1,
    loops = FALSE, 
    multiple = FALSE, 
    finite.size.correction = TRUE)
    
     rnd[[i]] <- set.vertex.attribute(rnd[[i]],"name",V(rnd[[i]]),V(gsO[[i]])$name)
     rnd[[i]] <- set.vertex.attribute(rnd[[i]],"GeneName",V(rnd[[i]]),V(gsO[[i]])$GeneName)

    rnd[[i]]      <- findLCC(rnd[[i]])
    names(rnd)[i] <- names(gsO)[i]

}

}



Nn  <- length(names(rnd))
oo  <- matrix(NA,nrow=1,ncol=Nn)
colnames(oo) <- names(rnd)

aa  <- matrix(NA,nrow=1,ncol=Nn)
colnames(aa) <- names(rnd)

tt  <- matrix(NA,nrow=1,ncol=Nn)
colnames(tt) <- names(rnd)

ee  <- matrix(NA,nrow=1,ncol=Nn)
colnames(ee) <- names(rnd)

for( i in 1:Nn ){
 #pFit <- FitDegree( as.vector(igraph::degree(graph=rnd[[i]])), NSIM , NTHREAD, NSEED )
 #oo[1,i] <- pFit[2]
 oo[1,i] <- 0

 val <- assortativity.degree(rnd[[i]],directed=F)
 aa[1,i] <- val

 val <- transitivity(rnd[[i]],type="global")
 tt[1,i] <- val

 val <- graphEntropy(rnd[[i]])
 ee[1,i] <- val

}

outfile <- file("pFIT.csv","w")
write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("assortMixing.csv","w")
write.table(aa, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("globalTrans.csv","w")
write.table(tt, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("graphEntropy.csv","w")
write.table(ee, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

run=0
if( run ){

#--- calculate edge overlaps between orignal and perutated graphs
indx  <- match(names(rnd),names(gsO) )
Nindx <- length(indx)

deg     <- matrix("",ncol=(1+Nindx),nrow=N)
deg[,1] <- ID

bet     <- matrix("",ncol=(1+Nindx),nrow=N)
bet[,1] <- ID

sl      <- matrix("",ncol=(1+Nindx),nrow=N)
sl[,1]  <- ID

cc      <- matrix("",ncol=(1+Nindx),nrow=N)
cc[,1]  <- ID

pr      <- matrix("",ncol=(1+Nindx),nrow=N)
pr[,1]  <- ID

sp      <- matrix("",ncol=(1+3*Nindx),nrow=N)
sp[,1]  <- ID

vars <- vector(length=3)
vars[1] <- "MIN"
vars[2] <- "MEAN"
vars[3] <- "SD"

str <- "ID"
for( i in 1:length(names(rnd)) ){
     str <- c(str,sprintf("%s_%s",vars,names(rnd)[i]))
}

colnames(deg) <- c("ID",names(rnd))
colnames(bet) <- c("ID",names(rnd))
colnames(sl)  <- c("ID",names(rnd))
colnames(cc)  <- c("ID",names(rnd))
colnames(pr)  <- c("ID",names(rnd))
colnames(sp)  <- str
	         
indA <- which(grepl("MIN_PPI",str)==TRUE)
indB <- which(grepl("MEAN_PPI",str)==TRUE)
indC <- which(grepl("SD_PPI",str)==TRUE)

for( i in 1:Nindx ){

    #ggO <- gsO[[indx[i]]]

    gg  <- rnd[[i]]

    ind <- match(V(gg)$name,ID)
    deg[ind,(i+1)] <- as.vector(igraph::degree(graph=gg))
    bet[ind,(i+1)] <- as.character(round(betweenness(gg),3))
    sl[ind,(i+1)]  <- as.character(round(Semilocal(gg),3))
    cc[ind,(i+1)]  <- as.character(round(transitivity(gg,"local"),3))
    pr[ind,(i+1)]  <- as.character(round(as.vector(page.rank(graph=gg,vids=V(gg),directed=F,options=igraph.arpack.default)$vector),6))

    res <- as.matrix(calShorestPaths(gg))
    sp[ind,indA[i]]  <- as.character(res[,1])
    sp[ind,indB[i]]  <- as.character(res[,2])
    sp[ind,indC[i]]  <- as.character(res[,3])

    }

outfile <- file("permuteDEGREE.csv","w")
write.table(deg, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("permuteBET.csv","w")
write.table(bet, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("permuteSL.csv","w")
write.table(sl, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("permuteCC.csv","w")
write.table(cc, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("permutePR.csv","w")
write.table(pr, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

outfile <- file("permuteSP.csv","w")
write.table(sp, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
close(outfile);

}