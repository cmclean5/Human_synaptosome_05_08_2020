rm(list=ls())

library(igraph)


#--- gobal parameters
SEED=1
a=1
ITS=10

#---Script required inputs
#args  <- commandArgs(TRUE);
#SEED  <- as.numeric(args[1]) #random seed no. 
#a     <- as.numeric(args[2]) #algorithm
#ITS   <- as.numeric(args[3]) #no: of iterations

## set random number seed
set.seed(as.numeric(SEED))

## select algorithm
alg = vector(length=2)
alg[1] = "Spectral"
alg[2] = "louvain2"

cat("Running DiseaseLoc.R: \n")
cat("SEED: ", SEED, "\n")
cat("ALG: ", alg[a], "\n")
cat("Permutations: ", ITS, "\n")


gene_com_prob = vector(length=2)
gene_com_prob[1] = "Spectral_Vprobs.csv"
gene_com_prob[2] = "louvain2_Vprobs.csv"


Cn = read.delim(gene_com_prob[a],sep="\t", header=T)
Cn = Cn[,2:length(colnames(Cn))]


#---load corresponding graph which was used to build the consensus matrices from
#subDIR = "full_FilteredPSD_ppi.gml" #Oksana_10_10_2019
subDIR = "PPI_comp.gml" #Oksana_13_11_2020
gg  <- igraph::read.graph(subDIR,format="gml")
ids <- V(gg)$name
n   <- length(ids)
Adj <- igraph::get.adjacency(gg)
#---

## results list
rmVatts = matrix(0,ncol=2, nrow=3)
rmVatts[1,1] = 0
rmVatts[2,1] = 1
rmVatts[3,1] = 1

rmVatts[1,2] = "GeneName"
rmVatts[2,2] = "probPSP"
rmVatts[3,2] = "probPRE"

RES <- list()

for( ri in 2:length(rmVatts[,1]) ){
    coms = table(igraph::get.vertex.attribute(gg, alg[a], V(gg)))
    nr   = length(coms)
    nc   = (3+ITS)
    tmp     = matrix(0, ncol=nc, nrow=nr)
    tmp[,1] = names(coms)
    tmp[,2] = as.vector(coms)
    tmp[,3] = log(sum(as.vector(coms))/as.numeric(coms))
    colnames(tmp) = c("C","Cn","Cscale",seq(1,ITS,1))
    RES[[(ri-1)]] = tmp
    names(RES)[(ri-1)] = rmVatts[ri,2]
}

#---

#--- path to study sets
#st1 = "ALLPapAllGene.csv" #Oksana_10_10_2019
#filein = read.table(st1, header=FALSE, sep="\t");

st1 = "AllGenesAllPapers.csv" #Oksana_13_11_2020
filein = read.table(st1, header=TRUE, sep="\t");

# filein$V2 == study compartment/region
# filein$V4 == Human Enretz.ID
# filein$V9 == study name

#remove any NA fields in Human.Entrez.ID column
filein = filein[!is.na(filein[,4]),]


## probability of seeing a region across all studies

regions <- unique(as.vector(filein[,2]))

## remove "Synaptosome"
regions <- regions[-which(regions=="Synaptosome")]

reg_prior = matrix(0,ncol=2, nrow=length(regions))
reg_prior[,1] = regions

R = length(unique(filein[,9]))
R = R - length(unique(filein[,9][filein[,2] == "Synaptosome"]))

for( i in 1:length(regions) ){

    r = length(unique(filein[,9][filein[,2] == reg_prior[i,1]]))
    #r = length(filein[,2][which(filein[,2] == reg_prior[i,1])])

    reg_prior[i,2] = R/r
        
}

reg_prior[,2] = as.numeric(reg_prior[,2])/sum(as.numeric(reg_prior[,2]))

#-------------------------

## probability of seeing a gene across all studies

SIDS = unique(filein[,4])
N    = length(SIDS)
SETS = list()


## create study sets
for( r in 1:length(regions) ){

    studies = unique(filein[,9][filein[,2]==regions[r]])
    M       = length(studies)
    
    mm = matrix(0, ncol=(1+M), nrow=N)
    colnames(mm) <- c("Human.Entrez.ID",studies)
    
    mm[,1] = SIDS

    for( s in 1:M ){
    
        indx = match(filein[,9],studies[s])
        indx = ifelse(!is.na(indx),TRUE,FALSE)

        sIDs  = filein[indx,4]
        indx2 = match(sIDs,mm[,1])
        indx2 = indx2[!is.na(indx2)]

        mm[indx2,(1+s)] = 1

    }

    SETS[[r]] = mm
    names(SETS)[r] = regions[r]
    
}




## weighting each study by log( (# network nodes)/sum(1's in study) )
SETSwe = list()

for( s in 1:length(SETS) ){

    sn = length(colnames(SETS[[s]]))

    sN = length(SETS[[s]][,1])

    temp = matrix(0, ncol=sn, nrow=1)
    colnames(temp) <- colnames(SETS[[s]])
       
    for( i in 2:sn ){

        ones = sum(as.numeric(SETS[[s]][,i]))

        temp[1,i] = log( sN/ones )
        
        
    }

    SETSwe[[s]] = temp
    names(SETSwe)[s] = names(SETS)[s]
    
}

##start timing
ptm <- proc.time()

cat("Running...\n")
for( its in 1:ITS ){

      cat("Permutation No: ", its, "\n")	

    ## randomise study sets
    for( r in 1:length(SETS) ){

        cn = length(colnames(SETS[[r]]))

        for( s in 2:cn ){
            SETS[[r]][,s] = sample(SETS[[r]][,s])        
        }
        
    }



    ## gene prob, given region
    gene_reg = matrix(ncol=(1+length(regions)), nrow=N)
    gene_reg[,1] = SIDS

    for( s in 1:length(SETS) ){

        cn = length(colnames(SETS[[s]]))
        WE = 0
    
        for( i in 1:N ){

            we = as.numeric(SETS[[s]][i,2:cn]) * as.numeric(reg_prior[s,2])
            we = sum(we * as.numeric(SETSwe[[s]][1,2:cn]))+1 ## +1 is the Laplacean prior
	    ## we = sum(as.numeric(SETS[[s]][i,2:cn]))+1 ## +1 is the Laplacean prior
            ## we = 1
            WE = WE + we 

            gene_reg[i,(s+1)] = we

        }

        gene_reg[,(s+1)] = as.numeric(gene_reg[,(s+1)])/WE
    
    }

#-----


    ## gene prob in network
    gene_ntwrk = gene_reg[match(ids,gene_reg[,1]),]
    gene_ntwrk[,2] = gene_ntwrk[,2]/sum(gene_ntwrk[,2])
    gene_ntwrk[,3] = gene_ntwrk[,3]/sum(gene_ntwrk[,3])

    #for( r in 2:length(colnames(gene_ntwrk)) ){
    #    gene_ntwrk[,r] = gene_ntwrk[,r]/sum(gene_ntwrk[,r])
    #}
    

   ##PSP    
   #if( is.null(igraph::get.vertex.attribute(gg,as.character(rmVatts[2,2]))) ){
   #   gg = igraph::set.vertex.attribute(gg,as.character(rmVatts[2,2]),V(gg),gene_ntwrk[,2])
   #  }	else {
   #   gg = igraph::remove.vertex.attribute(gg,as.character(rmVatts[2,2]))
   #   gg = igraph::set.vertex.attribute(gg,as.character(rmVatts[2,2]),V(gg),gene_ntwrk[,2])
   # }

   ##PRE    
   #if( is.null(igraph::get.vertex.attribute(gg,as.character(rmVatts[3,2]))) ){
   #   gg = igraph::set.vertex.attribute(gg,as.character(rmVatts[3,2]),V(gg),gene_ntwrk[,3])
   #  }	else {
   #   gg = igraph::remove.vertex.attribute(gg,as.character(rmVatts[3,2]))
   #   gg = igraph::set.vertex.attribute(gg,as.character(rmVatts[3,2]),V(gg),gene_ntwrk[,3])
   # }


    for( ri in 2:length(rmVatts[,1]) ){

        X1 = Adj %*% as.vector(gene_ntwrk[,ri]) #igraph::get.vertex.attribute(gg,rmVatts[ri,2], V(gg))
   
        Creg = t(Cn) %*% X1
   
        RES[[(ri-1)]][,(3+its)] = as.numeric(RES[[(ri-1)]][,3]) * as.numeric(Creg)
    
        if( !is.null(X1) )  { rm(X1)   }
        if( !is.null(Creg) ){ rm(Creg) }
    }

    ## proportions test
    ##we = as.numeric(RES[[1]][,(3+its)]) + as.numeric(RES[[2]][,(3+its)])
    ##RES[[1]][,(3+its)] = as.numeric(RES[[1]][,(3+its)])/we
    ##RES[[2]][,(3+its)] = as.numeric(RES[[2]][,(3+its)])/we

}

for( r in 1:length(RES) ){
   write.table(RES[[r]], sprintf("results_%s.csv",names(RES)[r]),sep="\t", row.names=F, col.names=T, quote=F)
}


pet <- proc.time() - ptm
cat("Finished! ", sprintf("time = %.3f \n", pet[[1]]), "\n")
