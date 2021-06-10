rm(list=ls())

options(stringsAsFactors=F)

require("homologene")

DB = homologeneData2

tax = homologene::taxData

ids <- read.delim("/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/Human_synaptosome_17_05_2019/GeneSets/Backgrounds/synaptic_genes_hum_7676_04122019.csv", sep="\t", header=F)
ids <- unique(ids[[1]])

oo <- matrix(0, ncol=(1+length(tax[,1])), nrow=length(ids))
colnames(oo) <- c("Human_Entrez_ID",tax[,1])
oo[,1] <- ids

for( i in 1:length(tax[,1]) ){

    temp <- homologene(genes=ids,inTax=9606, outTax=tax[i,1], db=DB)
    indx <- match(oo[,1],temp[,3])
    indx <- ifelse(is.na(indx),0,1)
    oo[,(1+i)] <- indx
    
}

outfile <- file("homologene.mapping.csv", "w")
write.table( oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);

##---
## Chen et al. Dissecting the Human Protein-Protein Interaction Network via Phylogenetic Decomposition
## Scientific Reports, 4: 7153 (2014).
##---
groups <- list()

groups[[1]]      <- c(33169,28985,318829,5141,4896,4932)
names(groups)[1] <- "G1" # FUNGUS

groups[[2]]      <- c(6239,7227,7165)
names(groups)[2] <- "G2" # FLIES

groups[[3]]      <- c(7955)
names(groups)[3] <- "G3" # ZEBRAFISH

groups[[4]]      <- c(9031)
names(groups)[4] <- "G4" # CHICKEN

groups[[5]]      <- c(9913,9615,10116,10090)
names(groups)[5] <- "G5" # RAT/MOUSE

groups[[6]]      <- c(9606,9598,9544)
names(groups)[6] <- "G6" # HUMAN/PRIMATES

ss <- rep("",length(colnames(oo)))

for( i in 1:length(groups) ){

    tmp <- match(colnames(oo),groups[[i]])
    tmp <- ifelse(is.na(tmp),"",names(groups)[i])

    indx <- tmp != ""
    ss[indx] = names(groups)[i]

}

oo2           <- oo
colnames(oo2) <- ss

tally     <- matrix(0,ncol=(1+length(groups)), nrow=length(ids))
tally[,1] <- ids
colnames(tally) <- c("Human_Entrez_ID",names(groups))

for( i in 1:length(tally[,1])){

    tmp <- oo2[i,]
    for( g in 1:length(groups) ){
        tally[i,(1+g)] <- sum(tmp[which(names(tmp) == names(groups)[g])])
    }    
   
}

label <- matrix("",ncol=2,nrow=length(tally[,1]))
label[,1] <- tally[,1]

for( i in 1:length(label[,1]) ){

    ind = which(tally[i,] != 0)
    ind = ind[-1]
    
    if( length(ind) != 0 ){
        #label[i,2] <- paste(colnames(tally)[ind],collapse=";")
        label[i,2] <- colnames(tally)[ind[1]]
    }
    
}

colnames(label) <- c("Human_Entrez_ID","Evo_Age_Group")
outfile <- file("synaptic_proteins_age_group.csv", "w")
write.table( label, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);
