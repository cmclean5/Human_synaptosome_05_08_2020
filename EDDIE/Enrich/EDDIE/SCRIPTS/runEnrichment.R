rm(list=ls())

options(stringsAsFactors=F)

getCI <- function(x, indx=1){
    strsplit(x,",")[[1]][indx]
}

buildIndx <- function(x, MAX=1){
	  seq(x,(x+MAX),1)
}

#---Script required inputs
args   <- commandArgs(TRUE);
RUNS   <- as.numeric(args[1]) #No: of randomisation runs
SEED   <- as.numeric(args[2]) #offset for random number seed
FILTER <- as.numeric(args[3]) #Should we filter on specific (C,GO:) pairs  

cat("R> Random permutations: ", RUNS,"\n")
cat("R> Random number seed : ", SEED, "\n")

#---READ IN PARAMETERS
#--- Pre-load Graph of interest, stored in file 'graphs.csv'
pramFILES <- "parameterFiles"

#---READ IN GRAPH
Graph   <- read.delim(sprintf("%s/graphs.csv",pramFILES),header=F,sep="\t",quote="")
subDIR  <- as.vector(Graph[which(as.vector(Graph[,1]) == 1)[1],2])

#---READ IN ALGORITHM 
algs <- read.delim(sprintf("%s/clusteringAlg.csv",pramFILES),header=F,sep="\t",quote="")
ALGS <- as.vector(algs[which(as.vector(algs[,1]) == 1),3])

#---READ IN ANNOTATION
anno <- read.delim(sprintf("%s/annotations.csv",pramFILES),header=F,sep="\t",quote="")
ANNO <- as.vector(anno[which(as.vector(anno[,1]) == 1),3])

#---READ IN FILTER LIST
cands <- read.delim(sprintf("%s/candidate_GOBP.csv",pramFILES),header=T,sep="\t",quote="")

if( FILTER == 1 ){
    cat("R> Filter on enriched community / Anno pairs.\n")
}

searchTERM = "GO"

RES <- list()

r=1
for( r in 1:RUNS ){

cat("R> Run = ", r, "...")

system(sprintf("./submitClustEnrch.sh %d", SEED))

k=1
for( a in 1:length(ALGS) ){

     if( FILTER == 1 ){
       candX = cands[cands[,1]==ALGS[a],]
     }

    str <- sprintf("RESULTS/%s/%s/permute_p_values_%s.csv",subDIR[1],ALGS[a],ANNO[1])

    cat("R> ",str,"\n")

     if( file.exists(str) && file.info(str)$size!=0 ){
    
         cat("\n")
         cat("R> running algorithm: ", ALGS[a], "...")

         #--- load functional enrichment file
         tt <- read.delim(str,sep="\t", header=T)

	 if( FILTER == 1 ){
           candA = unique(candX[,2]) 
           candC = unique(candX[,3])
           if( length(candA) > 0 && length(candC) > 0 ){
	      temp = tt[match(candC,tt[,1]),]
	      rm(tt)
              cn   = colnames(temp)
              fn   = grep(searchTERM,cn)
              MAXc = abs(fn[1] - fn[2]) -1
              SELc = match(unique(candA),cn)
              SELc = SELc[order(SELc)]
              indx = as.vector(sapply(SELc,buildIndx, MAX=MAXc))
              temp = temp[,c(1,2,indx)]
              tt   = temp
	      rm(temp,cn,fn,MAXc,SELc,indx,candA,candC)
	      }
         }

         fn <- colnames(tt)[grepl(searchTERM,colnames(tt))]
         FN <- length(fn)

         CN <- length(tt[,1])

         DF <- matrix("",nrow=(FN*CN), ncol=7)
         colnames(DF) <- c("Fn","C","Mu","OR","CIl","Pv","Ap")


         DF[,1] = rep(fn,CN)
         
         temp1 <- c()
         temp2 <- c()
         temp3 <- c()
         temp4 <- c()
         temp5 <- c()
         temp6 <- c()
         temp7 <- c()

         for( i in 1:length(tt[,1]) ){
             temp1 <- c(temp1, rep(tt[i,1],FN))
         
	     ov = as.vector(unlist(tt[i,grepl("actual",colnames(tt))]))
             temp2 <- c(temp2, ov)

	     or = as.vector(unlist(tt[i,grepl("OR",colnames(tt))]))
             temp3 <- c(temp3, or)

             ci = as.vector(unlist(tt[i,grepl("CI",colnames(tt))]))
             ci = gsub("\\[","",ci)
             ci = gsub("\\]","",ci)
             temp4 <- c(temp4,as.vector(sapply(ci, getCI, indx=1)))
             temp5 <- c(temp5,as.vector(sapply(ci, getCI, indx=2)))

	     pv = as.vector(unlist(tt[i,grepl("p.value", colnames(tt)) &
                                      !grepl("X.p_value",colnames(tt)) &
                                      !grepl("ALT",colnames(tt))]))
             temp6 <- c(temp6, pv)
             
	     pa = as.vector(unlist(tt[i,grepl("adjusted",colnames(tt)) &
                                      !grepl("ALT",colnames(tt))]))
             temp7 <- c(temp7, pa)

         }
   
         DF[,2] = temp1
 	 DF[,3] = temp2
         DF[,4] = temp3
         DF[,5] = temp4
         DF[,6] = temp6
         DF[,7] = temp7

	 rm(temp1, temp2, temp3, temp4, temp5, temp6, temp7)

         DF = as.data.frame(DF)

	 if( r == 1 ){
            RES[[k]] = DF
            names(RES)[k] = ALGS[a]
         } else {
             indx = which(ALGS[a] == names(RES))
             if( length(indx) != 0 ){
             	     RES[[indx[1]]] = cbind(RES[[indx[1]]],DF[,c(3,4,5,6,7)])
	       }
         }      

         k=k+1
         rm(DF)
     
         cat("...done.\n")
     }

   }

r=r+1

cat("..done.\n")

}

##save(RES, file="enrichment.RData", compress=TRUE)
##saveRDS(RES, file="enrichment.rds", compress=TRUE)

cat("Compress files...")

for( s in 1:length(RES) ){

  oo      <- as.data.frame(RES[[s]])
  CN      <- colnames(oo)
  CN      <- gsub("Mu.*" ,"Mu" ,CN)
  CN      <- gsub("OR.*" ,"OR" ,CN)
  CN      <- gsub("CIl.*","CIl",CN)	
  CN      <- gsub("Pv.*" ,"Pv" ,CN)
  CN      <- gsub("Ap.*" ,"Ap" ,CN)

  colnames(oo) <- CN

  write.table(oo,gzfile(sprintf("random_studies_%s.csv.gz",names(RES)[s]))) # write a compressed file

}

cat("done.\n")
rm(RES)
    
