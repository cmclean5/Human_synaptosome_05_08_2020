#--- Merge result files ---#

rm(list=ls())

#library(igraph)

#load corresponding graph
#---Directories
DIR    <- vector(length=3)
DIR[1] <- "../../datasets/"
DIR[2] <- "RESULTS"
DIR[3] <- "OUT"

files    <- vector(length=3)
files[1] <- "pFIT.csv"  
files[2] <- "assortMixing.csv"
files[3] <- "globalTrans.csv" 
files[4] <- "graphEntropy.csv" 

subdirs  = list.files(path=sprintf("../%s/",DIR[2]));
nstudies = length(subdirs);

#runSTEP1=FALSE
runSTEP1=TRUE

if( runSTEP1 ){

 for( f in 1:length(files) ){

    st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[1],files[f]);	
    if( file.exists(st1) && file.info(st1)$size!=0 ){
     tb = read.table(st1,header=T,sep="\t");
    }

     for( s in 2:nstudies ){

      st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[f]);	

      if( file.exists(st1) && file.info(st1)$size!=0 ){
      	  temp <- read.table(st1,header=T,sep="\t");
          tb <- rbind(tb,temp)
      }
     }

     outfile <- file(sprintf("%s/random_%s",DIR[3],files[f]),"w");
     write.table(tb, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
     close(outfile);

     rm(tb)
}

}


#runSTEPS2=TRUE
runSTEPS2=FALSE

if( runSTEPS2 ){

files    <- vector(length=6)
files[1] <- "permuteDEGREE.csv"  
files[2] <- "permuteBET.csv"
files[3] <- "permuteCC.csv"
files[4] <- "permuteSL.csv"
files[5] <- "permuteSP.csv"
files[6] <- "permutePR.csv"


for( f in 1:length(files) ){

st1=""
s=0
while( !(file.exists(st1) && file.info(st1)$size!=0) ){
       s=s+1    
       st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[f]);
  
       if( file.exists(st1) && file.info(st1)$size!=0 ){

       tb = read.table(st1,header=T,sep="\t");
  
       }

}

   CN <- colnames(tb)
   nmat <- length(CN)-1
   N <- length(tb[,1])

   lst <- list()
   for( m in 1:nmat ){
     lst[[m]]           <- matrix(0,ncol=(nstudies+1),nrow=N)
     colnames(lst[[m]]) <- c(CN[1],seq(1,nstudies,1))
     lst[[m]][,1]       <- as.character(tb[,1])
   }


for( s in 1:nstudies ){
  
  st1 = sprintf("../%s/%s/%s",DIR[2],subdirs[s],files[f]);
  
  if( file.exists(st1) && file.info(st1)$size!=0 ){

     temp = read.table(st1,header=T,sep="\t");
     temp = as.data.frame(temp)  
   
     for( m in 1:nmat ){

     #lst[[m]][,(s+1)] <- as.numeric(lst[[m]][,(s+1)]) + as.numeric(temp[,(m+1)])
     lst[[m]][,(s+1)] <- as.numeric(temp[,(m+1)])

     }

  }
}

for( m in 1:nmat ){
 oo <- as.data.frame(lst[[m]])
 outfile <- file(sprintf("%s/random_%s_%s",DIR[3],CN[(m+1)],files[f]),"w");
 write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);
}

rm(lst,tb,temp,oo)

}

}