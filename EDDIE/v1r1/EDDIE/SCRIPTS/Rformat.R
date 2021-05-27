#generate consensus file, using Spectral Mod. clustering results

#read subsampling from Rclustering.R script
tt1 <- read.table("consensusout.txt",sep="\t",skip=1)

#read Spetral clustering results
memb <- read.table("OUT/consensusout.txt",sep=" ", skip=1)

#store results
for( i in 1:length(tt1[,1]) ){
     indx <- which(tt1[i,2]==memb[,2])
     if( length(indx) != 0 ){
       tt1[i,3] <- memb[indx[1],3]
     }

}


#write consensus file
cc <- as.data.frame(tt1)
outfile <- file("consensusout.txt","w")
cat("#consensus",file=outfile,"\n")
write.table( cc, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
close(outfile);
