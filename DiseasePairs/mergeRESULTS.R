source('../setUp.R')

FILES <- "Disease_overlap_sig.csv"

tb1 <- read.table(sprintf("RESULTS/%s/%s/%s",subDIR[1],gdaDIR[gdas],FILES),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

for( s in 2:length(subDIR) ){
    
    temp <- read.table(sprintf("RESULTS/%s/%s/%s",subDIR[s],gdaDIR[gdas],FILES),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

    tb1 <- cbind(tb1,temp)
    
}

outfile <- file(sprintf("RESULTS/%s_Disease_overlap_sig.csv",gdaDIR[gdas]),"w");
write.table( tb1, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
close(outfile);
