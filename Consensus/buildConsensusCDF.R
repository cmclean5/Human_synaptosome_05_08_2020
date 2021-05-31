source('../setUp.R')

run=0

if( run ){

rm(list=ls())

#set default R options 
options(stringsAsFactors=F)

#-- unresoved problem in setUp.R, one of the library loaded causes a segmentation fault when loading consensus.txt files
#-- so will just load the single library we need to build/store the CDF results.
#library("data.table")

#Set the absolute path to this working directory
version="17_05_2019"
mainDir <- sprintf("/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/Human_synaptosome_%s",version)

#Get Path to all top-level directories in folder
DIRS <- list.dirs(mainDir,recursive=F)

#---dataDIR, where to read PPI files
dataDIR <- DIRS[grepl("datasets",DIRS)]

#---subDIR for loading and storing PPI Graphs
subDIR <- list.files(path=dataDIR)
subDIR <- subDIR[grep(".csv",subDIR)]
subDIR <- unlist(strsplit(subDIR,".csv"))

#---Location for randomisation files
rndDIR    <- vector(length=1)
rndDIR[1] <- sprintf("/disk/scratch/WORK/DATA/Human_synaptosome_%s",version)

#---Declare all clustering algorithms
#ALGS    <- vector(length=21)
#ALGS[1]  <- "fc"
#ALGS[2]  <- "sgG1"
#ALGS[3]  <- "sgG2"
#ALGS[4]  <- "sgG5"
#ALGS[5]  <- "Spectral"
#ALGS[6]  <- "louvain"
#ALGS[7]  <- "infomap"
#ALGS[8]  <- "lec"
#ALGS[9]  <- "wt"
#ALGS[10] <- "SVI"
#ALGS[11] <- "SPICi"
#ALGS[12] <- "Geodesic"
#ALGS[13] <- "CONSENSUS"
#ALGS[14] <- "louvain2"
#ALGS[15] <- "Spectral1per"
#ALGS[16] <- "Spectral25per"
#ALGS[17] <- "Spectral5per"
#ALGS[18] <- "fc2"
#ALGS[19] <- "lec2"
#ALGS[20] <- "wt2"
#ALGS[21] <- "sgG12"

#--- Set 
#--- Pre-load Graph of interest, stored in file 'graphs.csv'
pramFILES <- DIRS[grepl("parameterFiles",DIRS)]
Graph <- read.table(sprintf("%s/graphs.csv",pramFILES),header=F,sep="\t",quote="")
S     <- as.vector(Graph[which(as.vector(Graph[,1]) == 1)[1],2])
S     <- match(S,subDIR)

}
    
#Check or create out dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}


    
#---Clustering algorithms used
#set  <- c("fc","fc2","lec","lec2","louvain","louvain2","infomap","sgG1","sgG12","Spectral","wt","wt2")
set  <- c("fc","lec","louvain","infomap","sgG1","Spectral","wt")
algs <- ALGS[match(set,ALGS)]
#algs  <- ALGS[c(1,2,5,6,7,8,9,14,18,19,20,21)]
#algs <- ALGS[-c(3,4,10,11,12,13,15,16,17,22,23,24,25,26,27,28,29)]


#set CDF parameters
steps  = 100;
stepsi = 1.0/steps;
CDF <- vector(length=steps);

X   <- vector(length=steps)
for(x in 1:steps ){
    Xi=stepsi*x;
    X[x] <- Xi;
}      

#output results
oo <- matrix("",ncol=(1+length(algs)),nrow=steps)
colnames(oo) <- c("X",algs)
oo[,1]       <- X

#run over each algorithm
for( a in 1:length(algs) ){

    #print the algorithm we're running over
    cat("running over ", algs[a], "\n");

    st1 = sprintf("%s/%s/%s/consensusmatrix.txt.gz",rndDIR[1],subDIR[S],algs[a])
    
    if( file.exists(st1) && file.info(st1)$size!=0 ){
    
    #Read in consensus matrix
    CM = read.table(gzfile(st1) , header=FALSE, sep=","); 
    
    #reset CDF vector
    CDF <- rep(0,length(CDF))

    #build consensus matrix CDF
    Fn  <- ecdf(CM[lower.tri(CM)])

    CDF <- Fn(X[1:steps])

    #store CDF result
    oo[,(1+a)] <- CDF

    #output file
    outfile <- file(sprintf("%s/%s_consensusCDFs.csv",subDIR[S],subDIR[S]),"w")
    write.table(oo, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
    close(outfile);

    rm(CM,Fn);
    #gc();

    }
}


    
