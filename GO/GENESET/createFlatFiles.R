##
source('../../setUp.R')

##---Directories
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
#OUT[2] <- DIRS[grepl("Backgrounds",DIRS)]
OUT[2] <- "GOClips"#DIRS[grepl("GOClips",DIRS)]
OUT[3] <- "OUT"#DIRS[grepl("OUT",DIRS)]


##--- GO Clip
CLIPs <- vector(length=3)
CLIPs[1] <- "synapse"
CLIPs[2] <- "NIGO"
CLIPs[3] <- "myGOterms"

#set clip
#c=1

Path  <- OUT[2]
files <- list.files(path=Path);

ClippedGOterms <- c()

for( c in 1:length(CLIPs) ){

    file <- files[grepl(CLIPs[c],files)]

    if( length(file) != 0 ){
    
    str <- sprintf("%s/%s",Path, file)

    if( file.exists(str) && file.info(str)$size!=0 ){
    
        #fileName       <- strsplit(file[1],".csv")
        GOterms <- read.table(sprintf("%s/%s",Path,file[1]),sep="\t",quote="")
        GOterms <- as.vector(GOterms[[1]])
        ClippedGOterms <- c(ClippedGOterms, GOterms)

    }
    }
}

ClippedGOterms <- unique(ClippedGOterms)

##------------------------

PATH <- sprintf("%s/synaptic_genes_hum_5016_05082020/",OUT[3])

files <- list.files(path=PATH)

##--- cut-off
##    We'll only calculate the enrichment on GO terms, >= MIN & <= MAX, number of annotations
thresMIN <- 10
thresMAX <- 50000

foundterms <- vector()

for( i in 1:length(GOonto) ){

    strOO <- sprintf("flatfile.go.%s.csv",GOonto[i])

    file <- files[grepl(sprintf("flatfile_Human__%s",GOonto[i]),files)]

    str  <- sprintf("%s/%s",PATH,file)
    
    if( file.exists(str) && file.info(str)$size!=0 ){

        tt <- read.table( str, sep="\t", header=F, quote="")

        indx <- match( tt[,1], ClippedGOterms )

        oo <- tt[ifelse(is.na(indx),FALSE,TRUE),]

        GOterms <- table(oo[,1])
        foundterms <- c(GOterms,foundterms)
        indA    <- match(oo[,1],names(GOterms)[as.vector(GOterms) >= thresMIN & as.vector(GOterms) <= thresMAX])
        oo      <- oo[ifelse(is.na(indA),FALSE,TRUE),]
        
        write.table(oo,strOO,sep="\t",col.names=F,row.names=F,quote=F)
        
    }

}
