source('../setUp.R')

##---Check or create output plot dir
if( !file_test("-d","RESULTS") ){
    dir.create("RESULTS")
}


##---Group Disease Terms
Gterms <- list()

#Central Nervous System Disease, including neurdegenerative
CNSdis <- c("DOID:10652","DOID:1826","DOID:9255","DOID:12858","DOID:14330","DOID:2377")

Gterms[[1]]      <- CNSdis
names(Gterms)[1] <- "CNSdis"

#Diseases of Mental Health
MHdis  <- c("DOID:5419","DOID:3312","DOID:12849","DOID:0060041","DOID:1059")

Gterms[[2]]      <- MHdis
names(Gterms)[2] <- "MHdis"

#CardioVascular system disease
CVdis <- c("DOID:10763")

Gterms[[3]]      <- CVdis
names(Gterms)[3] <- "CVdis"

FILES <- "Disease_overlap_sig"

for( s in 1:length(subDIR) ){

    str <- sprintf("RESULTS/%s/%s/%s.csv",subDIR[s],gdaDIR[gdas],FILES)

    if( file.exists(str) && file.info(str)$size!=0 ){

        tt <- read.table(str,sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)
        
        #Filter for Diseases of size >= 10
        tt <- tt[tt[,4] >= 10 & tt[,8] >= 10,]

        #test, Filter identical terms, e.g. AD-AD
        tt <- tt[as.numeric(tt[,9]) != 0,]

        res <- list()
        k=1

        for( i in 1:length(names(Gterms)) ){
            for( j in i:length(names(Gterms)) ){

                if( i == j ){
                    temp <- tt[ifelse(is.na(match(tt[,1],Gterms[[i]])),FALSE,TRUE) & ifelse(is.na(match(tt[,5],Gterms[[j]])),FALSE,TRUE),]
                } else {
                    temp1 <- tt[ifelse(is.na(match(tt[,1],Gterms[[i]])),FALSE,TRUE) & ifelse(is.na(match(tt[,5],Gterms[[j]])),FALSE,TRUE),]
                    temp2 <- tt[ifelse(is.na(match(tt[,1],Gterms[[j]])),FALSE,TRUE) & ifelse(is.na(match(tt[,5],Gterms[[i]])),FALSE,TRUE),]
                    temp <- rbind(temp1,temp2)
                }

                temp <- as.data.frame(temp)
                res[[k]] <- temp
                names(res)[k] <- sprintf("%s_%s",names(Gterms)[i],names(Gterms)[j])
                k=k+1
            }
        }
        
        outfile <- file(sprintf("RESULTS/%s_HDOgrouped_%s.csv",subDIR[s],FILES[1]),"w")
        for( i in 1:length(names(res)) ){
            cat(names(res)[i],file=outfile,"\n")    
            write.table( data.frame(res[[i]]), file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
        }

        close(outfile);

    }
}

for( i in 1:length(names(Gterms)) ){
    Gterms[[i]] <- disl[match(Gterms[[i]],disn)]
}

