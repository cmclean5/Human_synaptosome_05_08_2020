##
source('../../setUp.R')

##---Directories
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- "Backgrounds"
OUT[3] <- "GOClips"

##SET output dir
outDIR <- "OUT"

if( !file_test("-d",outDIR) ){
    dir.create(outDIR)
}

TAX        <- TAXid[2]
Onto       <- GOonto[3]
STUDYtitle <- ""

if(TAX=="10090"){
    TAXmap=org.Mm.eg.db
    TAXmapStr="org.Mm.eg.db"    
    TAXtitle="Mouse"
}

if(TAX=="9606"){
    TAXmap=org.Hs.eg.db
    TAXmapStr="org.Hs.eg.db"    
    TAXtitle="Human"
 }

if(TAX=="7227"){
    TAXmap=org.Dm.eg.db
    TAXmapStr="org.Dm.eg.db"
    TAXtitle="Fly"
 }

geneID2TERM <- revmap(topGO::annFUN.org(feasibleGenes=NULL,whichOnto=Onto, mapping = TAXmapStr, ID = "Entrez"))
TERM2geneID <- revmap(geneID2TERM)

geneNames=names(geneID2TERM)

##--- Change Background set to synaptic gene list of 6899 genes.
##bgfile <- "synaptic_genes_hum_5387_11042019"
##bg <- read.table(sprintf("%s/%s/%s.csv",OUT[1],OUT[2],bgfile),sep="\t",header=T)[[1]]
##geneID2TERM <- geneID2TERM[names(geneID2TERM) %in% bg]
##TERM2geneID <- revmap(geneID2TERM)
##geneNames = names(geneID2TERM)

##cat("Running over synaptic Background: ", bgfile, "\n")
##------------------------

##--- GO Clips
CLIPs <- vector(length=2)
CLIPs[1] <- "synapse"
CLIPs[2] <- "NIGO"

##set clip
c=2

#runGOClip = TRUE
runGOClip = FALSE
CLIPtitle = ""

if( runGOClip ){
    CLIPtitle = CLIPs[c]
    
    Path  <- OUT[3]
    files <- list.files(path=Path);
    files <- files[grepl(CLIPtitle,files)]

    fileName       <- strsplit(files[1],".csv")[[1]]
    ClippedGOterms <- read.table(sprintf("%s/%s",Path,files[1]),sep="\t",quote="")
    ClippedGOterms <- as.vector(ClippedGOterms[[1]])
    
    TERM2geneID.filtered<-TERM2geneID[names(TERM2geneID) %in% ClippedGOterms]
    geneID2TERM<-revmap(TERM2geneID.filtered)
    
    geneNames=names(geneID2TERM)

    cat("Running over Clip: ", fileName, "\n")
    
}
##------------------------


##READ-IN THE DATASET FROM VENN/
Path  <- sprintf("%s/%s",OUT[1],OUT[2])
files <- list.files(path=Path);

files <- files[grepl(".csv",files)]
#files <- files[grepl("PPI_PSP.csv",files)]
#files <- files[!grepl(OUT[2],files)]

if( length(files) > 0 ){

    for( f in 1:length(files) ){
    
    cat("Running over file: ", files[f], "\n")

    fileName <- strsplit(files[f],".csv")
    
    str = sprintf("%s/%s",Path,files[f]);

    if( file.exists(str) && file.info(str)$size!=0 ){
    
        dd <- read.table(sprintf("%s/%s",Path,files[f]),sep="\t",header=T)[[1]]

        rm(geneList,GOdata,resultFisher,topnodes,allRes,myList)

        ##---Get Entrez and GeneName columns from input file
        myList <- as.vector(dd)
        myList <- unique(myList)
        myList <- na.omit(myList)
        
        gn <- mapIds(TAXmap,as.character(myList),column="SYMBOL",keytype="ENTREZID")

        geneList <- factor(as.integer(geneNames %in% myList))
        names(geneList) <- geneNames

        GOdata <- new("topGOdata", ontology=Onto, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2TERM)

        ##--- Run enrichment using topGO
        resultFisher        <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        resultFisher.elim   <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
        
        topnodes = resultFisher@geneData[[4]]

        if( resultFisher.elim@geneData[[4]] < topnodes ){
            topnodes = resultFisher.elim@geneData[[4]]
        }

        allRes <- GenTable(GOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topnodes )

        hits=lapply(allRes$GO.ID,function(x) intersect(sigGenes(GOdata),genesInTerm(GOdata, c(x))[[1]]))
        allRes$hits=hits

        ##---Get Gene.Names
        GeneList <- list()
        for( i in 1:length(allRes$hits) ){
          GeneList[[i]] <- as.vector(unique(as.vector(gn)[match(allRes$hits[[i]],names(gn))]))
    }


        ##---output data
        vv <- c(colnames(allRes),c("Bonferroni for classicFisher","Bonferroni for elimFisher","BY for classicFisher", "BY for elimFisher","GeneNames"));
        rr <- matrix("",  length(allRes[[1]]), length(vv) );
        colnames(rr) <- vv;
               
        for( i in 1:length(allRes[[1]])){
            for( j in 1:(length(allRes[1,]))){
                if( j == length(allRes[1,]) ){
                    rr[i,j] <- as.character( length(allRes[[j]][i][[1]]) );
                } else {
                    rr[i,j] <- as.character(allRes[[j]][i]);
                }
            }
        }
        ##---
        
    ##---Bonferroni correction
    Nlevels = length(nodes(GOdata@graph));

    INDCF   = which(colnames(rr)=="classicFisher")
    INDEF   = which(colnames(rr)=="elimFisher")        

    rr[,INDCF[1]] <- as.character(gsub("< ","",rr[,INDCF[1]]))
    rr[,INDEF[1]] <- as.character(gsub("< ","",rr[,INDEF[1]]))
            
    INDCFbc = which(colnames(rr)=="Bonferroni for classicFisher")
    INDEFbc = which(colnames(rr)=="Bonferroni for elimFisher")        
            
    for( i in 1:length(allRes[[1]])){
        rr[i,INDCFbc[1]] <- "";
        rr[i,INDEFbc[1]] <- "";        
                                        
        for( l in 1:length(alpha) ){

            ##classicFisher
            if( as.numeric(rr[i,INDCF[1]]) < (alpha[l]/Nlevels) ){
                rr[i,INDCFbc[1]] <- stars[l];
            }
       
        
            ##elimFisher
            if( as.numeric(rr[i,INDEF[1]]) < (alpha[l]/Nlevels) ){
                rr[i,INDEFbc[1]] <- stars[l];
            }
        
        }
    }

    ##---Benjamini & Yekutieli
    ##classic Fisher
    indBHF1 <- INDCFbc+2
    rr[,indBHF1[1]] <- p.adjust(as.numeric(rr[,INDCF[1]]),method="BY",n=Nlevels)
   
    ##elim Fisher
    indBHF2 <- INDEFbc+2
    rr[,indBHF2[1]] <- p.adjust(as.numeric(rr[,INDEF[1]]),method="BY",n=Nlevels)
          

    ##---Insert gene.names which hit the ontology term
    for( i in 1:length(allRes[[1]]) ){
        rr[i,(indBHF2[1]+1)] <- paste(GeneList[[i]],collapse=" ")
        
    }

       ##SET study output dir
        subDIR <- sprintf("%s/%s",outDIR,fileName)

        if( !file_test("-d",subDIR) ){
            dir.create(subDIR)
        }

    ##---enrichment file        
    outfile <- file(sprintf("%s/output_%s_%s_%s.csv",subDIR,TAXtitle,CLIPtitle,Onto), "w");
    write.table( rr, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
    close(outfile)

    ##---format enrichment file
    zz <- matrix(NA,  sum(allRes[[4]]), 3 );
    t=1;

    for(i in 1:length(hits) ){
    
        for(j in 1:(length(as.numeric(hits[[i]]))+1) ){

            if( !is.na(hits[[i]][j]) ){    
                zz[t,1] <- as.character(allRes[[1]][i]);
                zz[t,2] <- as.character(gsub(" ","_",allRes[[2]][i]));
                zz[t,3] <- as.character(hits[[i]][j]);
                t=t+1;
            }
        }     
    }

    ##---output flatfile of enrichment file above 
    outfile <- file(sprintf("%s/flatfile_%s_%s_%s.csv",subDIR,TAXtitle,CLIPtitle,Onto), "w");
    write.table( zz, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
    close(outfile);

        
     cat(" ****** ",length(unique(zz[,3]))," out of ", length(myList), " Entrez IDs found ****** \n " )

    outfile <- file(sprintf("%s/info_%s_%s_%s.csv",subDIR,TAXtitle,CLIPtitle,Onto), "w");

    sink(append=T,file=outfile)
    print(GOdata)
    sink();

    cat(" ****** NO: LEVELS FOR BONFERRONI CORRECTION ", Nlevels, file=outfile, append=T, sep="" );      

        close(outfile);

    }
    }
}



cat("\n")
cat("Done.\n")

##------------------------------------



