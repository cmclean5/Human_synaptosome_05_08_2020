#source("http://bioconductor.org/biocLite.R")
cat("##############################################################################\n")
cat("# The script runs with HDO.\n")
cat("# \n")
cat("##############################################################################\n")

#turn-off debugging
debuggingState(on=FALSE)

library(DBI)
library(topOnto) 

cat("Loading HDO objects from db...\n")
topOnto::initONT('HDO')                #HDO Ontology [Used]
#topOnto::initONT('NIGO')              #GO Ontology  [Used]

#topOnto::initONT('HDOCORENIGO')       #HDO Ontology
#topOnto::initONT('HDOSHAREGENENIGO')  #HDO Ontology
#topOnto::initONT('HDO150')            #HDO Ontology [Used]

#topOnto::initONT('MPO')                #MPO Ontology [Used]

#--- Mapping Human Entrez to Gene.Symbol
require('org.Hs.eg.db')
entrez2symbol<-revmap(as.list(org.Hs.egSYMBOL2EG))

#path to root of output dir.
mainDir <-  "/afs/inf.ed.ac.uk/user/c/cmclean5/Downloads/topOnto-master/OUT/Human_synaptosome_05_08_20/"


#path to gene-disease annotation files
annoDir <- "/afs/inf.ed.ac.uk/user/c/cmclean5/Downloads/topOnto-master/inst/extdata/annotation"

#path to gene list
geneDir <- "/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/Human_synaptosome_05_08_2020/GeneSets/"


#annotation file 
ss <- vector(length=1);
#annotation file related to HDO
ss[1] <- "human_gene2HDO"

#---Human_synaptosome pre and post 17/05/19
ff    <- vector(length=1)
#ff[1] <- "synaptic_genes_hum_4679_17052019"
ff[1] <- "synaptic_genes_hum_5016_05082020"


for( f in 1:length(ff) ){

   rm(myInterestingGenes)  

   #---For Human Synapotsome
   myInterestingGenes=read.table(file=sprintf("%s/%s.csv",geneDir,ff[f]),sep="\t",header=F)[[1]]
    
   myInterestingGenes=unique(myInterestingGenes)
   myInterestingGenes=na.omit(myInterestingGenes)

   myInterestingGenes <- myInterestingGenes[myInterestingGenes != ""]
       
    dir.create( file.path( mainDir, ff[f]) );

   for( s in 1:length(ss) ){

     rm(geneList,rank,resultKS,geneID2GO,hits,resultFisher.elim,GOdata,resultKS.elim,resultFisher,allRes,geneNames,outfile,Nlevels)

     dir.create( file.path( sprintf("%s/%s",mainDir,ff[f]),ss[s]) )
     
     ## HDO annotation file 
     geneID2GO    <- readMappings(file=sprintf("%s/%s",annoDir,ss[s]))
   
       
   TERM2geneID  <- revmap(geneID2GO)
       
     #my genes
     geneNames=names(geneID2GO)
     geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
     names(geneList) <- geneNames



     GOdata <- new("topONTdata", ontology = "HDO", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
     resultFisher    <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
     resultKS        <- runTest(GOdata, algorithm = "classic", statistic = "ks")
     resultFisher.elim   <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
     resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

     topnodes = resultFisher@geneData[[4]]

    if( resultFisher.elim@geneData[[4]] < topnodes ){
      topnodes = resultFisher.elim@geneData[[4]]
    }

    if( resultKS@geneData[[4]] < topnodes ){
      topnodes = resultKS@geneData[[4]]
    }

    if( resultKS.elim@geneData[[4]] < topnodes ){
      topnodes = resultKS.elim@geneData[[4]]
    }

   
    allRes <- GenTable(GOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = topnodes, show.gene.sig=T, use.symbol=T,entrez2symbol=entrez2symbol)


   
    Colnames <- c(names(allRes),"Rank","Bonferroni for classicFisher","Bonferroni for elimFisher","Benjamini & Yekutieli (BY) for classicFisher", "Benjamini & Yekutieli (BY) for elimFisher")
       
    #---output data
    rr <- matrix(NA, nrow=length(allRes[,1]), length(Colnames) );
    colnames(rr) <- Colnames;

    methodNames <- c("classicFisher","elimFisher","classicKS","elimKS")
    indF <- vector(length=length(methodNames))
    for( i in 1:length(indF) ){
        indF[i] = -1
        ind <- which(colnames(rr)==methodNames[i])
        if( length(ind) != 0 ){
            indF[i] = ind
        }
    }

    indSig  <- which(colnames(rr)=="sig.genes")
       
    indRank <- which(colnames(rr)=="Rank")

    indBCF1 <- which(colnames(rr)=="Bonferroni for classicFisher")
    indBCF2 <- which(colnames(rr)=="Bonferroni for elimFisher")

    indBHF1 <- which(colnames(rr)=="Benjamini & Yekutieli (BY) for classicFisher")
    indBHF2 <- which(colnames(rr)=="Benjamini & Yekutieli (BY) for elimFisher")
       
    #copy allRes data 
    for( i in 1:length(names(allRes)) ){
        rr[,i] <- c(allRes[i][[1]])
    }

    #---Check        
    rr[rr[,indF[1]] =="< 1e-30",indF[1]] ="1e-30"
    rr[rr[,indF[2]] =="< 1e-30",indF[2]] ="1e-30" 
       
    if( length(indSig) != 0 ){
        rr[,indSig] <- as.character(allRes[indSig][[1]])
    }
       
    #Rank
    rank = rep(0,length(rr[,1]));
    tally=0;
       if( length(indRank) != 0 ){
           for( r in 1:length(methodNames) ){
               
               if( indF[r] != -1 ){              
                   rank = rank + log10(as.numeric(rr[,indF[r]]))
                   tally=tally+1;
               }
           }
       }
    
       rr[,indRank] = rank;
       if( tally > 0 ){
           rr[,indRank] = exp( 1/tally * rank)
       }     
   
    #---Bonferroni correction
    Nlevels = length(nodes(GOdata@graph));
      
           
       alpha <- vector(length=3)
       alpha[1] <- 0.05
       alpha[2] <- 0.01
       alpha[3] <- 0.001

       stars    <- vector(length=3)
       stars[1] <- "*"
       stars[2] <- "**"
       stars[3] <- "***"
       
       #Classic Fisher
       if( length(indBCF1) != 0 ){
           rr[,indBCF1]                            <- ""
           for( x in 1:length(alpha) ){
               rr[as.numeric(rr[,indF[1]]) < as.numeric(alpha[x]/Nlevels),indBCF1]  <- stars[x]
           }          
       }
   
       #elim Fisher
       if( length(indBCF2) != 0 ){
           rr[,indBCF2]                            <- ""
           for( x in 1:length(alpha) ){
               rr[as.numeric(rr[,indF[2]]) < as.numeric(alpha[x]/Nlevels),indBCF2]  <- stars[x]
           }           
       }

    #for Benjamini & Hochberg corrected p.values, use p.adjust package
    #Because we cannot gaurantee independency use the Benjamini and Yekutieli test instead
    #Benjamini, Y. & Yekutieli, D. The control of the false discovery rate in multiple testing
    #                              under dependency, The Annals of Statistics, 29, 4, 1165-1188 (2001).
    rr[,indBHF1] <- p.adjust(as.numeric(rr[,indF[1]]),method="BY",n=Nlevels)#length(allRes$Annotated))   
    rr[,indBHF2] <- p.adjust(as.numeric(rr[,indF[2]]),method="BY",n=Nlevels)#length(allRes$Annotated))   
       
   
    outfile <- file(sprintf("%s/%s/%s/output.csv",mainDir,ff[f],ss[s]), "w");
    write.table( rr, file=outfile, append=T, row.names=F, col.names=T, sep="\t", quote=F);
    close(outfile);

    #---format file
    hh <- matrix(NA,  sum(allRes[[4]]), 3 );
    t=1;

    hits=lapply(allRes$TERM.ID,function(x) intersect(sigGenes(GOdata),genesInTerm(GOdata, c(x))[[1]])) 
       
    for(i in 1:length(hits) ){
      
      for(j in 1:(length(as.numeric(hits[[i]]))+1) ){

        if( !is.na(hits[[i]][j]) ){    
          hh[t,1] <- as.character(allRes[[1]][i]);
          hh[t,2] <- as.character(gsub(" ","_",allRes[[2]][i]));
          hh[t,3] <- as.character(hits[[i]][j]);
          #dd[t,3] <- as.numeric(hits[[i]][j]);
          t=t+1;
        }
      }     
    }

    outfile <- file(sprintf("%s/%s/%s/flatfile.csv",mainDir,ff[f],ss[s]), "w");
    write.table( hh, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
    close(outfile);

       #cat(" ****** UNIVERSE GENE ENTREZ ID SIZE   ", length(names(geneID2GO)),  " UNIVERSE HDO ID SIZE ", length(names(TERM2geneID)),       " ****** \n " )
       #cat(" ****** MY ENTREZ GENE LIST SIZE       ", length(myInterestingGenes)," MY GENE LIST IN UNIVERSE SIZE ", sum(geneList==1),        " ****** \n " )
       #cat(" ****** NO: OF HDO IDS IN MY GENE LIST ", length(unique(dd[,3])),                                                                " ****** \n " )
       
       outfile <- file(sprintf("%s/%s/%s/info.csv",mainDir,ff[f],ss[s]),"w")

       sink(append=T,file=outfile)
       print(GOdata)
       sink();
       
       cat(" ****** UNIVERSE GENE ENTREZ ID SIZE   ", length(names(geneID2GO)),  " UNIVERSE HDO ID SIZE ", length(names(TERM2geneID)),       " ****** \n", file=outfile, append=T, sep="" )
       cat(" ****** MY ENTREZ GENE LIST SIZE       ", length(myInterestingGenes)," MY GENE LIST IN UNIVERSE SIZE ", sum(geneList==1),        " ****** \n", file=outfile, append=T, sep="")
       cat(" ****** NO: OF HDO IDS IN MY GENE LIST ", length(unique(rr[,1])),                                                                " ****** \n", file=outfile, append=T, sep="" )
       cat(" ****** NO: LEVELS FOR BONFERRONI CORRECTION ", Nlevels, file=outfile, append=T, sep="" );
       
       close(outfile);

   }

}


