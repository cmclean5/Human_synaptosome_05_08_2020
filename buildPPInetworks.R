source('setUp.R')

#---Find Largest CC
findLCC <- function(GG){

    dec <- igraph::decompose.graph(GG)
    d=1
    CC=length(V(dec[[1]]))
    for( i in 1:length(dec) ){
        if(length(V(dec[[i]])) > CC){
            d=i
            CC=length(V(dec[[i]]))
        }
    }   

    GG  <- igraph::decompose.graph(GG)[[d]]
    return(GG)

}

findTERM <- function(eatt, TERMS){

    eatt  <- as.vector(eatt)    
    found <- rep(FALSE, length(eatt))
    
    T = length(TERMS)
    for( t in 1:T ){
        temp  <- grepl(eatt[t], eatt)
        found <- as.logical(found) | as.logical(temp)        
    }

    return(found)
    
}


filterPMIDs <- function(TERMS=NULL){

    filterIDs <- c()
    
    if( !is.null(TERMS ) && length(TERMS) != 0 ){
    
        pmids <- read.delim("/afs/inf.ed.ac.uk/user/c/cmclean5/ownCloud/Synaptic_proteome/anaysis_17_05_2019/mined_PPIs/pmid_keywords.csv",sep="\t",header=T)

        indX <- list()

        for( i in 1:length(TERMS) ){
            N = length(names(indX))
            indX[[N+1]] <- grepl(TERMS[i], pmids[,4])
            names(indX)[N+1] <- sprintf("%s_title",TERMS[i])
            N = length(names(indX))
            indX[[N+1]] <- grepl(TERMS[i], pmids[,5])
            names(indX)[N+1] <- sprintf("%s_keywords",TERMS[i])
        }

        exc <- indX[[1]]
        for( i in 2:length(names(indX)) ){
            exc <- exc | indX[[i]]
        }

        filterIDs <- as.vector(pmids[exc,1])
        
    }

    return(filterIDs)
    
}

#--- add attributes to igraph edges from it raw file
addEdgeAtts <- function(GG, gg){

    ATTS = names(edge.attributes(GG))
    
    if( !is.null(ATTS) ){

        ed = get.edgelist(gg)
        M  = length(E(gg))
        ED = get.edgelist(GG)

        VALUES = list()
        
        for( a in 1:length(ATTS) ){
            VALUES[[a]] = get.edge.attribute(GG,ATTS[a],E(GG))
            names(VALUES)[a] = ATTS[a]
        }
            
        cat("\n")
        cat("scanning edges...")
        RES    = matrix("",nrow=M, ncol=length(ATTS))
        
        for( e in 1:M ){
    
            indx = (ed[e,1] == ED[,1] & ed[e,2] == ED[,2]) | (ed[e,1] == ED[,2] & ed[e,2] == ED[,1])

            for( a in 1:length(ATTS) ){
                
                res = VALUES[[a]][indx]

                if( res != "" ){
                    res <- unique(res)
                    if( length(res) == 1 ){
                        RES[e,a] <- res
                    } else {
                        RES[e,a] <- paste(as.character(res),collapse=';') 
                    }
                }
            }
        }

        cat("done.\n")
        
        for( a in 1:length(ATTS) ){
            gg <- igraph::set.edge.attribute(gg,ATTS[a],E(gg),as.character(RES[,a]))
        }

    }
           
    return(gg)
    
}

#---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("datasets",DIRS)]
OUT[2] <- DIRS[grepl("GeneSets",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
#eldir <- sprintf("%s/%s",dataDIR,subDIR[S])
#if( !file_test("-d",eldir) ){
#    dir.create(eldir)
#}

gsdir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",gsdir) ){
    dir.create(gsdir)
}

grdir <- sprintf("%s/%s",OUT[3],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---

#---RAW DATASET TO READ IN
ff = read.table(sprintf("/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/Human_synaptosome_05_08_2020/datasets/_safe/%s.csv",subDIR[S]),sep="\t",header=T)

#--- build raw graph
GG <- graph.data.frame(ff[,1:2],directed=F)

run=0
if( run ){

#--- set and fill edge attributes
kw <- read.delim("/afs/inf.ed.ac.uk/user/c/cmclean5/ownCloud/Synaptic_proteome/analysis_17_05_2019/mined_PPIs/MINED_PUBMED_PPIs/pmid_keywords.csv",sep="\t",header=T)

GG = set.edge.attribute(GG,"METHOD",E(GG), as.character(ff[,3]))
GG = set.edge.attribute(GG,"TYPE",E(GG), as.character(ff[,7]))

PMIDS = ifelse(!grepl("unassigned",ff[,4]), sprintf("PMID:%s",ff[,4]), ff[,4])
GG = set.edge.attribute(GG,"PUBMED",E(GG), PMIDS)

YEARS = kw[match(gsub("PMID:","",E(GG)$PUBMED),kw[,1]),3]
YEARS = ifelse(is.na(YEARS),"na",YEARS)
GG = set.edge.attribute(GG,"YEAR",E(GG), YEARS)
#---

}
    
#--- build igraph, removing multiple edges and loops
gg <- igraph::simplify(GG,remove.multiple=T,remove.loops=T)

#---Print Graph of interest
cat("\n")
cat("build graph... ")

#---Find Largest CC
gg  <- findLCC(gg)

cat("done.")
cat("\n")

#---Write out the gene-set (Human Entrez IDs)
ids <- V(gg)$name;
write.table( ids, file=sprintf("%s/%s.csv",gsdir,subDIR[S]),append=T,row.names=F,col.names=F,sep="\t",quote=F)

#---Write out the gene-list (Human Entrez IDs)
ed  <- get.edgelist(gg)
ed  <- cbind(ed[,1],ed[,2], rep(1,length(ids)))
colnames(ed) <- c("entA","entB","We")
write.table( ed, file=sprintf("%s/%s.csv",OUT[1],subDIR[S]),append=T,row.names=F,col.names=T,sep="\t",quote=F)


#--- Add attributes to each igraph edge, from raw file
cat("\n")
cat("add edge attributes to graph... ")

gg <- addEdgeAtts(GG,gg)

cat("Finished.\n")


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")
