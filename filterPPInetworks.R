source('setUp.R')
require(igraph)
require(stringr)

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
        temp  <- grepl(TERMS[t], eatt)
        found <- as.logical(found) | as.logical(temp)        
    }

    return(found)
    
}

countTERM <- function(eatt, TERM){

    eatt  <- as.vector(eatt)
    tally <- rep("",length(eatt))
    
    for( p in 1:length(eatt) ){

        count = 0
        count = as.numeric(str_count(eatt[p],TERM))

        tally[p] = count
        
    }

    return(tally)
    
}

matching_edge <- function(g1,e,g2) {
  # given an edge e in g1, return the corresponding edge in g2
  name1 <- V(g1)[get.edges(g1,e)[1,1]]$name
  name2 <- V(g1)[get.edges(g1,e)[1,2]]$name
  len1 <- which(V(g2)$name == name1)
  len2 <- which(V(g2)$name == name2)
    if( length(len1) == 0 || length(len2) == 0 ){
        return(0)
    } else {
        return(get.edge.ids(g2,c(name1,name2),directed=FALSE,error=FALSE))
    }
}

maxLSi <- function( XX, BASE=0 ){

    XX  <- as.vector(as.numeric(XX))
    
    XXo <- XX[XX != 0]
    if( BASE == 2 ){
        return ( -sum( XXo * log2(XXo)) )
    }

    if( BASE == 10 ){
        return ( -sum( XXo * log10(XXo)) )
    }

    if( BASE == 0 ){
        return ( -sum( XXo * log(XXo)) )
    }
    
}


graphEntropy <- function(GG){

#--- initial entropy rate
V    <- length(V(GG))
E    <- length(E(GG))
ki   <- as.vector(igraph::degree(GG))
Kbar <- mean(ki)

#--- get adjacency matrix for graph
A    <- get.adjacency(GG)

#--- get leading eigenvalue and vector
R     <- eigen(A)
Rindx <- which.max(R$values)
gamma <- R$values[Rindx]
nu    <- R$vectors[,Rindx]

#--- calculate max entropy rate, maxSr
Pij   <- (A * nu) / (gamma * nu)
Pi    <- ki/(2*E)
maxSr <- sum(as.vector(Pi * apply(Pij,1,maxLSi,BASE=0)))


#--- calculate initial configuration
Norm <- as.numeric(V*Kbar) #as.numeric(2*E)
SRo  <- as.numeric(1/Norm)*sum(ki*log(ki))

if( maxSr == 0 ){ maxSr = 1 }

return( SRo/maxSr )
    
}

#---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- DIRS[grepl("Clustering",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
cldir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",cldir) ){
    dir.create(cldir)
}

grdir <- sprintf("%s/%s",OUT[3],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---

ownCloud <- "/afs/inf.ed.ac.uk/user/c/cmclean5/ownCloud/Synaptic_proteome/anaysis_11_04_2019/Graphs/"


#---READ IN PRIOR graph
seed <- read.graph(sprintf("%s/prior_20_11_2017/%s/%s.gml",ownCloud, subDIR[S],subDIR[S]),format="gml")

#---READ IN unfiltered graph
gg <- read.graph(sprintf("%s/unfiltered_PPIs/%s/%s.gml",ownCloud, subDIR[S],subDIR[S]),format="gml")
ids <- V(gg)$name

#--- Filter out interactions reported only from BioPlex    
#gg <- subgraph.edges(gg,which("PMID:28514442"!=as.vector(E(gg)$PUBMED)))
#gg <- findLCC(gg)

##seed <- gg %s% seed
##seed = findLCC(seed)

#---filter out single studies, which have only "MI:0915" TYPE and "MI:0004" METHOD terms
xx    = cbind(countTERM(E(gg)$PUBMED,"PMID:"),E(gg)$TYPE == "MI:0915", E(gg)$METHOD == "MI:0004")
eid   = as.numeric(xx[,1]) == 1 & as.logical(xx[,2]) == TRUE & as.logical(xx[,3]) == TRUE

#---High quality edges filtered from gg
highQ = subgraph.edges(gg,seq(1,length(E(gg)),1)[!eid])
edH   = get.edgelist(highQ)

#---Low quality edges, i.e. the remaining edges from gg
lowQ  = subgraph.edges(gg,seq(1,length(E(gg)),1)[eid])
edL   = get.edgelist(lowQ)

#---Find those genes not in the seed network, but in the unfiltered gg network     
gn = match(ids,V(seed)$name)
gn = ifelse(is.na(gn),TRUE, FALSE)
gn = ids[gn]

#---Find edges between these missing genes in the good and bad networks of unused edges
inda  = match(edH[,1],gn)
indb  = match(edH[,2],gn)
indx1 = !is.na(inda) | !is.na(indb)

extraH = subgraph.edges(highQ,seq(1,length(E(highQ)),1)[indx1])

indc  = match(edL[,1],gn)
indd  = match(edL[,2],gn)
indx2 = !is.na(indc) | !is.na(indd)

extraL = subgraph.edges(lowQ,seq(1,length(E(lowQ)),1)[indx2])

#--- add extra edges, and missing genes, to our seed network
new = seed %u% extraH %u% extraL

#--- Create the final graph
final = graph_from_edgelist(el=get.edgelist(new),directed=F)

#--- double-check
final <- igraph::simplify(final,remove.multiple=T,remove.loops=T)
final <- findLCC(final)
#---


#---add vertex attributes to final
if( is.null(igraph::get.vertex.attribute(final,"GeneName")) ){
    igraph::set.vertex.attribute(final,"GeneName",V(final),"")
    V(final)$GeneName = V(gg)$GeneName[match(V(final)$name,V(gg)$name)]
} else {
    final <- igraph::remove.vertex.attribute(final,"GeneName")
    igraph::set.vertex.attribute(final,"GeneName",V(final),"")
    V(final)$GeneName = V(gg)$GeneName[match(V(final)$name,V(gg)$name)]
}

if( is.null(igraph::get.vertex.attribute(final,"UniProt")) ){
    igraph::set.vertex.attribute(final,"UniProt",V(final),"")
    V(final)$UniProt = V(gg)$UniProt[match(V(final)$name,V(gg)$name)]
} else {
    final <- igraph::remove.vertex.attribute(final,"UniProt")
    igraph::set.vertex.attribute(final,"UniProt",V(final),"")
    V(final)$UniProt = V(gg)$UniProt[match(V(final)$name,V(gg)$name)]
}
#---

#---add edge attributes to final
atts = matrix("", ncol=4,nrow=length(E(final)))
for (e in E(final)) {
    eg <- matching_edge(final,e,gg)
    if( eg != 0 ){
        atts[e,1] = E(gg)[eg]$METHOD
        atts[e,2] = E(gg)[eg]$PUBMED
        atts[e,3] = E(gg)[eg]$TYPE
        atts[e,4] = E(gg)[eg]$YEAR
    }
}

if( is.null(igraph::get.edge.attribute(final,"METHOD")) ){
    final <- igraph::set.edge.attribute(final,"METHOD",E(final),as.character(atts[,1]))
} else {
    final <- igraph:: remove.edge.attribute(final,"METHOD")
    final <- igraph::set.edge.attribute(final,"METHOD",E(final),as.character(atts[,1]))
}

if( is.null(igraph::get.edge.attribute(final,"PUBMED")) ){
    final <- igraph::set.edge.attribute(final,"PUBMED",E(final),as.character(atts[,2]))
} else {
    final <- igraph:: remove.edge.attribute(final,"PUBMED")
    final <- igraph::set.edge.attribute(final,"PUBMED",E(final),as.character(atts[,2]))
}

if( is.null(igraph::get.edge.attribute(final,"TYPE")) ){
    final <- igraph::set.edge.attribute(final,"TYPE",E(final),as.character(atts[,3]))
} else {
    final <- igraph:: remove.edge.attribute(final,"TYPE")
    final <- igraph::set.edge.attribute(final,"TYPE",E(final),as.character(atts[,3]))
}

if( is.null(igraph::get.edge.attribute(final,"YEAR")) ){
    final <- igraph::set.edge.attribute(final,"YEAR",E(final),as.character(atts[,4]))
} else {
    final <- igraph:: remove.edge.attribute(final,"YEAR")
    final <- igraph::set.edge.attribute(final,"YEAR",E(final),as.character(atts[,4]))
}

#---


##---Write .gml graph to file
igraph::write.graph(final, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(final, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")



run=FALSE


### Potential problematic studies  ###
## "PMID:28611215", "PMID:29128334", "PMID:30021884", "PMID:30442766"
###                                ###

#---TEST
if (run ){

ownCloud <- "/afs/inf.ed.ac.uk/user/c/cmclean5/ownCloud/Synaptic_proteome/anaysis_11_04_2019/Graphs/"
     
#---READ IN unfiltered graph
gg  <- igraph::read.graph(sprintf("%s/unfiltered_PPIs/%s/%s.gml",ownCloud, subDIR[S],subDIR[S]),format="gml")
ids <- V(gg)$name

#--- Filter out interactions reported only from BioPlex    
#gg <- subgraph.edges(gg,which("PMID:28514442"!=as.vector(E(gg)$PUBMED)))
#gg <- findLCC(gg)
    
#---filter out single studies, which have only "MI:0915" TYPE and "MI:0004" METHOD terms
xx    = cbind(countTERM(E(gg)$PUBMED,"PMID:"),E(gg)$TYPE == "MI:0915", E(gg)$METHOD == "MI:0004")
eid   = as.numeric(xx[,1]) == 1 & as.logical(xx[,2]) == TRUE & as.logical(xx[,3]) == TRUE

#---Build our seed network from the filtered gg
seed  = subgraph.edges(gg,seq(1,length(E(gg)),1)[!eid])

#---Keep those remaining, unused, interactions from the filtered gg not used in seed network
lowQ = subgraph.edges(gg,seq(1,length(E(gg)),1)[eid])
edL  = get.edgelist(lowQ)

#---Find those genes not in the seed network, but in the unfiltered gg network     
gn = match(ids,V(seed)$name)
gn = ifelse(is.na(gn),TRUE, FALSE)
gn = ids[gn]

#---Find edges between these missing genes in the network of unused edges  
inda  = match(ed[,1],gn) & match(ed[,2],V(seed)$name)
indb  = match(ed[,2],gn) & match(ed[,1],V(seed)$name)
indx1 = !is.na(inda) | !is.na(indb)

indc  = match(ed[,1],gn)
indd  = match(ed[,2],gn)
indx2 = !is.na(indc) | !is.na(indd)

indx  = indx1 | indx2
    
extra = subgraph.edges(lowQ,seq(1,length(E(lowQ)),1)[indx])


#--- add extra edges, and missing gene, to our seed network
new = seed %u% extra

#--- Create the final graph
final = graph_from_edgelist(el=get.edgelist(new),directed=F)

#--- double-check
final <- findLCC(final)
final <- igraph::simplify(final,remove.multiple=T,remove.loops=T)
#---

#---add vertex attributes to final
    if( is.null(igraph::get.vertex.attribute(final,"GeneName")) ){
        igraph::set.vertex.attribute(final,"GeneName",V(final),"")
        V(final)$GeneName = V(gg)$GeneName[match(V(final)$name,V(gg)$name)]
    } else {
        final <- igraph::remove.vertex.attribute(final,"GeneName")
        igraph::set.vertex.attribute(final,"GeneName",V(final),"")
        V(final)$GeneName = V(gg)$GeneName[match(V(final)$name,V(gg)$name)]
    }

    if( is.null(igraph::get.vertex.attribute(final,"UniProt")) ){
        igraph::set.vertex.attribute(final,"UniProt",V(final),"")
        V(final)$UniProt = V(gg)$UniProt[match(V(final)$name,V(gg)$name)]
    } else {
        final <- igraph::remove.vertex.attribute(final,"UniProt")
        igraph::set.vertex.attribute(final,"UniProt",V(final),"")
        V(final)$UniProt = V(gg)$UniProt[match(V(final)$name,V(gg)$name)]
    }
#---


    if( is.null(igraph::get.edge.attribute(final,"METHOD")) ){
        final <- igraph::set.edge.attribute(final,"METHOD",E(final),as.character(E(new)$METHOD_1))
    } else {
        final <- igraph:: remove.edge.attribute(final,"METHOD")
        final <- igraph::set.edge.attribute(final,"METHOD",E(final),as.character(E(new)$METHOD_1))
    }

    if( is.null(igraph::get.edge.attribute(final,"PUBMED")) ){
        final <- igraph::set.edge.attribute(final,"PUBMED",E(final),as.character(E(new)$PUBMED_1))
    } else {
        final <- igraph:: remove.edge.attribute(final,"PUBMED")
        final <- igraph::set.edge.attribute(final,"PUBMED",E(final),as.character(E(new)$PUBMED_1))
    }

    if( is.null(igraph::get.edge.attribute(final,"TYPE")) ){
        final <- igraph::set.edge.attribute(final,"TYPE",E(final),as.character(E(new)$TYPE_1))
    } else {
        final <- igraph:: remove.edge.attribute(final,"TYPE")
        final <- igraph::set.edge.attribute(final,"TYPE",E(final),as.character(E(new)$TYPE_1))
    }

    if( is.null(igraph::get.edge.attribute(final,"YEAR")) ){
        final <- igraph::set.edge.attribute(final,"YEAR",E(final),as.character(E(new)$YEAR_1))
    } else {
        final <- igraph:: remove.edge.attribute(final,"YEAR")
        final <- igraph::set.edge.attribute(final,"YEAR",E(final),as.character(E(new)$YEAR_1))
    }
    
#---
 
##---Write .gml graph to file
#igraph::write.graph(final, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
#igraph::write.graph(final, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")
  
    
}



