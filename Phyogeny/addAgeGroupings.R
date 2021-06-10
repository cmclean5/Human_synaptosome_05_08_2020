source('../setUp.R')

#---OUT Dir
OUT    <- vector(length=2)
OUT[1] <- DIRS[grepl("Phyogeny",DIRS)]
OUT[2] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
grdir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---


#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
ids <- V(gg)$name;

#---READ IN EVO AGE GROUPS for proteins
ag <- read.delim(sprintf("%s/synaptic_proteins_age_group.csv",OUT[1]),sep="\t",header=T)

#---Add Gene Names
#--- mapping lists from Entrez IDs to gene names
if( is.null(igraph::get.vertex.attribute(gg,"AgeGroup")) ){
        
    igraph::set.vertex.attribute(gg,"AgeGroup",V(gg),"")
    V(gg)$AgeGroup = ag[match(ids,ag[,1]),2]

}


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")

