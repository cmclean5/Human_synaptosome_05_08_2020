##reset
rm(list=ls())

library(rEnrich)
library(igraph)

## load graph
gg = read.graph("PPI_Presynaptic.gml",format="gml")

## build cluster membership data.frame
alg  = c("louvain2","Spectral")
ids  = igraph::get.vertex.attribute(gg,"name",V(gg))
coms = igraph::get.vertex.attribute(gg,alg[2],V(gg))

membership = as.data.frame(cbind(ids,coms))

## load annotation flat file
anno = read.delim("flatfile.go.BP.csv",skip=1,sep="\t",header=F)
anno = as.data.frame(anno)

## load clustering and annotation data
rEnrich::load(x=membership,anno=anno)

## run enrichment analysis on loaded data
rEnrich::run()

## get enrichment values 
res = rEnrich::getResults(printTwoSided=1,
                          usePrintAlt=1,
                          usePrintID=1)

## save results
write.table(res,"permute_p_values_GOBP.csv", sep="\t", row.names=F, col.names=T, quote=F)
