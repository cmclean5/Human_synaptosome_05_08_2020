
source('../setUp.R')


#---OUT Dir
OUT    <- vector(length=1)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]


grdir <- sprintf("%s/%s",OUT[1],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---

#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")

oo = cbind(V(gg)$name,V(gg)$LABELADHTN46,V(gg)$KADHTN46)

ad       = grepl("AD",oo[,2])
htn      = grepl("HTN",oo[,2])
adANDhtn = grepl("AD&HTN",oo[,2])
other    = grepl("other",oo[,2])

adONLY   = ad & !adANDhtn
htnONLY  = htn & !adANDhtn
adORhtn  = ad | adANDhtn | htn


oo1 = cbind(rep(1,table(adONLY)[2]),  rep("ADonly",  table(adONLY)[2]),  oo[adONLY,1])
oo2 = cbind(rep(2,table(htnONLY)[2]), rep("HTNonly", table(htnONLY)[2]), oo[htnONLY,1])
oo3 = cbind(rep(3,table(adANDhtn)[2]),rep("ADandHTN",table(adANDhtn)[2]),oo[adANDhtn,1])
oo4 = cbind(rep(4,table(adORhtn)[2]), rep("ADorHTN", table(adORhtn)[2]), oo[adORhtn,1])
oo5 = cbind(rep(5,table(other)[2]),   rep("other",   table(other)[2]),   oo[other,1])

annos = rbind(oo1,oo2,oo3,oo4,oo5)

write.table(annos,"SBM_anno.csv",sep="\t",row.names=F,col.names=F,quote=F)

cc = cbind(oo[,1],oo[,3])

write.table(cc,"SBM_clustering.csv",sep="\t",row.names=F,col.names=F,quote=F)
