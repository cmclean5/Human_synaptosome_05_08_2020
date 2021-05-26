
source('../setUp.R')

removeVertexTerm <- function(GG,NAME){

    if( !is.null(igraph::get.vertex.attribute(GG,NAME)) ){
        GG <- igraph::remove.vertex.attribute(GG,name=NAME)
    }

    if( !is.null(igraph::get.vertex.attribute(GG,gsub("_","",NAME))) ){    
        GG <- igraph::remove.vertex.attribute(GG,name=gsub("_","",NAME))
    }

    return(GG)
    
}

addVertexData <- function(GG,TITLE="",OO, ATTRS=c("K","LABEL","Pmax")){

    for( i in 1:length(ATTRS) ){
    
        indx  <- which(ATTRS[i]==colnames(OO))

        val = as.vector(OO[,indx])
        
        if( is.numeric(val) ){        
            val   <- as.numeric(val)
        }

        if( is.character(val) ){
            val   <- as.character(val)
        }
        
        Vlab  <- sprintf("%s_%s",ATTRS[i],TITLE)
        GG    <- removeVertexTerm(GG,Vlab)
        GG    <- igraph::set.vertex.attribute(GG,Vlab,V(GG),val)

    }

    return(GG)
             
}

saveStudyData <- function(GG,STUDY,TITLE){

    N    = ncol(STUDY$NP)
    oo   = cbind(STUDY$DF,STUDY$VANNO,do.call(pmax,STUDY$NP[,2:N]))
    colnames(oo) = c(colnames(STUDY$DF),"LABEL","Pmax")
    oo$K = oo$K+1

    indx = match(V(GG)$name,oo$name)
    oo   = oo[indx,]

    GG = addVertexData(GG=GG,TITLE=TITLE,OO=oo)
    
    return(GG)
    
}

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
ids <- V(gg)$name;

run1=0
run2=0
run3=0
if( run1 ){

#---Load SBM study results
STUDY <- readRDS(sprintf("%s/%s/SBMstudies.RDS","RESULTS",subDIR[S]))

#--- add study to graph
for( i in 1:length(STUDY) ){
    gg = saveStudyData(GG=gg,STUDY=STUDY[[i]],TITLE=names(STUDY)[i])
}

}

if( run2 ){
#---Load SBM study results
STUDY <- readRDS(sprintf("%s/%s/SBMstudies_C46.RDS","RESULTS",subDIR[S]))

#--- add study to graph
gg = saveStudyData(GG=gg,STUDY=STUDY[[5]],TITLE=names(STUDY)[5])
}

if( run3 ){
#---Load SBM study results
STUDY <- readRDS(sprintf("%s/%s/SBMstudies_None_46.RDS","RESULTS",subDIR[S]))

#--- add study to graph
gg = saveStudyData(GG=gg,STUDY=STUDY[[6]],TITLE=names(STUDY)[6])
}

colorNodes <- 0
#colorNodes <- 1

if( colorNodes ){

    Ks <- c("31","43","44")
    
    cols <- matrix(NA,ncol=2,nrow=4)
    cols[1,1] = "#CCCCCC"; cols[1,2] = "#E5E5E5"; #grey
    cols[2,1] = "#FFA500"; cols[2,2] = "#FFD280"; #orange
    cols[3,1] = "#FF0000"; cols[3,2] = "#FF8080"; #red
    cols[4,1] = "#00BFFF"; cols[4,2] = "#80DFFF"; #blue

    oo    <- cbind(V(gg)$name,V(gg)$LABELADHTN46,V(gg)$KADHTN46,rep("",length(V(gg))))

    Kindx = match(oo[,3],Ks)
    Kindx = ifelse(is.na(Kindx),FALSE,TRUE)

    oo[grepl("other",oo[,2]) & Kindx == TRUE,4]  = cols[1,1]
    oo[grepl("other",oo[,2]) & Kindx == FALSE,4] = cols[1,2]

    oo[grepl("AD",oo[,2]) & !grepl("AD&HTN",oo[,2]) & Kindx == TRUE,4]  = cols[2,1]
    oo[grepl("AD",oo[,2]) & !grepl("AD&HTN",oo[,2]) & Kindx == FALSE,4] = cols[2,2]

    oo[grepl("HTN",oo[,2]) & !grepl("AD&HTN",oo[,2])& Kindx == TRUE,4]  = cols[3,1]
    oo[grepl("HTN",oo[,2]) & !grepl("AD&HTN",oo[,2])& Kindx == FALSE,4] = cols[3,2]

    oo[grepl("AD&HTN",oo[,2]) & Kindx == TRUE,4]  = cols[4,1]
    oo[grepl("AD&HTN",oo[,2]) & Kindx == FALSE,4] = cols[4,2]    
    
    gg    <- removeVertexTerm(gg,"color")
    gg    <- igraph::set.vertex.attribute(gg,"color",V(gg),oo[,4])

}


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")

