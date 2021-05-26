#add clustering results for C++ code
# Spectral
# Geodesic
# SVI
# SPICi

source('setUp.R')

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

#--Read PPI network graph
files <- list.files(sprintf("%s/%s",OUT[3],subDIR[S]))
files <- files[grepl(".gml" ,files,fixed=T)]

#---load graph
gg  <- igraph::read.graph(sprintf("%s/%s/%s",OUT[3],subDIR[S],files),format="gml")


#---Add Gene Names
#--- mapping lists from gene names to ids
ids = V(gg)$name

if( is.null(igraph::get.vertex.attribute(gg,"GeneName")) ){

    gn <- mapIds(org.Hs.eg.db,ids,column="SYMBOL",keytype="ENTREZID")

    igraph::set.vertex.attribute(gg,"GeneName",V(gg),"")
    V(gg)$GeneName = gn

}

addCnmin1  = 0
add0.1per  = 0
add0.25per = 0
add0.5per  = 0
add1per    = 0
add2.5per  = 0
add5per    = 0

CPP = c("Spectral")

#---loop through C++ algorithms
for( a in 1:length(CPP) ){    

    if( CPP[a] == "Geodesic" ){

        str <- sprintf("%s/%s/communities_cytoscape.txt",CPP[a],subDIR[S])

        if(file.exists(str)){
            cc <- read.table(str,sep="\t",header=F,skip=1)
            gg <- igraph::set.vertex.attribute(gg,CPP[a],V(gg),cc[match(V(gg)$name,cc[,1]),2])

        }
    }

    if( CPP[a] == "Spectral" ){

        str <- sprintf("%s/%s/communities_cytoscape.txt",CPP[a],subDIR[S])

        if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,CPP[a],V(gg),cc[match(V(gg)$name,cc[,1]),2])
            }
                
        
        if( addCnmin1 ){
        
            str <- sprintf("%s/%s_Cnmin1/communities_cytoscape.txt",CPP[a],subDIR[S])

            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,CPP[a],V(gg),cc[match(V(gg)$name,cc[,1]),2])

            }
        }
        
        if( add0.1per ) {
            
            str <- sprintf("%s/%s_0.1per/communities_cytoscape.txt",CPP[a],subDIR[S])        
        
            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,"Spectral_0.1per",V(gg),cc[match(V(gg)$name,cc[,1]),2])
            }
        }
        
        
        if( add0.25per ) {
            
            str <- sprintf("%s/%s_0.25per/communities_cytoscape.txt",CPP[a],subDIR[S])        
        
            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,"Spectral_0.25per",V(gg),cc[match(V(gg)$name,cc[,1]),2])
            }
        }
        
        
        if( add0.5per ) {
            
            str <- sprintf("%s/%s_0.5per/communities_cytoscape.txt",CPP[a],subDIR[S])        
        
            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,"Spectral_0.5per",V(gg),cc[match(V(gg)$name,cc[,1]),2])
            }
        }
        
        
        if( add1per ){
        
            str <- sprintf("%s/%s_1per/communities_cytoscape.txt",CPP[a],subDIR[S])        
        
            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,"Spectral_1per",V(gg),cc[match(V(gg)$name,cc[,1]),2])
            }
        }

        
        if( add2.5per ){
        
            str <- sprintf("%s/%s_2.5per/communities_cytoscape.txt",CPP[a],subDIR[S])        
        
            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,"Spectral_2.5per",V(gg),cc[match(V(gg)$name,cc[,1]),2])
            }
        }


        if( add5per ){
        
            str <- sprintf("%s/%s_5per/communities_cytoscape.txt",CPP[a],subDIR[S])        
        
            if(file.exists(str)){
                cc <- read.table(str,sep="\t",header=F, skip=1)
                gg <- igraph::set.vertex.attribute(gg,"Spectral_5per",V(gg),cc[match(V(gg)$name,cc[,1]),2])

            }
        }

        
    }

    

    if( CPP[a] == "SPICi" ){

        str <- sprintf("%s/%s/SPICi_communities_cytoscape.csv",CPP[a],subDIR[S])

        if(file.exists(str)){

            cc   <- read.table(str,sep="\t",header=F,fill=T)
            
            CC   <- data.frame(a=as.numeric(),b=as.numeric())
            for( i in 1:length(rownames(cc)) ){
                CC <- rbind(CC, data.frame(as.numeric(cc[i,]),as.numeric(rep(i,length(as.numeric(cc[2,]))))))
            }
            CC <- CC[!is.na(CC[,1]),]
            gg <- igraph::set.vertex.attribute(gg,CPP[a],V(gg),CC[match(V(gg)$name,CC[,1]),2])
        }        
    }


     if( CPP[a] == "SVI" ){

         str <- sprintf("%s/%s/network.gml",CPP[a],subDIR[S])

         if(file.exists(str)){
             svi <- read.graph(str,format="gml")
              gg <- igraph::set.vertex.attribute(gg,CPP[a],V(gg),(V(svi)[match(V(gg)$name,V(svi)$extid)]$group+1))
         }         
     }

    #write clustering to file
    cc <- igraph::get.vertex.attribute(gg,CPP[a],V(gg))
    if( !is.null(cc) ){
        oo <- cbind(V(gg)$name,cc)
        outfile <- file(sprintf("%s/%s_communities.csv",cldir,CPP[a]),"w")
        cat("#communities",file=outfile,"\n")
        write.table(oo, file=outfile, append=T, row.names=F, col.names=F, sep="\t", quote=F);
        close(outfile);
    }
    
}


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")

        

