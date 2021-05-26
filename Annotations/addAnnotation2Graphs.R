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


loopOverFiles <- function(GG, FILES, NAME, IDS, addIDS){

    for( f in 1:length(FILES) ){        
        
        GG <- removeVertexTerm(GG, NAME[f])

        if( addIDS ){
            GG <- removeVertexTerm(GG, sprintf("%s_ID",NAME[f]))
        }
        
        if( file.exists(sprintf("./%s",FILES[f])) ){            
            
            #--- Set Disease (geneRIF db) attributes in .gml graph
            igraph::set.vertex.attribute(GG, NAME[f],V(GG),"")

            annoF    <- read.table(FILES[f],sep="\t",skip=1,strip.white=T,quote="")
            annoFIDS <- as.character(annoF[,3])

            typeF   <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))  

            oo     <- matrix("",ncol=3,nrow=length(IDS))
            oo[,1] <- IDS
                        
            for( i in 1:length(IDS) ){

                ind1 = which(annoFIDS==IDS[i])

                Str1 <- "";
                Str2 <- "";                
                
                if( length(ind1) != 0 ){
                    
                    if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
                    else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }
                    
                    if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
                    else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }
                    
                    if( grepl(COLLAPSE[c], Str1) ){
                        Str1 <- strsplit(Str1,COLLAPSE[c])[[1]]
                        Str1 <- unique(Str1)
                        if( length(Str1) > 1 ){
                            Str1 <- paste(as.character(Str1),collapse=COLLAPSE[c])
                        }
                    }

                    if( grepl(COLLAPSE[c], Str2) ){
                        Str2 <- strsplit(Str2,COLLAPSE[c])[[1]]
                        Str2 <- unique(Str2)
                        if( length(Str2) > 1 ){
                            Str2 <- paste(as.character(Str2),collapse=COLLAPSE[c])
                        }
                    }                    
                
                    oo[i,2] <- Str1
                    oo[i,3] <- Str2
                
                }

            }
                
    GG <- igraph::set.vertex.attribute(GG,NAME[f],V(GG),as.character(oo[,2]))        

            if( addIDS ){
                GG <- igraph::set.vertex.attribute(GG,sprintf("%s_ID",NAME[f]),V(GG),as.character(oo[,3]))
                }
        
            
        }
    }


    return(GG)
    
}

#---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- DIRS[grepl("Clustering",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]
OUT[4] <- DIRS[grepl("Consensus",DIRS)]

#check Dir
grdir <- sprintf("%s/%s",OUT[3],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}


#---OPTIONS
runAnno     <- vector(length=13)
runAnno[1]  <- 0 #Add GeneNames
runAnno[2]  <- 0 #Add topOnto_ovg
runAnno[3]  <- 0 #Add topOnto_ov_P140papers
runAnno[4]  <- 0 #Add SCHanno synaptic functional groups
runAnno[5]  <- 0 #Add CHUA synaptic functional groups
runAnno[6]  <- 0 #Add InterPro Family and Domain synaptic functional groups
runAnno[7]  <- 0 #Add Core PSD and Pre-synpatic compartmental genes
runAnno[8]  <- 0 #Add Bridgeness Regions for each algorithm
runAnno[9]  <- 1 #Add GO MF
runAnno[10] <- 1 #Add GO BP
runAnno[11] <- 1 #Add GO CC
runAnno[12] <- 0 #Add celltypes
runAnno[13] <- 0 #Add pathways

    #--Read PPI network graph
    files <- list.files(sprintf("%s/%s",OUT[3],subDIR[S]))
    files <- files[grepl(".gml" ,files,fixed=T)]

    #---load graph
    gg  <- igraph::read.graph(sprintf("%s/%s/%s",OUT[3],subDIR[S],files),format="gml")

    #--- mapping lists from gene names to ids
    ids = V(gg)$name

 
if( runAnno[1] ){

    gn <- mapIds(org.Hs.eg.db,ids,column="SYMBOL",keytype="ENTREZID")

    gg <- removeVertexTerm(gg,"GeneName")
    
    igraph::set.vertex.attribute(gg,"GeneName",V(gg),"")
    V(gg)$GeneName = gn

}
    

if( runAnno[2] ){

    gg <- removeVertexTerm(gg,"TopOnto_OVG")
    gg <- removeVertexTerm(gg,"TopOnto_OVG_HDO_ID")    
     
    #--- Set Disease (geneRIF db) attributes in .gml graph
    igraph::set.vertex.attribute(gg,"TopOnto_OVG",V(gg),"")
    igraph::set.vertex.attribute(gg,"TopOnto_OVG_HDO_ID",V(gg),"")

    #par    <- read.table("flatfile_human_gene2HDO.parentTerm.csv",sep="\t",skip=1,strip.white=T,quote="")
    dis    <- read.table("flatfile_human_gene2HDO.csv",sep="\t",skip=1,header=F,strip.white=T,quote="")

    #dis    <- rbind(dis,par)
    
    disIDS <- dis[,3]

for( i in 1:length(ids) ){

  ind1 = which(disIDS==ids[i])
  
    Str1 <- "";
    Str2 <- "";
  
  if( length(ind1) != 0 ){

    disv <- as.vector(dis[ind1,1]); 

    indx <- match(disv,disn)

    for( j in 1:length(disv) ){

        if( !is.na(indx[j]) ){
        
            if( Str1 == "" ) { Str1 <- as.character(dtype[indx[j]]) }
            else {
                Str1 <- paste(c(Str1,as.character(dtype[indx[j]])),collapse=COLLAPSE[c]) }
        
            if( Str2 == "" ) { Str2 <- as.character(disn[indx[j]]) }
            else {
                Str2 <- paste(c(Str2,as.character(disn[indx[j]])),collapse=COLLAPSE[c]) }
        }

    }
  }

    Str1 = paste(unique(strsplit(Str1,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])
    Str2 = paste(unique(strsplit(Str2,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])
    
    V(gg)[i]$TopOnto_OVG = as.character(Str1);
    V(gg)[i]$TopOnto_OVG_HDO_ID = as.character(Str2);

  
}
}
#---

if( runAnno[3] ){

    gg <- removeVertexTerm(gg,"TopOnto_OV_PAPERS")
    gg <- removeVertexTerm(gg,"TopOnto_OV_PAPERS_HDO_ID")
    
    #--- Set Disease (geneRIF db) attributes in .gml graph
    igraph::set.vertex.attribute(gg,"TopOnto_OV_PAPERS",V(gg),"")
    igraph::set.vertex.attribute(gg,"TopOnto_OVG_PAPERS_HDO_ID",V(gg),"")

    par    <- read.table("flatfile_human_gene2HDO.parentTerm.csv",sep="\t",skip=1,strip.white=T,quote="")
    dis    <- read.table("flatfile_human_ov_PAPERS.csv",sep="\t",skip=1,strip.white=T,quote="")

    dis    <- rbind(dis,par)
    
    disIDS <- dis[,3]

for( i in 1:length(ids) ){

  ind1 = which(disIDS==ids[i])
  
    Str1 <- "";
    Str2 <- "";
  
    if( length(ind1) != 0 ){

    disv <- as.vector(dis[ind1,1]); 

    indx <- match(disv,disn)

    for( j in 1:length(disv) ){

        if( !is.na(indx[j]) ){
        
            if( Str1 == "" ) { Str1 <- as.character(dtype[indx[j]]) }
            else {
                Str1 <- paste(c(Str1,as.character(dtype[indx[j]])),collapse=COLLAPSE[c]) }
        
            if( Str2 == "" ) { Str2 <- as.character(disn[indx[j]]) }
            else {
                Str2 <- paste(c(Str2,as.character(disn[indx[j]])),collapse=COLLAPSE[c]) }
        }

    }
  }

    Str1 = paste(unique(strsplit(Str1,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])
    Str2 = paste(unique(strsplit(Str2,COLLAPSE[c])[[1]]),collapse=COLLAPSE[c])
    
    V(gg)[i]$TopOnto_OV_PAPERS= as.character(Str1);
    V(gg)[i]$TopOnto_OV_PAPERS_HDO_ID = as.character(Str2);

  
}   
}

if( runAnno[4] ){

    gg <- removeVertexTerm(gg,"SCHanno")
    
    #--- Set Family attributes in .gml graph
    igraph::set.vertex.attribute(gg,"SCHanno",V(gg),"")

    anno    <- read.table("SCH_flatfile.csv",sep="\t",skip=1,strip.white=T)
    annoIDS <- as.character(anno[,3])

    type <- unique(unlist(strsplit(as.character(unique(anno[,2])),",")))


for( i in 1:length(ids) ){

  ind1 = which(annoIDS==ids[i])

  Str <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str <- as.character(anno[ind1[1],2]) }
      else { Str <- paste(as.character(anno[ind1,2]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$Schanno = as.character(Str);

  
}
}
#---

if( runAnno[5] ){

    gg <- removeVertexTerm(gg,"chua")    
    
    #--- Set Family attributes in .gml graph
    igraph::set.vertex.attribute(gg,"chua",V(gg),"")

    anno    <- read.table("flatfile_chua.csv",sep="\t",skip=1,strip.white=T)
    annoIDS <- as.character(anno[,3])

    type <- unique(unlist(strsplit(as.character(unique(anno[,2])),",")))

for( i in 1:length(ids) ){

  ind1 = which(annoIDS==ids[i])

  Str <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str <- as.character(anno[ind1[1],2]) }
      else { Str <- paste(as.character(anno[ind1,2]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$chua = as.character(Str);

  
}
}

if( runAnno[6] ){

    gg <- removeVertexTerm(gg,"InterProFamilyID")
    gg <- removeVertexTerm(gg,"InterProFamily")

    gg <- removeVertexTerm(gg,"InterProDomainID")
    gg <- removeVertexTerm(gg,"InterProDomain")    
    
    #--- Set interproFamily attributes in .gml graph
    igraph::set.vertex.attribute(gg,"InterProFamilyID",V(gg),"")
    igraph::set.vertex.attribute(gg,"InterProFamily",V(gg),"")
    igraph::set.vertex.attribute(gg,"InterProDomainID",V(gg),"")
    igraph::set.vertex.attribute(gg,"InterProDomain",V(gg),"")

    annoF    <- read.table("flatfile.interpro.Family.csv",sep="\t",skip=1,strip.white=T)
    annoFIDS <- as.character(annoF[,3])

    typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))

    annoD    <- read.table("flatfile.interpro.Domain.csv",sep="\t",skip=1,strip.white=T)
    annoDIDS <- as.character(annoD[,3])

    typeD <- unique(unlist(strsplit(as.character(unique(annoD[,2])),",")))

 
for( i in 1:length(ids) ){

  ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$InterProFamilyID = as.character(Str2);  
  V(gg)[i]$InterProFamily   = as.character(Str1);


  ind1 = which(annoDIDS==ids[i])

  Str1 <- "";
  Str2 <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoD[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoD[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoD[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoD[ind1,1]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$InterProDomainID = as.character(Str2);  
  V(gg)[i]$InterProDomain   = as.character(Str1);

    
}

}

#---ADD MEASURES
if( runAnno[7] ){    

    gg <- removeVertexTerm(gg,"COREPRE")
    gg <- removeVertexTerm(gg,"CORESPD")

    
    #---ADD VIP gene lists
    VIP <- vector(length=2)
    VIP[1] <- "SynsysBaits.csv"
    VIP[2] <- "CorePSD95Complex.csv"
    set1 <- read.table(sprintf("%s/%s",OUT[4],VIP[1]),sep="\t",header=F)[[1]]
    set2 <- read.table(sprintf("%s/%s",OUT[4],VIP[2]),sep="\t",header=F)[[1]]

    igraph::set.vertex.attribute(gg,"COREPRE",V(gg),"")
    for( i in 1:length(ids) ){

        ind1 = which(set1==ids[i])

        Str <- "";
        
        if( length(ind1) != 0 ){
            Str <- "YES"
        }
        
        V(gg)[i]$COREPRE = as.character(Str);

  
    }

    igraph::set.vertex.attribute(gg,"CORESPD",V(gg),"")
    for( i in 1:length(ids) ){
        
        ind1 = which(set2==ids[i])
        
        Str <- "";
        
        if( length(ind1) != 0 ){
            Str <- "YES"
        }

        V(gg)[i]$COREPSD = as.character(Str);
        
        
    }

}
    
if( runAnno[8] ){

    str <- sprintf("%s/%s/REGIONS/",OUT[4],subDIR[S])

    files <- list.files(str)

    fn <- gsub(".csv","",files)

    for( i in 1:length(files) ){

        gg <- removeVertexTerm(gg,fn[i])
        #gg <- removeVertexTerm(gg,gsub("_","",fn[i]))
        
        if( file.exists(sprintf("%s/%s",str,files[i])) ){

            ff <- read.table(sprintf("%s/%s",str,files[i]),sep="\t",header=F)

            gg <- igraph::set.vertex.attribute(gg,fn[i],V(gg),ff[match(ff[,1],ids),2])

        }
    }

}


if( runAnno[9] ){

    gg <- removeVertexTerm(gg,"GO_MF")
    gg <- removeVertexTerm(gg,"GO_MF_ID")
    
    #--- Set Disease (geneRIF db) attributes in .gml graph
    igraph::set.vertex.attribute(gg,"GO_MF",V(gg),"")
    igraph::set.vertex.attribute(gg,"GO_MF_ID",V(gg),"")

    annoF    <- read.table("flatfile.go.MF.csv",sep="\t",skip=1,strip.white=T,quote="")
    annoFIDS <- as.character(annoF[,3])

    typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))    

for( i in 1:length(ids) ){

  ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$GO_MF_ID = as.character(Str2);  
  V(gg)[i]$GO_MF    = as.character(Str1);

}
    
}

if( runAnno[10] ){

    gg <- removeVertexTerm(gg,"GO_BP")
    gg <- removeVertexTerm(gg,"GO_BP_ID")
    
    #--- Set Disease (geneRIF db) attributes in .gml graph
    igraph::set.vertex.attribute(gg,"GO_BP",V(gg),"")
    igraph::set.vertex.attribute(gg,"GO_BP_ID",V(gg),"")

    annoF    <- read.table("flatfile.go.BP.csv",sep="\t",skip=1,strip.white=T,quote="")
    annoFIDS <- as.character(annoF[,3])

    typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))    

for( i in 1:length(ids) ){

  ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$GO_BP_ID = as.character(Str2);  
  V(gg)[i]$GO_BP    = as.character(Str1);

}
    
}

if( runAnno[11] ){

    gg <- removeVertexTerm(gg,"GO_CC")
    gg <- removeVertexTerm(gg,"GO_CC_ID")
    
    #--- Set Disease (geneRIF db) attributes in .gml graph
    igraph::set.vertex.attribute(gg,"GO_CC",V(gg),"")
    igraph::set.vertex.attribute(gg,"GO_CC_ID",V(gg),"")

    annoF    <- read.table("flatfile.go.CC.csv",sep="\t",skip=1,strip.white=T,quote="")
    annoFIDS <- as.character(annoF[,3])

    typeF <- unique(unlist(strsplit(as.character(unique(annoF[,2])),",")))    

for( i in 1:length(ids) ){

  ind1 = which(annoFIDS==ids[i])

    Str1 <- "";
    Str2 <- "";
  
  if( length(ind1) != 0 ){

      if( length(ind1) == 1 ){ Str1 <- as.character(annoF[ind1[1],2]) }
      else { Str1 <- paste(as.character(annoF[ind1,2]),collapse=COLLAPSE[c]) }

      if( length(ind1) == 1 ){ Str2 <- as.character(annoF[ind1[1],1]) }
      else { Str2 <- paste(as.character(annoF[ind1,1]),collapse=COLLAPSE[c]) }
      
  }

  V(gg)[i]$GO_CC_ID = as.character(Str2);  
  V(gg)[i]$GO_CC    = as.character(Str1);

}
    
}


if( runAnno[12] ){

    files <- list.files("./")
    files <- files[grepl("celltypes_",files)]

    fn <- gsub(".csv","",files)
    fn <- gsub("celltypes_","",fn)
    fn <- sprintf("CellTypes_%s",fn)

    gg <- loopOverFiles(gg, files, fn, ids, FALSE)
    
}

if( runAnno[13] ){

    files <- list.files("./")
    files <- files[grepl("Pathways_",files)]

    fn <- gsub(".csv","",files)
    fn <- gsub("Pathways_","",fn)
    fn <- sprintf("PathWays_%s",fn)

    gg <- loopOverFiles(gg, files, fn, ids, TRUE)
    
}


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")


