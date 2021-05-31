source('../setUp.R');

#Calculate the Median absolute difference
MAD <- function( X ){

    X <- as.numeric(X)
    
    Xmd <- median(X)

    MAD <- median(abs(X-Xmd))

    return(MAD)
}

scale <- function(X, VALUE=NULL){

    X = as.numeric(as.vector(X))

    xmin <- min(X,na.rm=T)
    xmax <- max(X,na.rm=T)
   
    if( is.null(VALUE) ){

        X  <- X-xmin
        X  <- ifelse(!is.na(X), X/(xmax-xmin), NA) 

        return(X)
    }

    VALUE = as.numeric(as.vector(VALUE)[1])
    VALUE = VALUE-xmin
    VALUE = ifelse(!is.na(VALUE), VALUE/(xmax-xmin), NA) 
    return(VALUE)
}

Counts <- function(x){

    if( x=="" ){ return(0) }
    else {
        return( str_count(x,";")+1 )
    }
}

#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("POWERlawFIT",DIRS)]
OUT[3] <- DIRS[grepl("SVI",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

log.adj <- function( X ){

    X = as.numeric(X)
    if( X == 0 ){ return(0) }
    else        { return(-X*log(X)) }
    
}

#---Set Options
runBridge  <- vector(length=2)
runBridge[1] <- 0  #Calculate Bridgeness
runBridge[2] <- 1  #Plot Bridgeness


#---declare clustering algorithms in graph, and with a corresponding consensus matrix
set <- c("Spectral")
alg <- ALGS[match(set,ALGS)]

#CorePSD95Complex
corePSD <- read.delim("CorePSD95Complex.csv", sep="\t", header=F)[[1]]
corePSD <- unique(corePSD)


#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

if(runBridge[1]){

    N    <- length(V(gg)$name)
    
    mm   <- readRDS(sprintf("%s/originalNtwrkMeasures.rds",OUT[2]))

    INDX <- which(names(mm)==subDIR[S])

    INDR <- match(V(gg)$name,mm[[INDX]][,1])
    
    CN  <- c('ENTREZ.ID','GENE.NAME','DEGREE','Closeness','Bet','CC','Cl','Clnorm','SP','PR',alg,sprintf("BRIDGE_%s",alg),sprintf("PROB_%s",alg),sprintf("MIX_%s",alg),'BR_CONSENSUS','BR_CONSENSUS_MAD',"BR_CONSENSUS_ADJ","CorePSD95")

    COREIndx = which(CN=="CorePSD95")
    
    meas     <- matrix(0, nrow=N, ncol=length(CN))
    colnames(meas) <- CN

    meas[,1] <- as.character(V(gg)$name)
    meas[,2] <- as.character(V(gg)$GeneName)

    meas[,3]  <- as.character(mm[[INDX]][INDR,2])
    meas[,4]  <- as.numeric(round(closeness(gg,mode="all",normalized=T),3))
    meas[,5]  <- as.numeric(mm[[INDX]][INDR,3])
    meas[,6]  <- as.numeric(mm[[INDX]][INDR,4])
    meas[,7]  <- as.numeric(mm[[INDX]][INDR,5])
    meas[,8]  <- as.character( (as.numeric(meas[,7]) - min(as.numeric(meas[,7])))/(max(as.numeric(meas[,7])) - min(as.numeric(meas[,7]))) )
    meas[,9]  <- as.numeric(mm[[INDX]][INDR,6])
    meas[,10] <- as.numeric(mm[[INDX]][INDR,7])

    meas[match(corePSD,meas[,1]),COREIndx] = 1

      
#START filling meas after PageRank column
FROM <- which(CN=="PR")

N <- length(V(gg))
M <- length(E(gg))

## run over each algorithm
    for( a in 1:length(alg) ){

        cat("calculating Bridgeness for: ", alg[a], "\n")
    
        ##--- build reference matrix from the graph
        refin     <- matrix("",ncol=3,nrow=length(V(gg)$name))
        refin[,1] <- igraph::get.vertex.attribute(gg,"name",V(gg))
        refin[,2] <- igraph::get.vertex.attribute(gg,"name",V(gg))
        refin[,3] <- igraph::get.vertex.attribute(gg,alg[a],V(gg))

        ##--- store alg. cluster results
        meas[,(FROM+a)] <- as.numeric(refin[,3])

        ##--- path to consensus matrix
        st1 = sprintf("%s/%s/%s/consensusmatrix.txt.gz",rndDIR[1],subDIR[S],alg[a]);
        
        ##Read in consensus matrix
        filein = read.table(gzfile(st1), header=FALSE, sep=",");
        dimnames(filein)[2] <- dimnames(filein)[1]
        
        ##format reference matrix
        ##the reference matrix with the correct row.names(as a data.frame)
        refmat = as.matrix(refin);
        ref    = refmat[,2:3];
        rm           <- data.frame(ref);
        rownames(rm) <- rm$X1;
        rm$X1        <- NULL;
        names(rm)    <- 'cm';

        ##format consensus matrix
        ##the consensus matrix you may have made (as a numeric matrix) 
        conmat = as.matrix(filein);
        cm           <- data.frame(conmat);
        names(cm)    <- rownames(rm);
        rownames(cm) <- rownames(rm);
        cm           <- as.matrix(cm);

        rm(filein,refin)

        ##get consensus matrix indices for each edge in edge list
        indA <- match(igraph::get.edgelist(gg)[,1],rownames(cm))
        indB <- match(igraph::get.edgelist(gg)[,2],rownames(cm))

        dat  <- data.frame(indA,indB)
    
        ##get community numbers for each vertex in edge list for the algorithm
        elA    <- igraph::get.vertex.attribute(gg,alg[a],V(gg))[match(igraph::get.edgelist(gg)[,1],V(gg)$name)]
        elB    <- igraph::get.vertex.attribute(gg,alg[a],V(gg))[match(igraph::get.edgelist(gg)[,2],V(gg)$name)]
    
        ed      <- matrix(ncol=6,nrow=length(E(gg)))
        ed[,1]  <- igraph::get.edgelist(gg)[,1]
        ed[,2]  <- igraph::get.edgelist(gg)[,2]
        ed[,3]  <- elA
        ed[,4]  <- elB
        ed[,5]  <- apply(dat,1,function(x,mat) mat[x[1],x[2]], mat=cm)
        ed[,6]  <- (elA-elB)
    
        Cmax  <- max(igraph::get.vertex.attribute(gg,alg[a],V(gg)))
    
        rm(cm)

        coms <- matrix(ncol=5,nrow=Cmax)

        for( m in 1:Cmax ){        

            Mtot = length(which(elA==m | elB==m))
            
            Min  = length(which(elA==m & elB==m))

            Mout = Mtot - Min

            Cn   = length(which(igraph::get.vertex.attribute(gg,alg[a],V(gg))==m ))

            Kw=max(as.numeric(meas[,3]))
        
            coms[m,1] <- m
            coms[m,2] <- 2*Min
            coms[m,3] <- Mout
            coms[m,4] <- Cn
            coms[m,5] <- Kw
        
        }

        ##store the probabilities of gene belonging to each community
        Vprobs           <- matrix(0,nrow=N,ncol=(2+Cmax))
        CNv              <- c("ENTREZ.ID","GENE.NAME",seq(1,Cmax,1))
        colnames(Vprobs) <- CNv    
            
        ##loop over each vertex in the graph
        for( i in 1:length(V(gg)) ){

            ##get edges belonging to the i'th veretx
            ind <- which(ed[,1] == V(gg)$name[i] | ed[,2] == V(gg)$name[i])

            ##get community belonging to the i'th vertex       
            c <- igraph::get.vertex.attribute(gg,alg[a],V(gg))[i]

            ##reorder edge communities, so ed[,3] equals current community no: 'c'
            for( k in 1:length(ind) ){
                if( ed[ind[k],6] != 0 && ed[ind[k],4] == c ){
                    ed[ind[k],4] <- ed[ind[k],3]
                    ed[ind[k],3] <- c
                }
            }
        
            ##prob of i'th vertex being found in community Cn relative to a random model
            ##Lancichinetti et al, Statistical significance of communities in networks (2010). Physics Review E, 81, 046110.  
            K   <- length(ind)
            Kin <- 0
            indX <- which(ed[ind,3] == c & ed[ind,4] == c)
            if( length(indX) != 0 ){
                Kin <- length(indX)
            }
            Kout <- (K - Kin)

            ##eqn 1 in (Lancichinetti et al, 2010)    
            t1 <- choose( (coms[c,3] - Kout), Kin )
            t2 <- choose( (2*M - sum(coms[c,2:3]) - K - coms[c,3] -Kout), (K-Kin) )
            t3 <- choose( (2*M - sum(coms[c,2:3]) - K), K)

            prob <- (t1*t2)/t3

            if( is.na(prob) ){ prob <- 0 }

            ##PROB_
            meas[i,(FROM+2*length(alg)+a)] <- (1-prob)

            ##MIX_
            meas[i,(FROM+3*length(alg)+a)] <- Kout/K

            ##number of communities i'th vertex is connected too (via it's edges)
            cc <- unique(ed[ind,4])

            ##use sum of consensus values to calculate the likelihood of i'th
            ##vertex beloning to to k'th community. 
            prob <- vector(length=length(cc))
            for( k in 1:length(cc) ){                
                prob[k] = sum(as.numeric(ed[which(ed[ind,4]==cc[k]),5]))/length(ind)
            }

            ##normalise
            prob <- prob/sum(prob)
        
            ##calculate bridgeness of i'th vertex
            ##Fuzzy communities and the concept of bridgeness in complex networks, T. Nepusz, arXiv, 2007
            b    <- sum( (prob - 1/Cmax) * (prob - 1/Cmax))

            Kzero <- Cmax - length(cc)
            b = b + sum(rep((1/(Cmax*Cmax)),times=Kzero))
        
            ##store values
            ##BRIDGE_
            meas[i,(FROM+length(alg)+a)]  <- 1-sqrt( Cmax/(Cmax-1) * b )

            ##store vertex probs across communities
            INDC <- match(ed[ind,4],colnames(Vprobs))
            for( k in 1:length(INDC) ){
                Vprobs[i,INDC[k]] <- (as.numeric(ed[ind[k],5]) + as.numeric(Vprobs[i,INDC[k]]));
            }
        
        }#inner for

        ##normalise vertex probs per cluster
        for(i in 1:length(V(gg)) ){
            Vprobs[i,3:length(colnames(Vprobs))] <- as.vector(Vprobs[i,3:length(colnames(Vprobs))])/sum(as.vector(Vprobs[i,3:length(colnames(Vprobs))]))
        }        
        
        Vprobs[,1] <- as.character(V(gg)$name)
        Vprobs[,2] <- as.character(V(gg)$GeneName)

        ##print vertex probs to file
        write.table(Vprobs,sprintf("%s/%s_Vprobs.csv",subDIR[S],alg[a]),row.names=F,col.names=T,sep="\t",quote=F)
    
    }#for
 

    
    outfile <- file(sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]),"w")
    write.table(meas, file=outfile, append=T, row.names=F, col.names=T, sep="\t",quote=F);
    close(outfile);

}


if( runBridge[2] ){

    ##Read-in all measures    
    meas <- read.table(sprintf("%s/%s_Measures.csv",subDIR[S],subDIR[S]),sep="\t",header=T)

    Xmeas = "Cl"
    #Xmeas = "entropy"

    if( Xmeas == "entropy" ){
    
        ##--- Read-in Entropy values
        ent <- read.table(sprintf("../EntropyRate/%s/SignalEntropyRate.csv",subDIR[S]),sep="\t",header=T)

        X=ent[,4]
        Xl = X
        SRmax = 0.74024#---Initial Entropy value
        IntCon=scale(X,VALUE=SRmax)
        X=scale(X)
        xtck = seq(min(X), max(X), length.out=9)
        xlab = round(seq(min(Xl),max(Xl),length.out=length(xtck)),3)
        ytck = seq(0,1,length.out=9)
        ylab = seq(0,1,length.out=9)
        Xlab = "Entropy"
        Ylab = "Bridgeness (B)"

        MainDivSize=0.8
        xmin <- 0
        xmax <- 1
        ymin <- 0
        ymax <- 1

        #---Set Scheme in Regions.R for B V. Centrality Plot Regions
        Scheme = 3
        source(sprintf("%s/Regions.R",pramFILES))
        
    }

    if( Xmeas == "Cl" ){
         indB <- which(colnames(meas)=="Cl")
         X    <- as.numeric(as.vector(meas[,indB]))
         X    <- scale(X)
         Xlab <- "Local Centrality (Cl)" 
         Ylab <- "Bridgeness (B)"

         MainDivSize=0.8
         xmin <- 0
         xmax <- 1
         ymin <- 0
         ymax <- 1

         #---Set Scheme in Regions.R for B V. Centrality Plot Regions
         Scheme = 1
         source(sprintf("%s/Regions.R",pramFILES))
         
    }

    
    ##Build Consensus Bridgness 
    indC   <- match(sprintf("BRIDGE_%s",alg),colnames(meas))
    if( length(indC) == 1 ){
        val    <- median(meas[,indC])
        mad    <- MAD( meas[,indC])
        valAdj <- (val-meas[,indC])/mad
    } else {
        val    <- apply(meas[,indC],1,median)
        mad    <- apply(meas[,indC],1,MAD)

        XX     <- as.vector(unlist(meas[,indC]))
        XXmed  <- median(XX)
        XXmad  <- MAD(XX)
        valAdj <- (val-XXmed)/XXmad
    }
        

    
    ##---
    ## Plot Bridging (B) V. Semi-local Centrality (SL)
    ##---
    
    ##---Check or create output plot dir
    plotDIR <- sprintf("%s/PLOTS/",subDIR[S])

    if( !file_test("-d",plotDIR) ){
        dir.create(plotDIR)
    }

    ##---Check or create output REGIONS dir
    regDIR <- sprintf("%s/REGIONS/",subDIR[S])

    if( !file_test("-d",regDIR) ){
        dir.create(regDIR)
    }

    

    ##--- set rndm seed and for geom_text_repel
    set.seed(42)

    Force1    <- vector(length=3)
    Force1[1] <- 10
    Force1[2] <- 5
    Force1[3] <- 10

    Force3    <- vector(length=3)
    Force3[1] <- 10
    Force3[2] <- 1
    Force3[3] <- 10
    ##----    
    
    cons <- matrix(0,ncol=3,nrow=length(meas[,1]))
    colnames(cons) <- c(colnames(meas)[c(1,2)],"reg")
    cons[,1] <- meas[,1]              #Entrez ID
    cons[,2] <- as.character(meas[,2])#Gene Symbol
    cons[,3] <- rep(0,length(meas[,1]))
    
    alg <- c(alg)

    ##---Read-in Core PSD95 interactor genes
    CORE  <- read.table("CorePSD95Complex.csv",sep="\t",header=F)[[1]]
    indx   = match(V(gg)$name,CORE)
    group  = ifelse( is.na(indx), 0,1)
    
    ##---Read-in my gene list
    SPEC   <- read.table("BrProteins_2020.csv",sep="\t",header=F)[[1]]
    indx    = match(V(gg)$name,SPEC)
    group2  = ifelse( is.na(indx), 0,1)

    ##Also set these as VIPs
    VIPs = SPEC    
    ##---    

    
    ##---B V. SL Plot Regions
    for( a in 1:length(alg) ){

        dd <- data.frame(meas)
    
        str <- sprintf("BRIDGE_%s",alg[a])
   
        indA    <- which(colnames(dd)==str)
        bridge  <- dd[,indA]

        ## get Bridgeness for our algorithm
        Y <- as.numeric(as.vector(bridge))    
    
        cat("Plotting Br V. ", Xmeas ," for: ", alg[a], "\n")

        ##---Quadrant (or Region) each gene lies in
        cons[,3] <- rep(0,length(cons[,1])) 
        cons     <- fillRegions( cons, X, Y, 3 );      
 
    
        ##---store Br. V. SL Region data for alg.
        oo <- data.frame(cons[,c(1,3)])
        outfile <- file(sprintf("%s/%s_REGION.csv",regDIR,alg[a]),"w")
        cat("#region",file=outfile,"\n")
        write.table(oo, file=outfile, append=T, row.names=F, col.names=F, sep="\t",quote=F);
        close(outfile);
    
        genes1 <- ifelse( cons[,3]==quad[1],cons[,2],"")
        genes2 <- ifelse( cons[,3]==quad[2],cons[,2],"")
        genes3 <- ifelse( cons[,3]==quad[3],cons[,2],"")
        genes4 <- ifelse( cons[,3]==quad[4],cons[,2],"")
        genes5 <- ifelse( cons[,3]==quad[5],cons[,2],"")
        genes6 <- ifelse( cons[,3]==quad[6],cons[,2],"")

        VIPsGN <- cons[match(VIPs,cons[,1]),2]   
      
        genes1 <- ifelse(!is.na(match(genes1,VIPsGN)),genes1,"")
        ##genes2 <- ifelse(!is.na(match(genes2,VIPsGN)),genes2,"")
        genes3 <- ifelse(!is.na(match(genes3,VIPsGN)),genes3,"")
        genes4 <- ifelse(!is.na(match(genes4,VIPsGN)),genes4,"")
        genes5 <- ifelse(!is.na(match(genes5,VIPsGN)),genes5,"")
        genes6 <- ifelse(!is.na(match(genes6,VIPsGN)),genes6,"")        

        GenesUL <- genes1
        GenesUR <- genes2
        GenesLL <- genes3
        GenesLR <- genes4

        if( Scheme == 2 ){            
            GenesUL <- genes4
            GenesUR <- genes6
            GenesLL <- genes2
            GenesLR <- genes3
        }


   
        if( Xmeas == "entropy" ) {##Br V entropy plot 

           df = cbind(X,Y,group,group2)
            df = as.data.frame(df)
            df_core = df[df[,3]==1,]
            df_sp   = df[df[,4]==1 & df[,3] !=1,]
            df_bs   = df[df[,3] != 1 & df[,4] != 1,]
        
            ##plot
            gplot <- ggplot(df,aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)) ))+
                ##base
                geom_point(data=df_bs,
                           aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)), alpha=(X*Y)),colour="magenta",shape=16,show.legend=F)+
                ##Special
                geom_point(data=df_sp,
                           aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y))),colour="magenta",alpha=1,shape=16,show.legend=F)+
                ##core PSD95
                geom_point(data=df_core,
                           aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y))),colour="deepskyblue3",alpha=1,shape=16,show.legend=F)+
            
    ##lr
    ##geom_text_repel(aes(label=as.vector(GenesLR)),force=5, xlim=c(0.5,1.0),ylim=c(0.0,0.5),color="black",point.padding=NA,size=rel(3.8),show.legend=F)+ 

    ##ur
    ##geom_text_repel(aes(label=as.vector(GenesUR)),color='black',fontface='bold',segment.color='grey50',xlim = c(0.3,1.0), ylim = c(0.3,1.0),point.padding=NA,force=5,size=rel(3.8),show.legend=F)+  

    ##prim-loc-2
        geom_label_repel(aes(label=as.vector(GenesLL)),force=5, xlim=c(0.0,0.5),ylim=c(0.1,0.5),color="black",fontface="bold",segment.color='grey50',fill='white',box.padding=0.0,point.padding=NA,label.padding=0.1,size=rel(3.0),show.legend=F)+

    ##ul
        geom_label_repel(aes(label=as.vector(GenesUL)),fontface='bold',color='black',fill='white',box.padding=0.1,point.padding=NA,label.padding=0.2,segment.color='grey50',xlim = c(0.0,0.6), ylim = c(0.5, 1.0),force=10,size=rel(3.8),show.legend=F)+
    
        labs(x=Xlab,y=Ylab,title=sprintf("%s",alg[a]))+
            coord_cartesian(xlim = c(xmin,xmax), ylim=c(ymin,ymax))+
                ##xlim(c(0,1))+ylim(c(0,1))+
                scale_y_discrete(expand=c(0,0), limit=ytck, labels=ylab)+
                scale_x_discrete(expand=c(0,0), limit=xtck, labels=xlab)+
                theme(            
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(1.5)),
                    legend.text=element_text(face="bold",size=rel(1.5)),
                    legend.key=element_blank())+
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+

    geom_vline(xintercept=IntCon,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
        geom_hline(yintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
            
                Gseg1 +
                Gseg2 +
                Gseg3 +
                GAnno1 +
                GAnno2 +
                GAnno3 +
                GAnno4 +
                GAnno5 +      
    
       png(sprintf("%s/%s_%s_BrV%s.png",plotDIR,subDIR[S],alg[a],Xmeas),width=WIDTH,height=HEIGHT,units="px")
       print(gplot)
       dev.off()

    }

        if( Xmeas == "Cl" ) {##Br V Cl plot 

            df = cbind(X,Y,group,group2)
            df = as.data.frame(df)
            df_core = df[df[,3]==1,]
            df_sp   = df[df[,4]==1 & df[,3] !=1,]
            df_bs   = df[df[,3] != 1 & df[,4] != 1,]       

            ##baseColor="magenta"
            baseColor="royalblue2"

            ##PSDColor="royalblue2"
            PSDColor="magenta"
        
            ##plot
            gplot <- ggplot(df,aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)) ))+          
                
                ##base
                geom_point(data=df_bs,
                           aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y)), alpha=(X*Y)),colour=baseColor,shape=16,show.legend=F)+
                ##Special
                geom_point(data=df_sp,
                           aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y))),colour=baseColor,alpha=1,shape=16,show.legend=F)+
                ##core PSD95
                geom_point(data=df_core,
                           aes(x=as.numeric(as.vector(X)),y=as.numeric(as.vector(Y))),colour=PSDColor,alpha=1,shape=16,show.legend=F)+
            
           
                ##ul
                geom_label_repel(aes(label=as.vector(GenesUL)),fontface='bold',color='black',fill='white',box.padding=0.1,point.padding=NA,label.padding=0.15,segment.color='black',xlim = c(0.0,0.5), ylim = c(0.5, 1.0),force=1,size=rel(3.8),show.legend=F)+
   
                ##ur
                geom_text_repel(aes(label=as.vector(GenesUR)),color='black',fontface='bold',segment.color='grey50',xlim = c(0.45,1.0), ylim = c(0.45,1.0),point.padding=NA,force=0.5,size=rel(3.8),show.legend=F)+  

                ##prim-loc-2
                geom_label_repel(aes(label=as.vector(GenesLL)),force=0.5, xlim=c(0.0,0.5),ylim=c(0.1,0.5),color="black",fontface="bold",segment.color='black',fill='white',box.padding=0.0,point.padding=NA,label.padding=0.1,size=rel(3.0),show.legend=F)+

                ##lr
                geom_text_repel(aes(label=as.vector(GenesLR)),force=5, xlim=c(0.5,1.0),ylim=c(0.0,0.5),color="black",point.padding=NA,size=rel(3.8),show.legend=F)+ 
    
        labs(x=Xlab,y=Ylab,title=sprintf("%s",alg[a]))+
            scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) + 
                scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax))+
                theme(            
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(1.5)),
                    legend.text=element_text(face="bold",size=rel(1.5)),
                    legend.key=element_blank())+
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+

        geom_vline(xintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
        geom_hline(yintercept=0.5,colour="grey40",size=MainDivSize,linetype=2,show.legend=F)+
            
        Gseg1 +
        Gseg2 +
        Gseg3 +

        GAnno1 +
        GAnno2 +
        GAnno3 +
        GAnno4 +
        GAnno5 +

       ##svg(sprintf("%s/%s_%s_BrV%s.svg",plotDIR,subDIR[S],alg[a],Xmeas),width=WIDTHin,height=HEIGHTin)
       ##print(gplot)
       ##dev.off()
    
       png(sprintf("%s/%s_%s_BrV%s.png",plotDIR,subDIR[S],alg[a],Xmeas),width=WIDTH,height=HEIGHT,units="px")
       print(gplot)
       dev.off()

    }



     
    
}

}


    
