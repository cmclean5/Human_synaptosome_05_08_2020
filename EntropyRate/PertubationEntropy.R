# PPI network Permutation analysis
# Based on Figure 4C in: A. E. Teschendorff et al. Increased signaling entropy in cancer requires the scale-free property of protein interaction networks, Scientific Reports, 5:9646.

source('../setUp.R')

#---Directories needed
OUT <- vector(length=3)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("EntropyRate",DIRS)]

#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
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


#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

# ER random graph
#rnd <- igraph::sample_gnm( n=length(V(gg)), m=length(E(gg)), directed=F, loops=F)
#rnd <- igraph::set.vertex.attribute(rnd, "name",     V(rnd), V(gg)$name)
#rnd <- igraph::set.vertex.attribute(rnd, "GeneName", V(rnd), V(gg)$GeneName)
#gg  <- findLCC(rnd)
#---

#--- initial entropy rate
V    <- length(V(gg))
E    <- length(E(gg))
ki   <- as.vector(igraph::degree(gg))
Kbar <- mean(ki)

#--- get adjacency matrix for graph
A    <- get.adjacency(gg)

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


#--- perturbated each PPI node/gene

#--- expression values
xx    <- vector(length=2)
xx[1] <- 2  #active
xx[2] <- 16 #inactive

#--- perturbation expression values
lambda    <- vector(length=2)
lambda[1] <- 14   #active
lambda[2] <- -14  #inactive

#--- Norm for PI'
NORM      <- vector(length=2)
NORM[1]   <- 0
NORM[2]   <- 0

SRprime <- cbind(V(gg)$name, V(gg)$GeneName, ki, rep("",V), rep("",V))

for( v in 1:V ){

    #--- name of gene to perturb
    GN     <- as.character(SRprime[v,2])
    GNindx <- which(V(gg)$GeneName==GN)

    #--- PI'
    PIprime <- cbind( rep("",V), rep("",V) )

    #--- LS'
    LSprime <- cbind( rep("",V), rep("",V) )

    #--- reset NORM
    NORM[1] = 0; NORM[2] = 0;
    
    #--- calculate norm for PI'
    for( s in 1:length(lambda) ){
        X               <- rep(xx[s], V)
        X[GNindx[1]]    <- X[GNindx[1]] + lambda[s]
        NORM[s]         <- X %*% A %*% X
    }
    

    #--- find all neighors to v, i.e. N(v)
    Nv <- V(gg)$name[neighbors(gg,GNindx,mode="all")]

    oo <- cbind( ki, !(V(gg)$name %in% Nv) )
    

    #--- PI' when v is not N(v)
    for( s in 1:length(lambda) ){
        PIprime[,s] <- ifelse(oo[,2] == 1, (1/NORM[s] * xx[s] * xx[s] * as.numeric(oo[,1])), ".")
    }

    #--- PI' when v is v
    for( s in 1:length(lambda) ){

        X   <- as.numeric(xx[s])
        lam <- as.numeric(lambda[s])
        DEG <- as.numeric(oo[GNindx[1],1])

        PIprime[GNindx[1],s] <- ((X + lam) * DEG * X) / NORM[s] 
        
    }
    
    #--- PI' when v is N(v) 
    for( s in 1:length(lambda) ){
        PIprime[,s] <- ifelse(oo[,2] == 0, (1/NORM[s] * xx[s] * ( xx[s] + lambda[s] + (as.numeric(oo[,1]) - 1) * xx[s])),PIprime[,s])
    }


    #--- LS' when v is not N(v)
    for( s in 1:length(lambda) ){
        X <- as.numeric(xx[s])
        LSprime[,s] <- ifelse(oo[,2] == 1, (-log(X) + log(X*as.numeric(oo[,1]))),".")
    }

    #--- LS' when v is N(v)     
    Ni <- grep(0,oo[,2])

    for( i in 1:length(Ni) ){

        DEGi <- as.numeric(oo[Ni[i],1])
        SUM  <- DEGi-1     
        
        for( s in 1:length(lambda) ){

            X   <- as.numeric(xx[s])
            lam <- as.numeric(lambda[s])
                        
            dem <- X + lam + (DEGi -1) * X
            
            pij <- X / dem
            pi1 <- (X + lam) / dem

            #cat("v",v, ": i ", i, ": pij ", pij, ": pi1 ", pi1,"\n")
            
            LSi <- pij * log(pij)
            LSi <- - SUM * LSi - pi1 * log(pi1)
                    
            LSprime[Ni[i],s]  <- as.character(LSi)
    
        }

    }
    
    SRprime[v,4] <- sum( as.numeric(PIprime[,1]) * as.numeric(LSprime[,1]) )
    SRprime[v,5] <- sum( as.numeric(PIprime[,2]) * as.numeric(LSprime[,2]) )

}


#--- PLOT
subTIT    <- rep("",length=S)
#subTIT    <- vector(length=3)
#subTIT[1] <- ""#"Presynaptic"
#subTIT[2] <- ""#"PSD"
#subTIT[3] <- ""#"PSD"

colours <- c('lawngreen','firebrick2')

SRprime[,4] <- as.numeric(SRprime[,4])/maxSr
SRprime[,5] <- as.numeric(SRprime[,5])/maxSr

colnames(SRprime) <- c("ENTREZ.ID","GENE.NAME","DEGREE","UP","DOWN")
#colnames(SRprime) <- c("ENTREZ.ID","GENE.NAME","DEGREE","SR_UP","SR_DOWN")

write.table(SRprime,sprintf("%s/SignalEntropyRate.csv",subDIR[S]),sep="\t",row.names=F,col.names=T,quote=F)

IntCon = SRo/maxSr
outfile <- file(sprintf("%s/%s_initial_configuration.csv",subDIR[S],subDIR[S]),"w")
cat(sprintf("%.5f",IntCon), file=outfile,"\n")
close(outfile);


DF1 <- SRprime[,c(1,2,3,4)]
DF1 <- cbind(DF1,rep("SR_UP",length(SRprime[,1])))

DF2 <- SRprime[,c(1,2,3,5)]
DF2 <- cbind(DF2,rep("SR_DOWN",length(SRprime[,1])))

#--- Bottom 1% UP, i.e. OVER-EXPRESSED
XX  <- as.numeric(SRprime[,4])
MIN <- min(XX)
MAX <- max(XX)
XX2 <- (XX-MIN)/(MAX-MIN)
oo2 <- cbind(SRprime[,2],XX2)
oo2 <- oo2[order(as.numeric(oo2[,2])),]
ii  <- floor(1/100 * V)
GN  <- oo2[1:ii,1]
DF3 <- SRprime[match(GN,SRprime[,2]),c(1,2,3,4)]
DF3 <- cbind(DF3,rep("1%",length(GN)))

write.table(DF3,sprintf("%s/SignalEntropyRate_1percent_UP.csv",subDIR[S]),sep="\t",row.names=F,col.names=T,quote=F)
#---

DF  <- rbind(DF1,DF2)
colnames(DF) <- c("ENTREZ.ID","GENE.NAME","DEGREE","SR","GROUP")

DF <- as.data.frame(DF)

gplot <- ggplot(DF,aes(x=log(as.numeric(as.vector(DF$DEGREE))),y=as.numeric(as.vector(DF$SR)), colour=DF$GROUP) )+
    geom_point()+
    labs(x="log(k)",y="SR",title=subTIT[S])+
    guides(color=guide_legend(override.aes=list(fill=NA,size=4)),
           fill  = FALSE,
           group = FALSE,
           alpha = FALSE)+
    theme(            
        axis.title.x=element_text(face="bold",size=rel(2)),
        axis.text.x =element_text(face="bold",size=rel(2)), 
        axis.title.y=element_text(face="bold",size=rel(2)),
        axis.text.y =element_text(face="bold",size=rel(2)), 
        legend.title=element_text(face="bold",size=rel(1.5)),
        legend.text=element_text(face="bold",size=rel(1.5)),
        legend.position="top",
        legend.key=element_blank())+    
    theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
          panel.grid.minor = element_line(colour="grey40",size=0.1),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(linetype="solid",fill=NA))+
    scale_color_manual("",breaks=levels(factor(DF$GROUP)),values=c(colours))+
    geom_hline(yintercept=SRo/maxSr,colour="black",size=2,linetype=2,show.legend=F)
   #geom_hline(yintercept=SRo/maxSr,colour="grey40",size=2,linetype=2,show.legend=F)

png(sprintf("%s/SignalEntropyRate.png",subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot)
dev.off()


## test

run=0

if( run ){

    OBS=0.66801

    ERmn=0.9890
    ERsd=0.0032#0.0005
        
    PLobsmn=0.9127
    PLobssd=0.0032
    
    PL2mn=0.312 #(colour == turquoise3)
    PL2sd=0.0036    

    FillColour <- vector(length=3)
    FillColour[1] <- "black"
    FillColour[2] <- "tan2"
    FillColour[3] <- "steelblue3"    
    #FillColour[4] <- "turquoise3"

    #names(FillColour) = c("A","B","C","D")
    
    
    Leg <- vector(length=3)
    Leg[1] <- c(expression(paste(PSP)))
    Leg[2] <- c(expression(paste("E-R")))
    Leg[3] <- c(expression(paste("P-L")))
    #Leg[1] <- c(expression(paste(PPI[obs])))
    #Leg[2] <- c(expression(paste("E-R")))
    #Leg[3] <- c(expression(paste(alpha,"=",PPI[obs])))
    ##Leg[4] <- c(expression(paste(alpha,"=",2)))

    #names(Leg) = c("A","B","C","D")
    

    Alpha <- vector(length=3)
    Alpha[1] <- 1
    Alpha[2] <- 1
    Alpha[3] <- 1
    #Alpha[4] <- 1
    
    ## x-range log(max(ki) = log(232) = 5.5    
    y=seq(0,1,0.001)
    x=seq(0,5.5, 5.5/length(y))
    x=x[1:length(y)]
    d = data.frame(x=x,y=y)

    N=length(x)
    d2 = cbind(x=x,y=rep(OBS,length=N),group=rep("A",length=N),up=rep(0,length=N),low=rep(0,length=N))
    d3 = cbind(x=x,y=rep(ERmn,length=N),group=rep("B",length=N),up=rep((ERmn+ERsd),length=N), low=rep((ERmn-ERsd),length=N) )
    d4 = cbind(x=x,y=rep(PLobsmn,length=N),group=rep("C",length=N),up=rep((PLobsmn+PLobssd),length=N), low=rep((PLobsmn-PLobssd),length=N))
    d5 = cbind(x=x,y=rep(PL2mn,length=N),group=rep("D",length=N),up=rep((PL2mn+PL2sd),length=N), low=rep((PL2mn-PL2sd),length=N))

    ##df = rbind(d2,d3,d4,d5)
    df = rbind(d2,d3,d4)
    df = as.data.frame(df)

    #The only way i can make this work is swapping the "C" in df4 to "D" and vice versa.
    
    #df$group <- factor(df$group, levels=c("A","B","C","D"), labels=c("A","B","C","D"))

    gplot <-  ggplot(df, aes(x=round(as.numeric(as.vector(x)),3), y=round(as.numeric(as.vector(y)),3), colour=group, fill=group))+
        geom_line(size=1)+
        scale_colour_manual(name = "", values=FillColour, breaks=levels(factor(df$group)), labels=Leg )+
        guides(color=guide_legend(override.aes=list(fill=NA,size=3)),
               fill  = FALSE,
               group = FALSE,
               alpha = FALSE)+
        labs(x="log(k)",y="SR",title="")+        
        theme(            
            axis.title.x=element_text(face="bold",size=rel(2)),
            axis.text.x =element_text(face="bold",size=rel(2)), 
            axis.title.y=element_text(face="bold",size=rel(2)),
            axis.text.y =element_text(face="bold",size=rel(2)), 
            legend.title=element_blank(),#element_text(face="bold",size=rel(1.5)),
            legend.text=element_text(face="bold",size=rel(2)),
            legend.text.align=0,
            legend.justification=c(1,0),
            legend.position=c(0.95,0.01))+
        #geom_ribbon(aes(ymin=round(as.numeric(as.vector(low)),3), ymax=round(as.numeric(as.vector(up)),3)), alpha=.3, linetype=0)+
            
    theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(x))) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
        #coord_cartesian(xlim = c(min(x), max(x)), ylim=c(0, 1))
       
    
    png(sprintf("%s/graphEntropies.png",subDIR[S]),width=floor(WIDTH/2),height=HEIGHT,units="px");
print(gplot)
dev.off()
  
}
