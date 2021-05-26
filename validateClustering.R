source('setUp.R')
require(cowplot)

#
setClass(Class="VAL",representation(
                         RI="numeric",
                         ARI="numeric",
                         NMI="numeric"
                     )
         )			   

#inverse community frequency
icf <- function( N, GG, ALGS ){

    XX = get.vertex.attribute(GG,ALGS)
    XX = table(XX)

    MM = matrix(NA,ncol=3,nrow=length(XX))
    MM[,1] = names(XX)
    MM[,2] = as.vector(XX)
    MM[,3] = log( (as.numeric(N))/(as.numeric(MM[,2])) )

    return(MM)
    
}

icfPlots <- function(V1, DIR, subDIR, algNAME, shortTIT, MIN, MAX, BINS=30 ){

    DF = as.data.frame(V1)
        
    gplot1 <- ggplot(DF)+geom_histogram(bins=BINS,aes(x=DF$V1),colour="brown1",fill="cadetblue") +
        labs(x=c("ICF"),y="Count",title=sprintf("%s",algNAME))+
            xlim(c(MIN,MAX))+
            #scale_x_continuous(limit = c(MIN, MAX))+
        theme(
            title=element_text(face="bold",size=rel(1.5)),
            axis.title.x=element_text(face="bold",size=rel(2.0)),
            axis.title.y=element_text(face="bold",size=rel(2.0)),
            legend.title=element_text(face="bold",size=rel(2.0)),
            legend.text=element_text(face="bold",size=rel(2.0)),
            legend.key=element_blank())+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))
                
        png(sprintf("%s/%s_%s_%s.png",DIR, shortTIT, algNAME, subDIR),width=WIDTH,height=HEIGHT,units="px")
        print(gplot1)
        dev.off()

    return(gplot1)

}

compare <- function( RES1, RES2, VAR="sqrt" ){

    val1 = RI(RES1,RES2)
    val2 = ARI(RES1,RES2)
    val3 = NMI(RES1,RES2,variant=VAR)
    
    return(new("VAL",
               RI=val1,
               ARI=val2,
               NMI=val3))
}

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }

  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

plotMetric <- function(OO, DIR, subDIR, TIT, shortTIT, valSIZE=3){

    cormat    <- round(cor(OO),2)
    cormat    <- reorder_cormat(cormat)
    upper_tri <- get_upper_tri(cormat)

    write.table( cormat, file=sprintf("%s/%s_cor_%s.csv",DIR,subDIR,shortTIT), append=F, row.names=T, col.names=T, sep="\t", quote=F)
    
    # Melt the correlation matrix
    melted_cormat <- reshape2::melt(upper_tri, na.rm=T)
    df = melted_cormat    
    colnames(df) <- c("V1","V2","value")
    
    # Create a ggheatmap
    gplot1 <- ggplot(df, aes(df$V2, df$V1, fill = df$value))+
        labs(x="",y="",title=TIT)+
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, limit = c(-1,1), space = "Lab", 
                             name="Pearson\nCorrelation") +
        theme_minimal()+ # minimal theme
        theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                         size = 12, hjust = 1))+
        coord_fixed()+
        geom_text(aes(df$V2, df$V1, label = df$value), color = "black", size = valSIZE)+

    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(0.5, 0.7),
        legend.direction = "horizontal")+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                     title.position = "top", title.hjust = 0.5))
    
    png(sprintf("%s/%s_%s.png",DIR,subDIR,shortTIT),width=WIDTH,height=HEIGHT,units="px");
    print(gplot1)
    dev.off()

    return(gplot1)
}

#---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- DIRS[grepl("Clustering",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]

#---Check or create output dir
#eldir <- sprintf("%s/%s",dataDIR,subDIR[S])
#if( !file_test("-d",eldir) ){
#    dir.create(eldir)
#}

cldir <- sprintf("%s/%s",OUT[2],subDIR[S])
if( !file_test("-d",cldir) ){
    dir.create(cldir)
}

grdir <- sprintf("%s/%s",OUT[3],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---


#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
ids <- V(gg)$name;

set  <- c("fc","fc2","lec","lec2","SVI","louvain","louvain2","infomap","sgG1","sgG12","Spectral","wt","wt2")
algs <- ALGS[match(set,ALGS)]
#algs <- ALGS[-c(11,12,13)]
#algs <- ALGS[-c(11,12,13,15,16,17,22,23,24,25,26,27,28,29)]
N    <- length(algs)   

run <- vector(length=2)
run[1] = 1
run[2] = 0

if( run[1] ){

    oo1 <- matrix(NA, nrow=(N), ncol=(N))
    colnames(oo1) <- algs
    rownames(oo1) <- algs

    oo2 <- matrix(NA, nrow=(N), ncol=(N))
    colnames(oo2) <- algs
    rownames(oo2) <- algs

    oo3 <- matrix(NA, nrow=(N), ncol=(N))
    colnames(oo3) <- algs
    rownames(oo3) <- algs


    for( i in 1:N ){
        res1 = igraph::get.vertex.attribute(gg,algs[i])    
        for( j in i:N ){
            res2 = igraph::get.vertex.attribute(gg,algs[j])
            rr = compare(res1,res2,"sqrt")
            
            oo1[i,j] = rr@RI
            oo1[j,i] = oo1[i,j]
            
            oo2[i,j] = rr@ARI
            oo2[j,i] = oo2[i,j]
        
            oo3[i,j] = rr@NMI
            oo3[j,i] = oo3[i,j]
        
        }
    }
   
    write.table( oo1, file=sprintf("%s/RI.csv",OUT[2]), append=F, row.names=T, col.names=T, sep="\t", quote=F)
    write.table( oo2, file=sprintf("%s/ARI.csv",OUT[2]), append=F, row.names=T, col.names=T, sep="\t", quote=F)
    write.table( oo3, file=sprintf("%s/NMI.csv",OUT[2]), append=F, row.names=T, col.names=T, sep="\t", quote=F)

#--- plotting
#size=2
    size=3
    pt1 = plotMetric(oo1, OUT[2], subDIR[S], "Rand Index", "RI", valSIZE=size)
    pt2 = plotMetric(oo2, OUT[2], subDIR[S], "Adjusted Rand Index", "ARI",valSIZE=size)
    pt3 = plotMetric(oo3, OUT[2], subDIR[S], "Normalized Mutual Information", "NMI",valSIZE=size)

}

if( run[2] ){

    RES    = list()
    GPLOTS = list()
    MIN    = 0
    MAX    = 0
    for( i in 1:N ){
        RES[[i]] = icf(length(ids), gg, algs[i] )
        names(RES)[i] = algs[i]
        xx = as.numeric(as.vector(unlist(RES[[i]])[,3]))
        if( i == 1 ){
            MIN = min(xx)
            MAX = max(xx)           
        } else {
            if( min(xx) < MIN ){ MIN = min(xx) }
            if( max(xx) > MAX ){ MAX = max(xx) }
        }
        
    }

    for( i in 1:N ){
        
        xx = as.numeric(as.vector(unlist(RES[[i]])[,3]))
        
        GPLOTS[[i]] = icfPlots(xx, OUT[2], subDIR[S], algs[i], "ICF", MIN, MAX, BINS=20)
        names(GPLOTS)[i] <- LETTERS[i]
        
    }

    p <- plot_grid(plotlist=GPLOTS, labels = "AUTO", label_size=50, label_fontface="bold")
    ggsave(sprintf("%s/%s_Plots.png",OUT[2],"ICF"), p, width=20, height=20, device="png")

    
}
