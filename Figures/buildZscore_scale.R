source('../setUp.R')
require(cowplot)

scale <- function(X){

    X = as.numeric(as.vector(X))

    xmin <- min(X,na.rm=T)
    xmax <- max(X,na.rm=T)

    return((X-min(X))/(max(X)-min(X)))
    
}

getCDF <- function(X){

    X = as.vector(X)
    
    mn = floor(min(X))
    mx = ceiling(max(X))
    #mx = floor(max(X))

    cdf = vector(length=length(X))

    z = seq(mn,mx,by=0.01)
    P = ecdf(X)
    p=P(z)

    for(i in 1:length(X)){
        #cdf[i] = last(p[z <= X[i]])
        cdf[i] = last(p[z <= DF$zscore[i]])
    }

    return(cdf)

}


#---Directories needed
OUT <- vector(length=2)
OUT[1] <- DIRS[grepl("Graphs",DIRS)]
OUT[2] <- DIRS[grepl("Consensus",DIRS)]

#---Check or create plot dir
pltdir <- "PLOTS"
if( !file_test("-d",pltdir) ){
    dir.create(pltdir)
}


#---Check or create output dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}

#---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[1],subDIR[S],subDIR[S]),format="gml")

ids <- V(gg)$name
n   <- length(ids)

alg = vector(length=1)
alg[1] = "Spectral"
#alg[2] = "louvain2"

DF = read.delim(sprintf("%s/%s/%s_meanClusterProb.csv", OUT[2],subDIR[S],alg[1]),sep="\t", header=T)

oo1 = matrix(0, ncol=3, nrow=nrow(DF))
oo1[,1] = 1
oo1[,2] = DF$C
oo1[,3] = DF$meanCProb

df = as.data.frame(oo1)
colnames(df) <- c("V1","V2","value")

mx  = 1#max(df$value)
mn  = 0#min(df$value)
mdp = 0.5#mean(df$value)

indx  = order(df$value,decreasing=F)
label = round(df$V2[indx],2)

## color low  = royalblue2 = #436EEE
## color high = firebrick2 = #EE2C2C

colLOW ="#FFFFCC"
colMID ="#00EEEE"
colHIGH="#436EEE"

#--- plotting
gplot1 <- ggplot(df, aes(df$V2, df$V1, fill = df$value))+
    labs(x="Community (C)",y="",title="")+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = colLOW, high = colHIGH, mid = colMID, 
                         midpoint = mdp, limit = c(mn,mx),
                         breaks = round(seq(0,0.98,length.out=5),2),
                         labels = round(seq(0,1,length.out=5),2),
                         space = "Lab", 
                         #name="Regional ecdf(z-score)")+
                         #name="Regional diversity")+
                         name="mean Cluster Probability")+
    scale_x_continuous(breaks=seq(min(df$V2),max(df$V2),5))+
    coord_fixed(ratio=2)+
    theme(
        #axis.title.x = element_blank(),
        axis.title.x=element_text(face="bold",size=rel(2.2)),        
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.title=element_text(face="bold",size=rel(2.2)),
        legend.text=element_text(face="bold",size=rel(2.2)),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.ticks = element_blank(),
        legend.box.margin=margin(rep(20,4)),
        legend.justification = c(1, 0),
        legend.position = c(0.74, 1.0),
        legend.direction = "horizontal")+
guides(fill = guide_colorbar(barwidth = 20, barheight = 1,
                                     title.position = "top", title.hjust = 0.5))

ggsave(gplot1,file=sprintf("%s/%s_scale.svg",pltdir,alg[1]), width=2*WIDTHcm, height=HEIGHTcm, units="cm");

png(sprintf("%s/%s_scale.png",pltdir,alg[1]),width=2*WIDTH,height=HEIGHT,units="px");
print(gplot1)
dev.off()


#png("ecdf_for_zscore.png",width=WIDTH,height=HEIGHT,unit="px")
#plot(xy.coords(x=oo1[,3],y=oo1[,4]),xlab="z-score",ylab="cumulative probability")
#dev.off()
