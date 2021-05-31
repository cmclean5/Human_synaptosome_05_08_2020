source('../setUp.R')

#Check or create out dir
if( !file_test("-d",subDIR[S]) ){
    dir.create(subDIR[S])
}
    
#---Clustering algorithms used 
#set  <- c("fc","fc2","lec","lec2","louvain","louvain2","infomap","sgG1","sgG12","Spectral","wt","wt2")
set  <- c("fc","lec","louvain","infomap","sgG1","Spectral","wt")
algs <- ALGS[match(set,ALGS)]
#algs  <- ALGS[c(1,2,5,6,7,8,9,14,18,19,20,21)]
#algs <- ALGS[-c(3,4,10,11,12,13,15,16,17,22,23,24,25,26,27,28,29)]


#set CDF parameters
steps  = 100;
stepsi = 1.0/steps;
CDF <- vector(length=steps);

X   <- vector(length=steps)
for(x in 1:steps ){
    Xi=stepsi*x;
    X[x] <- Xi;
}      


#--- load the consensusCDF file generated from running 'buildConsensusCDF.R'
oo <- read.delim(sprintf("%s/%s_consensusCDFs.csv",subDIR[S],subDIR[S]),header=T,sep="\t")
    
   
#---proportion of ambiguously clustered pairs (PAC)
# Senbabaoglu Y., Michailidis G. and Li J. Critical Limitations of consensus clustering in class discovery. Scientific Reports, 4, 6207, 2014. doi:10.1038/srep06207

#set parameters
X1  = 10 #0.1
X2  = 90 #0.9
PAC = -1

N   = length(oo[1,])

PAC    <- as.numeric(oo[X2,2:N]) - as.numeric(oo[X1,2:N])
ALG    <- colnames(oo)[2:N]
MINCDF <- as.numeric(oo[1,2:N])

STATS <- matrix("",ncol=4,nrow=length(ALG))

#print PAC results for each algorithm
cat("ALG Name \t=\t PAC value \t=\t CDFmin\n")
for( i in 1:length(ALG) ){
    cat(ALG[i], "\t=\t", PAC[i], "\t=\t", MINCDF[i],"\n")
    
    Col = "grey"#royalblue"
    #if( MINCDF[i] >= 0.1 ){
    #    Col = "grey70"
    #}

    STATS[i,1] = ALG[i]
    STATS[i,2] = PAC[i]
    STATS[i,3] = MINCDF[i]
    STATS[i,4] = Col
    
}
cat("\n")


#find algorithm with minimum PAC value
PACmin <- which.min(PAC)

#print algorithm with min PAC value
cat("Algorithm ", ALG[PACmin], " has min PAC value of ", PAC[PACmin],"\n")


#--- order algorithms from lowest PAC values
STATS <- STATS[order(as.numeric(STATS[,2])),]
indx  <- match(c("X",STATS[,1]),colnames(oo))

oo2   <- oo[,indx]
colnames(oo2) <- colnames(oo)[indx]
colnames(oo2)[colnames(oo2) == "Spectral1per"]  = "Spectral_1per"
STATS[STATS[,1]=="Spectral1per",1] = "Spectral_1per"
colnames(oo2)[colnames(oo2) == "Spectral25per"] = "Spectral_2.5per"
STATS[STATS[,1]=="Spectral25per",1] = "Spectral_2.5per"
colnames(oo2)[colnames(oo2) == "Spectral5per"]  = "Spectral_5per"
STATS[STATS[,1]=="Spectral5per",1] = "Spectral_5per"
colnames(oo2)[colnames(oo2) == "sgG12"]         = "sgG1_2"
STATS[STATS[,1]=="sgG12",1] = "sgG1_2"

oo <- oo2

colnames(STATS) <- c("Alg","PAC","minCDF","Colour")
write.table(STATS, file=sprintf("%s/%s_summarySTATS.csv",subDIR[S],subDIR[S]), append=F, row.names=F, col.names=T, sep="\t",quote=F);
#----

#---
# Plot consensus CDF for each algorithm
#---

CN <- colnames(oo)[2:N]
X  <- length(oo[,1])


#---reshape the oo (i.e. the output) dataframe to use in ggplot 
df <- data.frame(a=as.character(),b=as.character(),c=as.numeric(),d=as.numeric(),e=as.numeric(),f=as.numeric())

for( i in 2:N ){

    label <- as.character(rep(CN[i-1],X))
    col   <- as.character(rep(STATS[i-1,4]),X)
    alpha <- as.character(rep(0.7),X)
    xx    <- as.numeric(oo[,1])
    cdf   <- as.numeric(oo[,i])
    size  <- as.character(rep(1.8),X)
    
    df    <- rbind(df,data.frame(label,col,xx,cdf,alpha,size))
    
}

colnames(df) <- c("ALG","COL","X","CDF","ALPHA","LSIZE")

df     <- df[order(match(df[,1],colnames(oo)[2:N])),]
df$ALG <- factor(df$ALG, levels=colnames(oo)[2:N])
df     <- as.data.frame(df)

df$COL[df$ALG == "Spectral"]="royalblue"
df$COL[df$ALG == "louvain2"]="magenta"

df$ALPHA[df$ALG == "Spectral"]= 1.0
df$ALPHA[df$ALG == "louvain2"]= 1.0

df$LSIZE[df$ALG == "Spectral"]= 2.0
df$LSIZE[df$ALG == "louvain2"]= 2.0
#---

#---Check or create plot dir
plotDIR <- sprintf("%s/PLOTS/",subDIR[S])

if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}

#---set ploting colours
#colOrder<- c(levels(factor(df$ALG)))
#alpha   <- rep("", length(colOrder))

#STATS <- STATS[order(as.numeric(STATS[,2])),]

#for( i in 1:length(STATS[,1]) ){

#    STATS[i,4] = "0.6"
#    if( STATS[i,1] == "Spectral_1per" ){ STATS[i,4] = "1" }
#    if( STATS[i,1] == "louvain2" )     { STATS[i,4] = "1" }
    
#}

#alpha <- STATS[match(colOrder,STATS[,1]),4]
#---

#---Generate CDF plot
SIZEa=2
SIZEb=2
library(viridis)
library(gghighlight)

gplot <- ggplot(df,aes(x=(as.numeric(df$X)),y=as.numeric(as.vector(df$CDF)),colour=df$ALG))+
    geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
    labs(x="c",y="CDF(c)")+
    theme(axis.title.x=element_text(face="bold",size=rel(2.0)),
          axis.title.y=element_text(face="bold",size=rel(2.0)),
          legend.text=element_text(face="bold",size=rel(1.5)),
          )+
    scale_color_viridis("",discrete = TRUE, option = "D")+
    theme(panel.grid.major = element_line(colour = "grey40"),
    panel.grid.minor = element_line(colour="grey40",size=0.1),
    panel.background = element_rect(fill = "white"),
     panel.border = element_rect(linetype="solid",fill=NA))+
    geom_vline(xintercept=(as.numeric(X1/steps)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
    geom_vline(xintercept=(as.numeric(X2/steps)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
    guides(color = guide_legend(override.aes = list(size=4)),
           alpha = FALSE,
           size  = FALSE)


png(sprintf("%s/%sCDF.png",plotDIR,subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot)
dev.off()

gplot2 <- ggplot(df,aes(x=log10(as.numeric(df$X)),y=as.numeric(as.vector(df$CDF)),colour=df$ALG))+
    geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
    labs(x="log(c)",y="CDF(c)")+
    theme(axis.title.x=element_text(face="bold",size=rel(2.0)),
          axis.title.y=element_text(face="bold",size=rel(2.0)),
          legend.text=element_text(face="bold",size=rel(1.5)),
          )+
    scale_color_viridis("",discrete = TRUE, option = "D")+
    theme(panel.grid.major = element_line(colour = "grey40"),
    panel.grid.minor = element_line(colour="grey40",size=0.1),
    panel.background = element_rect(fill = "white"),
     panel.border = element_rect(linetype="solid",fill=NA))+
    geom_vline(xintercept=log10(as.numeric(X1/steps)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
    geom_vline(xintercept=log10(as.numeric(X2/steps)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
    guides(color = guide_legend(override.aes = list(size=4)),
           alpha = FALSE,
           size  = FALSE)


png(sprintf("%s/%sCDFlog.png",plotDIR,subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot2)
dev.off()

colours = df$COL[match(levels(factor(df$ALG)),df$ALG)]

gplot3 <- ggplot(df,aes(x=log10(as.numeric(df$X)),y=as.numeric(as.vector(df$CDF)),colour=df$ALG))+
    geom_line(size=as.numeric(as.vector(df$LSIZE)),alpha=as.numeric(as.vector(df$ALPHA)))+
    labs(x="log(c)",y="CDF(c)")+
    theme(axis.title.x=element_text(face="bold",size=rel(2.0)),
          axis.title.y=element_text(face="bold",size=rel(2.0)),
          legend.text=element_text(face="bold",size=rel(1.5)),
          )+
    #scale_color_viridis("",discrete = TRUE, option = "D")+
     scale_color_manual("",breaks=c(levels(factor(df$ALG))),values=c(colours))+
    theme(panel.grid.major = element_line(colour = "grey40"),
    panel.grid.minor = element_line(colour="grey40",size=0.1),
    panel.background = element_rect(fill = "white"),
     panel.border = element_rect(linetype="solid",fill=NA))+
    geom_vline(xintercept=log10(as.numeric(X1/steps)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
    geom_vline(xintercept=log10(as.numeric(X2/steps)),colour="grey10",size=SIZEb,linetype=2,show.legend=F)+
    guides(color = guide_legend(override.aes = list(size=4)),
           alpha = FALSE,
           size  = FALSE)

png(sprintf("%s/%sCDFlogv2.png",plotDIR,subDIR[S]),width=WIDTH,height=HEIGHT,units="px");
print(gplot3)
dev.off()
