
source('../setUp.R')

##---Semi-local Centrality (Cl)
##   Identifying influential nodes in complex networks, D. Chen et al., Physica A, 2012
Semilocal <- function(gg){

N    <- length(V(gg)$name)
meas <- matrix(0, nrow=N, ncol=3)

for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- igraph::neighbors(gg,v=ids,mode="all")

    if( length(neig) > 0 ){
  
        for( w in 1:length(neig) ){
            neig <- c(neig,igraph::neighbors(gg,v=as.character(V(gg)$name[neig[w]]),mode="all"))    
        }

        neig <- unique(neig)

        meas[i,1] <- length(neig)-1

  }
    
}

for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)

    neig <- igraph::neighbors(gg,v=ids,mode="all")
  
    meas[i,2] <- sum(as.numeric(meas[neig,1]))
  
}


for( i in 1:N ){

    ids <- as.character(V(gg)[i]$name)
  
    neig <- igraph::neighbors(gg,v=ids,mode="all")
  
    meas[i,3] <- sum(as.numeric(meas[neig,2]))
  
}

 return(as.numeric(meas[,3]))

}


##calculate the mean and sd of the shortest paths for each gene
calShorestPaths <- function(gg){

 N    <- length(V(gg)$name)
 meas <- matrix(0, nrow=N, ncol=3)
		
 for( i in 1:N ){
      sp <- as.numeric(igraph::shortest.paths(gg,i))
      sp <- sp[-i]
      sp <- sp[!sp == Inf]
      meas[i,1] <- min(sp)
      meas[i,2] <- round(mean(sp),3)
      meas[i,3] <- round(sd(sp),3)
 } 		

 return(meas)

}


filterZeroDegree <- function( DIR, SUB ){

    ft <- list()

    for( d in 1:length(DIR) ){
        tt <- read.table(sprintf("%s/random_%s_permuteDEGREE.csv",DIR[d],SUB),sep="\t",header=T)

        Nc <- length(tt[1,])
        tt <- tt[,2:Nc]
        
        tmp <- apply(tt, 2, function(x) {ifelse(x == 0, NA, 1)})

        ft[[d]] <- tmp    
        
        names(ft)[d] <- sprintf("RANDOM%d",d)    
        
    }

    return(ft)
    
}


readRandomDataFiles <- function( DIR, SUB, CENT, C, FILTER , COLN ){

    ranDF <- list()

    for( d in 1:length(DIR) ){
    
        if( C == 6 ){
            tt <- read.table(sprintf("%s/random_MEAN_%s_permute%s.csv",DIR[d],SUB,CENT[C]),sep="\t",header=T)
        } else {
            tt <- read.table(sprintf("%s/random_%s_permute%s.csv",DIR[d],SUB,CENT[C]),sep="\t",header=T)
        }
        

        Nc <- length(tt[1,])   

        tt <- tt[,2:Nc]

        if( !is.null(FILTER) ){        
            tt <- tt*FILTER[[d]]
        }

        ##test
        tt  <- tt[,1]
        
        tmp <- as.numeric(unlist(tt))
        
        tmp <- tmp[!is.na(tmp)]        

        tmp <- data.frame(x=as.numeric(tmp),group=COLN[d])

        ranDF[[d]] <- tmp    
        
        names(ranDF)[d] <- sprintf("RANDOM%d",d)    

        rm(tt,Nc,tmp)
        
    }

    return(ranDF)
    
}


formatLogLogPlot <- function( X, GROUP ){

    X = as.vector(X)

    mm <- ecdf(X)

    df <- data.frame(x=sort(X),y=1-mm(sort(X)),group=GROUP)
   
    return(df)
    
}

filterLogLog <- function( df, xMAX, yMIN ){

    if( !is.null(xMAX) ){
        df <- df[ df$x <= xMAX, ]
    }

    if( !is.null(yMIN) ){
        df <- df[ df$y >= yMIN, ]
    }

    indx <- which(df$y==0)

    if( length(indx) != 0 ){
        df <- df[-indx,]
    }
    
    return(df)
}

##Calculate the Median absolute difference
MAD <- function( X ){

    X <- as.numeric(X)
    
    Xmd <- median(X)

    MAD <- median(abs(X-Xmd))

    return(MAD)
}

##---Directories
OUT     <- vector(length=3)
OUT[1]  <- DIRS[grepl("Graphs",DIRS)]
OUT[2]  <- sprintf("%s/POWERlawFIT/RANDOM1",rndDIR[1]) #alpha == Inf, i.e. Erdos-Renyi model
OUT[3]  <- sprintf("%s/POWERlawFIT/RANDOM2",rndDIR[1]) #obs alpha from fits,
OUT[4]  <- sprintf("%s/POWERlawFIT/RANDOM3",rndDIR[1]) #alpha == 2 

subTIT    <- vector(length=3)
subTIT[1] <- ""#"Presynaptic"
subTIT[2] <- ""#"PSD"
subTIT[3] <- ""#"PSD"

##---Set Options
runCen   <- vector(length=7)
runCen[1] <- 1   #Generate Centrality RSD file
runCen[2] <- 0   #Plots over single study
runCen[3] <- 0   #PLOT1
runCen[4] <- 0   #PLOT2
runCen[5] <- 0   #PLOT3
runCen[6] <- 0   #PLOT4
runCen[7] <- 0  #Superimpose each network centrality  


if( runCen[1] ){

gsO = setNames(vector("list",length(subDIR)),nm=subDIR)

for( i in 1:length(subDIR) ){

    str <- sprintf("%s/%s/%s.gml",OUT[1],subDIR[i],subDIR[i])

    if( file.exists(str) && file.info(str)$size!=0 ){    
        gg <- igraph::read.graph(str,format="gml")
        gsO[[i]] <- gg
    } else {
        gsO[[i]] <- NULL
    }

}

measO <- list()

for( i in 1:length(subDIR) ){    
    
    gg <- gsO[[i]]

    if( !is.null(gg) ){
    
        ID <- V(gg)$name
        N  <- length(ID)
    
        CN  <- c("ID","DEG","BET","CC","SL","mnSP","PR","sdSP")
    
        tmp <- matrix(0,nrow=N,ncol=length(CN))
        colnames(tmp) <- CN

        tmp[,1] <- ID
        tmp[,2] <- as.vector(igraph::degree(graph=gg))
        tmp[,3] <- as.character(round(betweenness(gg),3))
        tmp[,4] <- as.character(round(transitivity(gg,"local"),3))
        tmp[,5] <- as.character(round(Semilocal(gg),3))
    
        res <- as.matrix(calShorestPaths(gg))
        tmp[,6]  <- as.character(res[,2])
        tmp[,7]  <- as.character(round(as.vector(page.rank(graph=gg,vids=V(gg),directed=F,options=igraph.arpack.default)$vector),6))
        tmp[,8]  <- as.character(res[,3])
    
        measO[[i]]      <- tmp
        names(measO)[i] <- names(gsO)[i]
    } else {
        measO[[i]]      <- NULL
        names(measO)[i] <- ""
    }
}

    saveRDS(measO,"originalNtwrkMeasures.rds")

} else {

    measO <- readRDS("originalNtwrkMeasures.rds")

}
    

cent    <- vector(length=7)
cent[1] <- ""
cent[2] <- "DEGREE"
cent[3] <- "BET"
cent[4] <- "CC"
cent[5] <- "SL"
cent[6] <- "SP"
cent[7] <- "PR"

Xlabs    <- vector(length=7)
Xlabs[1] <- ""
Xlabs[2] <- "log(k)"   #"Degree"
Xlabs[3] <- "log(Bet.)"#"Betweenness"
Xlabs[4] <- "log(CC)"  #"Clustering Co."
Xlabs[5] <- "log(SL)"  #"Semi-local"
Xlabs[6] <- "log(SP)"  #"Mean Shortest Path"
Xlabs[7] <- "log(PR)"  #"PageRank"

Colours <- vector(length=4)
Colours[1] <- "darkmagenta"
Colours[2] <- "deepskyblue4"
Colours[3] <- "deepskyblue4"
Colours[4] <- "deepskyblue4"

FillColour <- vector(length=4)
FillColour[1] <- "black"#"magenta"
FillColour[2] <- "tan2"#"deepskyblue4"
FillColour[3] <- "steelblue3"
FillColour[4] <- "turquoise3"

Group <- vector(length=4)
Group[1] <- "A"
Group[2] <- "B"
Group[3] <- "C"
Group[4] <- "D"

Leg <- vector(length=4)
Leg[1] <- c(expression(paste(PPI[obs])))
Leg[2] <- c(expression(paste("E-R")))
Leg[3] <- c(expression(paste(alpha,"=",PPI[obs])))
Leg[4] <- c(expression(paste(alpha,"=",2)))

Smoothing=4

Alpha <- vector(length=4)
Alpha[1] <- 0.8
Alpha[2] <- 0.3
Alpha[3] <- 0.3
Alpha[4] <- 0.3

Sample=NULL
Zero=NULL
Log=NULL
Cname="x"

Xmax    <- vector(length=7)
Xmax[1] <- 1
Xmax[2] <- 75
Xmax[3] <- 1000
Xmax[4] <- 1
Xmax[5] <- 10000
Xmax[6] <- 7
Xmax[7] <- 0.005


if( runCen[2] ){

    study <- subDIR[S]

    inds <- which(names(measO)==study)
    
    ##---check or create plot Dir
    plotDIR <- sprintf("PLOTS/%s",study)

    if( !file_test("-d",plotDIR) ){
        dir.create(plotDIR)
    } 
    
    #Filter (mask) all zero degree nodes from measures
    Filter <- filterZeroDegree(OUT[2:4], study);
    
    for( c in 2:3){#length(cent) ){

        CC <- as.numeric(measO[[inds]][,c])
        CC <- CC[!is.na(CC)]
        
        CC <- data.frame(x=as.numeric(CC),group=Group[1])

        RR <- readRandomDataFiles( OUT[2:4], study, cent, c, Filter, Group[2:4] )

        CCp  <- formatLogLogPlot(CC$x,levels(as.factor(CC$group)))
        RR1p <- formatLogLogPlot(RR[[1]]$x,levels(as.factor(RR[[1]]$group)))
        RR2p <- formatLogLogPlot(RR[[2]]$x,levels(as.factor(RR[[2]]$group)))
        RR3p <- formatLogLogPlot(RR[[3]]$x,levels(as.factor(RR[[3]]$group)))

        dat  <- rbind(CCp,RR1p,RR2p,RR3p)

        CCp  <- filterLogLog( CCp,  NULL, NULL )
        
        MAX  <- max(CCp$x)
        MIN  <- min(CCp$y)
        
        RR1p <- filterLogLog( RR1p, MAX,  MIN  )
        RR2p <- filterLogLog( RR2p, MAX,  MIN  )
        RR3p <- filterLogLog( RR3p, MAX,  MIN  )

        dat2 <- rbind(CCp,RR1p,RR2p,RR3p)
        ##write.table(dat2, "dat2.csv", sep="\t", row.names=F, col.names=T, quote=F)
        
        Alpha[2:4] <- 0.3
        
        ## Non - Log Plots
        #PLOT1
        if( runCen[3] ){
            gplot <- ggplot( dat, aes(dat$x, fill=dat$group, colour=dat$group, alpha=dat$group)) +
                geom_density(aes(y=..scaled..),adjust=4)+        
                scale_colour_manual(name = "", values=Colours,   labels=Leg) +
                scale_alpha_manual(name = "",  values=Alpha,     labels=Leg) +
                scale_fill_manual(name = "",   values=FillColour,labels=Leg)+
                guides(color = FALSE,
                       alpha = FALSE,
                       fill  = guide_legend(order = 1),
                       group = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+                
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="Density",title="")+
                xlim(c(0,Xmax[c]))            

            png(sprintf("%s/%s_%s_Density.png",plotDIR,study,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot)
            dev.off()
        }#PLOT1

        #PLOT2
        if( runCen[4] ){
            gplot2 <- ggplot( dat, aes(log10(dat$x), fill=dat$group, colour=dat$group, alpha=dat$group)) +
                geom_density(aes(y=..scaled..),adjust=4)+            
                scale_colour_manual(name = "", values=Colours,   labels=Leg) +
                scale_alpha_manual(name = "",  values=Alpha,     labels=Leg) +
                scale_fill_manual(name = "",   values=FillColour,labels=Leg)+
                guides(color = FALSE,
                       alpha = FALSE,
                       fill  = guide_legend(order = 1),
                       group = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="Density",title="")

            png(sprintf("%s/%s_%s_DensityLogx.png",plotDIR,study,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot2)
            dev.off()
        }#PLOT2

        ##Log-Log-Plots       
        Alpha[2:4] <- 0.6

        #PLOT3 - no filtering
        if( runCen[5] ){
            
            dat$x <- log10(dat$x)
            dat$y <- log10(dat$y)
            
            gplot3 <- ggplot( dat, aes(x=dat$x, y=dat$y, colour=dat$group, alpha=dat$group)) +
                geom_line(size=2)+
                scale_colour_manual(name = "", values=FillColour, labels=Leg) +
                scale_alpha_manual(name = "", values=Alpha, labels=Leg) +            
                guides(color = guide_legend(order = 1, override.aes = list(alpha=1,size=1,shape=0)),
                       fill  = FALSE,
                       group = FALSE,
                       alpha = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+        
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+
                labs(x=Xlabs[c],y="",title="")

            png(sprintf("%s/%s_%s_LogLog.png",plotDIR,study,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot3)
            dev.off()
        }#PLOT3
        
        #PLOT4 - filtering
        if( runCen[6] ){
            
            #dat2$x <- log10(dat2$x)
            #dat2$y <- log10(dat2$y)

            Alpha[2:4] <- 0.6
            
            gplot4 <- ggplot( dat2, aes(x=dat2$x, y=dat2$y, colour=dat2$group, alpha=dat2$group)) +
                geom_line(size=2)+
                scale_colour_manual(name = "", values=FillColour, labels=Leg) +
                scale_alpha_manual(name = "", values=Alpha, labels=Leg) +            
                guides(color = guide_legend(order = 1, override.aes = list(alpha=1,size=2.5,shape=0)),
                       fill  = FALSE,
                       group = FALSE,
                       alpha = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0,                    
                    legend.key=element_blank())+        
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="",title="")+
                scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                              labels = trans_format("log10", math_format(10^.x))) +
                scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                              labels = trans_format("log10", math_format(10^.x))) +
                annotation_logticks() 

            png(sprintf("%s/%s_%s_LogLogZoom.png",plotDIR,study,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot4)
            dev.off()
            
        }#PLOT4

    }

}

if( runCen[7] ){

    ##---Check or create plot dir
    plotDIR <- "PLOTS/NETWORKS"

    if( !file_test("-d",plotDIR) ){
        dir.create(plotDIR)
    }
     
    INDX <- match(subDIR,names(measO))
    
    FillColour[1] <- "deeppink"
    FillColour[2] <- "cyan3"
    FillColour[3] <- "limegreen"
    FillColour[4] <- ""

    Colours[1] <- "deeppink"
    Colours[2] <- "cyan3"
    Colours[3] <- "limegreen"
    Colours[4] <- ""
    
    Leg[1] <- c(expression(paste(PPI[pre])))
    Leg[2] <- c(expression(paste(PPI[psp])))
    Leg[3] <- c(expression(paste(PPI[pspRED])))
    Leg[4] <- c("")
    
    TIT=""
    
    for( c in 2:length(cent) ){

        CC1 <- as.numeric(measO[[INDX[1]]][,c])
        CC1 <- CC1[!is.na(CC1)]
        
        CC1 <- data.frame(x=as.numeric(CC1),group=Group[1])

        CC2 <- as.numeric(measO[[INDX[2]]][,c])
        CC2 <- CC2[!is.na(CC2)]
        
        CC2 <- data.frame(x=as.numeric(CC2),group=Group[2])

        CC3 <- as.numeric(measO[[INDX[3]]][,c])
        CC3 <- CC3[!is.na(CC3)]
        
        CC3 <- data.frame(x=as.numeric(CC3),group=Group[3])

        CC1p  <- formatLogLogPlot(CC1$x,levels(CC1$group))
        CC2p  <- formatLogLogPlot(CC2$x,levels(CC2$group))
        CC3p  <- formatLogLogPlot(CC3$x,levels(CC3$group))
        
        dat  <- rbind(CC1p,CC2p,CC3p)

        CC1p  <- filterLogLog( CC1p,  NULL, NULL )
        
        MAX  <- max(CC1p$x)
        MIN  <- min(CC1p$y)
        
        CC2p <- filterLogLog( CC2p, MAX,  MIN  )
        CC3p <- filterLogLog( CC3p, MAX,  MIN  )

        dat2 <- rbind(CC1p,CC2p,CC3p)
        
        ## Non - Log Plots
        #PLOT1
        Alpha[1:4] <- 0.3
        if( runCen[3] ){
            gplot <- ggplot( dat, aes(dat$x, fill=dat$group, colour=dat$group, alpha=dat$group)) +
                geom_density(aes(y=..scaled..),adjust=4)+            
                scale_colour_manual(name = "", values=Colours,   labels=Leg) +
                scale_alpha_manual(name = "",  values=Alpha,     labels=Leg) +
                scale_fill_manual(name = "",   values=FillColour,labels=Leg)+
                guides(color = FALSE,
                       alpha = FALSE,
                       fill  = guide_legend(order = 1),
                       group = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="Density",title=TIT)+
                xlim(c(0,Xmax[c]))


            png(sprintf("%s/PPI_%s_Density.png",plotDIR,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot)
            dev.off()
        }#PLOT1
            
        #PLOT2
        if( runCen[4] ){
            gplot2 <- ggplot( dat, aes(log10(dat$x), fill=dat$group, colour=dat$group, alpha=dat$group)) +
                geom_density(aes(y=..scaled..),adjust=4)+            
                scale_colour_manual(name = "", values=Colours,   labels=Leg) +
                scale_alpha_manual(name = "",  values=Alpha,     labels=Leg) +
                scale_fill_manual(name = "",   values=FillColour,labels=Leg)+
                guides(color = FALSE,
                       alpha = FALSE,
                       fill  = guide_legend(order = 1),
                       group = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="Density",title=TIT)

            png(sprintf("%s/PPI_%s_DensityLogx.png",plotDIR,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot2)
            dev.off()
        }#PLOT2

        ##Log-Log-Plots       

        Alpha[1:4] <- 0.8

        #PLOT3 - no filtering
        if( runCen[5] ){
            
            dat$x <- log10(dat$x)
            dat$y <- log10(dat$y)

            gplot3 <- ggplot( dat, aes(x=dat$x, y=dat$y, colour=dat$group, alpha=dat$group)) +
                geom_line(size=2)+
                scale_colour_manual(name = "", values=FillColour, labels=Leg) +
                scale_alpha_manual(name = "", values=Alpha, labels=Leg) +            
                guides(color = guide_legend(order = 1, override.aes = list(alpha=1,size=1,shape=0)),
                       fill  = FALSE,
                       group = FALSE,
                       alpha = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+        
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="",title=TIT)

            png(sprintf("%s/PPI_%s_LogLog.png",plotDIR,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot3)
            dev.off()
        }#PLOT3

        
        #PLOT4 - filtering
        if( runCen[6] ){

            dat2$x <- log10(dat2$x)
            dat2$y <- log10(dat2$y)

            gplot4 <- ggplot( dat2, aes(x=dat2$x, y=dat2$y, colour=dat2$group, alpha=dat2$group)) +
                geom_line(size=2)+
                scale_colour_manual(name = "", values=FillColour, labels=Leg) +
                scale_alpha_manual(name = "", values=Alpha, labels=Leg) +            
                guides(color = guide_legend(order = 1, override.aes = list(alpha=1,size=1,shape=0)),
                       fill  = FALSE,
                       group = FALSE,
                       alpha = FALSE)+            
                theme(
                    plot.title = element_text(face="bold",size=rel(2.0)),
                    axis.title.x=element_text(face="bold",size=rel(2.5)),
                    axis.title.y=element_text(face="bold",size=rel(2.5)),
                    legend.title=element_text(face="bold",size=rel(2.0)),
                    legend.text=element_text(face="bold",size=rel(2.5)),
                    legend.text.align=0)+                
                theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                      panel.grid.minor = element_line(colour="grey40",size=0.1),
                      panel.background = element_rect(fill = "white"),
                      panel.border = element_rect(linetype="solid",fill=NA))+    
                labs(x=Xlabs[c],y="",title=TIT)

            png(sprintf("%s/PPI_%s_LogLogZoom.png",plotDIR,cent[c]),width=WIDTH,height=HEIGHT,units="px")
            print(gplot4)
            dev.off()
        }#PLOT4    

    }

}
