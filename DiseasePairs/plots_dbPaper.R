#run 'mergeRESULTS.R' first

source('../setUp.R')
require(cowplot)

plots <- function( INDX, FF, TESTS, DA, DB, STUDY ){

    obsZs <- as.numeric(FF[INDX,12])    
    DF    <- as.numeric(as.vector(TESTS[INDX,]))
    DF    <- as.data.frame(DF)
    names(DF) <- "V1"

    MIN = sign(min(DF)) * (abs(min(DF))+1)
    MAX = sign(max(DF)) * (abs(max(DF))+1)

    if( abs(MIN) < 5 ){
        MIN = -5
    }

    if( abs(MAX) < 5 ){
        MAX = 5
    }
        
    if( obsZs <= MIN ){
        MIN = sign(obsZs) * (abs(obsZs)+1)
    }

    if( obsZs >= MAX ){
        MAX = sign(obsZs) * (abs(obsZs)+1)
    }
       
        
    BINS=20
    ARROW_Y_HEIGHT=550
        
    tit = sprintf("(%s) %sv%s",STUDY,DA,DB)

    annoT = grobTree(textGrob(sprintf("z-score = %.3f",round(as.numeric(obsZs),4)),
                              x=0.05, y=0.35, hjust=0,
                              gp=gpar(col="black", fontface="bold", fontsize=20)))
    
    gplot1 <- ggplot(DF,aes(x=DF$V1))+geom_histogram(bins=BINS,colour="brown1",fill="cadetblue") +
        annotation_custom(annoT)+
        labs(x=c(expression(paste("Separation ",S[AB]))),y="Count",title=sprintf("%s Vs %s",as.character(gsub("_"," ",FF[INDX,2])),as.character(gsub("_"," ",FF[INDX,6]))))+
        scale_x_continuous(limit = c(MIN, MAX))+
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
              panel.border = element_rect(linetype="solid",fill=NA))+
        geom_segment(data=data.frame(obsZs), aes(x=obsZs,y=ARROW_Y_HEIGHT,xend=obsZs,yend=0.1),size=2,colour="orange2",arrow=arrow(type = "closed"))
        
            #png(sprintf("%s/%s_test.png",plotdir, tit),width=WIDTH,height=HEIGHT,units="px")
            #print(gplot1)
            #dev.off()

            #plot(gplot1)

    return(gplot1)

}

#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

plotdir <- "PLOTS/"

if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}


FILES    <- vector(length=2)
FILES[1] <- "Disease_overlap_sig"
FILES[2] <- "random_zscores"

ff    <- list()
tests <- list()

#--- load presynaptic results
ff[[1]] = read.table(sprintf("RESULTS/PPI_Presynaptic/ovg/%s.csv",FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)
names(ff)[1] = "PPI_Pre"

#--- load presynaptic randomised zscores
tests[[1]] = read.table(sprintf("RESULTS/PPI_Presynaptic/ovg/%s.csv",FILES[2]),sep="\t",header=F,stringsAsFactors=F,quote="", check.names=F)
names(tests)[1] = "PPI_Pre"

#--- load psp results
ff[[2]] = read.table(sprintf("RESULTS/PPI_PSP/ovg/%s.csv",FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)
names(ff)[2] = "PPI_PSP"

#--- load psp randomised zscores
tests[[2]] = read.table(sprintf("RESULTS/PPI_PSP/ovg/%s.csv",FILES[2]),sep="\t",header=F,stringsAsFactors=F,quote="", check.names=F)
names(tests)[2] = "PPI_PSP"

#--set disease pair of interest
disA = vector(length=3)
disA[1] = "Epi"
disA[2] = "Epi"
disA[3] = "Epi"
#disA[4] = "Epi"
#disA[5] = "MS"

disB = vector(length=3)
disB[1] = "MS"
disB[2] = "ASD"
disB[3] = "BD"
#disB[4] = "SCH"
#disB[5] = "HTN"

GPLOTS <- list()
g=1

for( i in 1:length(disA) ){

    dA = disA[i]
    dB = disB[i]

    for( j in 1:length(ff) ){
    
        #--- locate disease pair in results file
        indx = (ff[[j]][,3] == dA & ff[[j]][,7] == dB) | (ff[[j]][,3] == dB & ff[[j]][,7] == dA)

        #--- Plot 1: Z-score of disA versus disB for 1000 randomised studies on the network model.
        if( sum(indx) == 1 ){
            
            GPLOTS[[g]] = plots(indx, ff[[j]], tests[[j]], dA, dB, names(ff)[j])
            names(GPLOTS)[g] <- LETTERS[g]
            g=g+1

        }
    }
}
    

p <- plot_grid(plotlist=GPLOTS, labels = "AUTO", label_size=50, label_fontface="bold", nrow=3, ncol=2)
ggsave(sprintf("%s/dbPlots.png",plotdir), p, width=20, height=30, device="png")

