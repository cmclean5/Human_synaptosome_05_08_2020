#run 'mergeRESULTS.R' first

source('../setUp.R')
require(cowplot)

makeBold <- function(src, bolder) {
    require(purrr)
    if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
    src_levels <- levels(src)                                 # retrieve the levels in their order
    temp <- bolder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
    if (all(temp)) {                                         # if so
        b_pos <- purrr::map_int(bolder, ~which(.==src_levels)) # then find out where they are
        b_vec <- rep("plain", length(src_levels))               # make'm all plain first
        b_vec[b_pos] <- "bold"                                  # make our targets bold
        b_vec                                                   # return the new vector
    } else {
        stop("All elements of 'bolder' must be in src")
    }
}


#---Check or create plots dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

plotdir <- sprintf("PLOTS/%s",subDIR[S])

if( !file_test("-d",plotdir) ){
    dir.create(plotdir)
}


run1=0

if( run1 ){

FILES    <- vector(length=2)
FILES[1] <- "Disease_overlap_sig"
FILES[2] <- "random_zscores"

#--- load results
ff <- read.table(sprintf("RESULTS/%s/ovg/%s.csv",subDIR[S],FILES[1]),sep="\t",header=T,stringsAsFactors=F,quote="", check.names=F)

#--- load randomised zscores
tests <- read.table(sprintf("RESULTS/%s/ovg/%s.csv",subDIR[S],FILES[2]),sep="\t",header=F,stringsAsFactors=F,quote="", check.names=F)

#--set disease pair of interest
disA = vector(length=5)
disA[1] = "AD"
disA[2] = "AD"
disA[3] = "AD"
disA[4] = "PD"
disA[5] = "MS"

disB = vector(length=5)
disB[1] = "HTN"
disB[2] = "PD"
disB[3] = "MS"
disB[4] = "HTN"
disB[5] = "HTN"

GPLOTS <- list()
g=1

for( i in 1:length(disA) ){

    dA = disA[i]
    dB = disB[i]
     
    #--- locate disease pair in results file
    indx = (ff[,3] == dA & ff[,7] == dB) | (ff[,3] == dB & ff[,7] == dA)

    #--- Plot 1: Z-score of disA versus disB for 1000 randomised studies on the network model.
    if( sum(indx) == 1 ){

        obsZs <- as.numeric(ff[indx,12])    
        DF    <- as.numeric(as.vector(tests[indx,]))
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
        
        tit = sprintf("%sv%s",dA,dB)

        annoT = grobTree(textGrob(sprintf("z-score = %.3f",round(as.numeric(obsZs),4)),
                                  x=0.05, y=0.35, hjust=0,
                                  gp=gpar(col="black", fontface="bold", fontsize=20)))
    
        gplot1 <- ggplot(DF,aes(x=DF$V1))+geom_histogram(bins=BINS,colour="brown1",fill="cadetblue") +
            annotation_custom(annoT)+
            labs(x=c(expression(paste(" ",S[AB]))),y="Count",title=sprintf("%s Vs %s",as.character(gsub("_"," ",ff[indx,2])),as.character(gsub("_"," ",ff[indx,6]))))+
            ##labs(x=c(expression(paste("Separation ",S[AB]))),y="Count",title=sprintf("%s Vs %s",as.character(gsub("_"," ",ff[indx,2])),as.character(gsub("_"," ",ff[indx,6]))))+
            scale_x_continuous(limit = c(MIN, MAX))+
            theme(
                title=element_text(face="bold",size=rel(1.5)),
                axis.title.x=element_text(face="bold",size=rel(1.5)),
                axis.title.y=element_text(face="bold",size=rel(1.5)),
                legend.title=element_text(face="bold",size=rel(2.0)),
                legend.text=element_text(face="bold",size=rel(2.0)),
                legend.key=element_blank())+
            theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
                  panel.grid.minor = element_line(colour="grey40",size=0.1),
                  panel.background = element_rect(fill = "white"),
                  panel.border = element_rect(linetype="solid",fill=NA))+
            geom_segment(data=data.frame(obsZs), aes(x=obsZs,y=ARROW_Y_HEIGHT,xend=obsZs,yend=0.1),size=2,colour="orange2",arrow=arrow(type = "closed"))
        
        png(sprintf("%s/%s_%s.png",plotdir, tit, subDIR[S]),width=WIDTH,height=HEIGHT,units="px")
        print(gplot1)
        dev.off()

        GPLOTS[[g]]      <- gplot1
        names(GPLOTS)[g] <- LETTERS[g]
        g=g+1
        
    }
}

p <- plot_grid(plotlist=GPLOTS, labels = "AUTO", label_size=50, label_fontface="bold")
ggsave(sprintf("%s/Plots.png",plotdir), p, width=20, height=20, device="png")

}

run2=0

if( run2 ){


sabIDX <- which(grepl("sAB",colnames(ff),fixed=T)==TRUE)

zsIDX  <- which(grepl("zScore",colnames(ff),fixed=T)==TRUE)

BonIDX <- which(grepl("Bonferroni",colnames(ff),fixed=T)==TRUE)

pvIDX  <- which(grepl("pvalue",colnames(ff),fixed=T)==TRUE)

adIDX  <- which(grepl("p.adjusted",colnames(ff),fixed=T)==TRUE)
    
qvIDX  <- which(grepl("q-value",colnames(ff),fixed=T)==TRUE)
    
DDpairs <- sprintf("%s-%s",ff[,3],ff[,7])

#requires all three networks 
#for( s in 1:length(subDIR) ){
#s=2

oo     <- matrix(0, ncol=7, nrow=length(DDpairs) )
colnames(oo) <- c("pairs","sAB","zscore","pvalue","p.adjusted","qvalue","label")
oo[,1] <- as.character(DDpairs)
oo[,2] <- as.character(ff[,sabIDX[1]])
oo[,3] <- as.character(ff[,zsIDX[1]])
oo[,4] <- as.character(ff[,pvIDX[1]])
oo[,5] <- as.character(ff[,adIDX[1]])
oo[,6] <- as.character(ff[,qvIDX[1]])
oo[,7] <- rep("",length(DDpairs))

oo <- oo[oo[,2] != 0 & oo[,3] != 0,]
oo <- oo[!grepl("AUT",oo[,1]),]
    
df <- as.data.frame(oo)
df <- df[df$zscore < 0,]

ReOrder = order(as.numeric(df$pvalue),decreasing=T)
df = df[ReOrder,]   
df$pairs      <- as.vector(factor(df$pairs))
df$sAB        <- as.vector(factor(df$sAB))
df$zscore     <- as.vector(factor(df$zscore))
df$pvalue     <- as.vector(factor(df$pvalue))
df$p.adjusted <- as.vector(factor(df$p.adjusted))
df$qvalue     <- as.vector(factor(df$qvalue))
df$pairs      <- factor(df$pairs, levels=df$pairs)
    
write.table(df,sprintf("%s_test.csv", subDIR[S]), sep="\t", col.names=T, row.names=F, quote=F)

    #plot q.values
    XMIN= 0
    XMAX=-log10(min(as.numeric(df$qvalue)))
    if( XMAX < 50 ){ XMAX = 50 } 

    #Highlight Epi-XX pairs
    #target <- c("AD-HTN","HTN-PD","AD-PD")

    #plot qvalues
 gplot <- ggplot(df,aes(y=as.factor(df$pairs),x=-log10(as.numeric(as.vector(df$qvalue))) ))+
        geom_point(alpha=I(0.8),colour="deepskyblue3", shape=17, size=rel(3), show.legend=F)+
        scale_size_identity()+labs(y="Disease pairs",x="-log10(q-value)")+
     theme(legend.key=element_blank())+
     geom_vline(xintercept=-log10(0.05),colour="grey25",linetype="longdash", alpha=0.5, size=rel(1),show.legend=F)+
     theme_bw()+
     scale_x_continuous(limit = c(XMIN, XMAX))+
     theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
           panel.grid.minor = element_line(colour="grey40",size=0.1),
           panel.background = element_rect(fill = "white"))+
     theme(axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.5)),
           axis.text.y = element_text(size=rel(1.5)),
           axis.title.x=element_text(face="bold",size=rel(2.0)),
           axis.title.y=element_text(face="bold",size=rel(2.0)),
           legend.title=element_text(face="bold",size=rel(1.5))
           )+
     geom_text_repel(aes(label=as.vector(df$label)),color="grey40",force=2.0, point.padding = unit(0.8,"lines"),fontface="bold",size=rel(5.0),show.legend=F)
    
png(sprintf("%s/disease_pairs_%s.png",plotdir, subDIR[S]),width=WIDTH,height=2*HEIGHT,units="px");
print(gplot)
dev.off()

}

run3 = 1

if( run3 ){

    df1 = read.table("PPI_Presynaptic_test.csv",sep="\t", header=T)
    df1$label = "pre"
    df2 =  read.table("PPI_PSP_test.csv",sep="\t", header=T)
    df2$label = "PSP"
    df3 =  read.table("PPI_PSP_consensus_test.csv",sep="\t", header=T)
    df3$label = "PSPc"

    #groups <- c("PSP")
    #groups <- c("pre","PSP")
    groups <- c("pre","PSP","PSPc")

    if( length(groups) == 1 ){
        DF = rbind(df2[,c(1,4,5,6,7)])
        df = DF
        colours <- c('royalblue2')
        shapes  <- c(17)
    }
    
    if( length(groups) == 2 ){        
        DF = rbind(df1[,c(1,4,5,6,7)],df2[,c(1,4,5,6,7)])
        df = DF
        colours <- c('firebrick2','royalblue2')
        shapes  <- c(16,17)
    }

    if( length(groups) == 3 ){        
        DF = rbind(df1[,c(1,4,5,6,7)],df2[,c(1,4,5,6,7)],df3[,c(1,4,5,6,7)])
        df = DF
        colours <- c('firebrick2','royalblue2','green2')
        shapes  <- c(16,17,15)
    }
     

    #Highlight Epi-XX pairs
    #target <- c("SCH-Epi","ASD-Epi","Epi-MS","BD-Epi")
    target <- c("AD-HTN","HTN-PD","AD-PD")
    #target <- unique(as.vector(df$pairs[grepl("HTN",as.vector(df$pairs))]))

    ##Remove HTN
    ##df <- df[!grepl("HTN",df[,1]),]
    
    ReOrder = order(as.numeric(df$qvalue),decreasing=T)
    df = df[ReOrder,]
    df$pairs      <- as.vector(factor(df$pairs))
    df$pvalue     <- as.vector(factor(df$pvalue))
    df$p.adjusted <- as.vector(factor(df$p.adjusted))
    df$qvalue     <- as.vector(factor(df$qvalue))
    df$label      <- as.vector(factor(df$label))
    df$pairs = factor(df$pairs, levels=unique(df$pairs))

    XMIN= 0
    XMAX=-log10(min(as.numeric(df$qvalue)))
    if( XMAX < 50 ){ XMAX = 50 }    

    #plot qvalues
    gplot <- ggplot(df,aes(y=as.factor(df$pairs),x=-log10(as.numeric(as.vector(df$qvalue))),group=df$label) )+
        geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(3))+
        #geom_errorbarh(aes(xmin=(-log10(as.numeric(as.vector(df$qvalue)))), xmax=-log10(as.numeric(as.vector(df$qvalue)))), height=.2, width=.1, alpha=I(0.8))+
    scale_size_identity()+labs(y="Disease pairs",x="-log10(q-value)")+
        theme(legend.key=element_blank())+
        geom_vline(xintercept=-log10(0.05),colour="grey25",linetype="longdash", alpha=0.5, size=rel(1),show.legend=F)+
        theme_bw()+
        scale_x_continuous(limit = c(XMIN, XMAX))+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"))+
        theme(axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.5)),
              axis.text.y = element_text(face=makeBold(df$pairs,target),size=rel(1.5)),
              axis.title.x=element_text(face="bold",size=rel(2.0)),
              axis.title.y=element_text(face="bold",size=rel(2.0)),
              legend.title=element_text(face="bold",size=rel(1.7)),
              legend.text=element_text(face="bold",size=rel(1.7)),
              legend.position="top"
              )+
        guides(colour = guide_legend(override.aes = list(shape = shapes, size=rel(7))),
               size   = FALSE,
               shape  = FALSE)+
        scale_color_manual("Compartment",breaks=levels(factor(df$label)),values=c(colours))

    #ggsave(gplot,file=sprintf("%s/%s_scale.svg",pltdir,alg[1]), width=WIDTHcm, height=2*HEIGHTcm, units="cm");
    
    png("PLOTS/disease_pairs_overlaps.png",width=WIDTH,height=2*HEIGHT,units="px");
    print(gplot)
    dev.off()


    #plot p.adjusted
    ReOrder = order(as.numeric(df$p.adjusted),decreasing=T)
    df = df[ReOrder,]
    df$pairs      <- as.vector(factor(df$pairs))
    df$pvalue     <- as.vector(factor(df$pvalue))
    df$p.adjusted <- as.vector(factor(df$p.adjusted))
    df$qvalue     <- as.vector(factor(df$qvalue))
    df$label      <- as.vector(factor(df$label))
    df$pairs = factor(df$pairs, levels=unique(df$pairs))

    XMIN= 0
    XMAX=-log10(min(as.numeric(df$p.adjusted)))
    if( XMAX < 50 ){ XMAX = 50 }    

    #plot qvalues
    gplot <- ggplot(df,aes(y=as.factor(df$pairs),x=-log10(as.numeric(as.vector(df$p.adjusted))),group=df$label) )+
        geom_point(aes(colour=df$label, shape=df$label),alpha=I(0.8), size=rel(3))+
        scale_size_identity()+
        labs(y="Disease pairs",x="-log10(p.adjusted)")+
        theme(legend.key=element_blank())+
        geom_vline(xintercept=-log10(0.05),colour="grey25",linetype="longdash", alpha=0.5, size=rel(1),show.legend=F)+
        theme_bw()+
        scale_x_continuous(limit = c(XMIN, XMAX))+
        theme(panel.grid.major = element_line(colour = "grey40",size=0.2),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"))+
        theme(axis.text.x = element_text(angle = 40, hjust = 1, face="bold", size=rel(1.5)),
              axis.text.y = element_text(face=makeBold(df$pairs,target),size=rel(1.5)),
              axis.title.x=element_text(face="bold",size=rel(2.0)),
              axis.title.y=element_text(face="bold",size=rel(2.0)),
              legend.title=element_text(face="bold",size=rel(1.7)),
              legend.text=element_text(face="bold",size=rel(1.7)),
              legend.position="top"
              )+
        guides(colour = guide_legend(override.aes = list(shape = shapes, size=rel(7))),
               size   = FALSE,
               shape  = FALSE)+
        scale_color_manual("Compartment",breaks=levels(factor(df$label)),values=c(colours))

    #ggsave(gplot,file=sprintf("%s/%s_scale.svg",pltdir,alg[1]), width=WIDTHcm, height=2*HEIGHTcm, units="cm");
    
    png("PLOTS/disease_pairs_overlaps_padj.png",width=WIDTH,height=2*HEIGHT,units="px");
    print(gplot)
    dev.off()
    
}
