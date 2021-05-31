#source('../setUp.R')
rm(list=ls())

##Set the absolute path to this working directory
version="05_08_2020"
mainDir <- sprintf("/afs/inf.ed.ac.uk/user/c/cmclean5/WORK/DATA/Human_synaptosome_%s",version)

##Get Path to all top-level directories in folder
DIRS <- list.dirs(mainDir,recursive=F)

##---dataDIR, where to read PPI files
dataDIR <- DIRS[grepl("datasets",DIRS)]

##---subDIR for loading and storing PPI Graphs
subDIR <- list.files(path=dataDIR)
subDIR <- subDIR[grep(".csv",subDIR)]
subDIR <- unlist(strsplit(subDIR,".csv"))

##--- Set 
##--- Pre-load Graph of interest, stored in file 'graphs.csv'
pramFILES <- DIRS[grepl("parameterFiles",DIRS)]
Graph <- read.table(sprintf("%s/graphs.csv",pramFILES),header=F,sep="\t",quote="")
S     <- as.vector(Graph[which(as.vector(Graph[,1]) == 1)[1],2])
S     <- match(S,subDIR)


library(igraph);
library(lattice);
library(methods);
library(poweRlaw);
library(scales);
library(grid);
library(latex2exp);
library(stringr);

##---WIDTH and HEIGHT for plots
WIDTH=480
HEIGHT=480

##set igraph as S4
setClass("poweRlaw")

## Gene Set Analysis (GSA) object
setClass(Class="law",representation(
                           fit="displ",
                           p="numeric",
			   alpha="numeric",
                           SDxmin="numeric",
                           SDalpha="numeric")
         )			   

##log sequence of numbers
lseqBy <- function(from=1, to=100000, by=1, length.out=log10(to/from)+1) {
    tmp <- exp(seq(log(from), log(to), length.out = length.out))
    tmp[seq(1, length(tmp), by)]  
}

changeSciNot <- function(n) {

  n <- format(n, scientific = TRUE)

  oo <- strsplit(as.character(n),"e")

  out <- vector(length=length(oo))

  out[1] <- TeX(sprintf("10^{%s}",sub("-0?","-",oo[[1]][2])))
  
  for( i in 2:length(oo) ){

      if(grepl("-",oo[[i]][2])){
          out[i] <- TeX(sprintf("10^{%s}",sub("-0?","-",oo[[i]][2])))
      } else {
          out[i] <- TeX(sprintf("10^{%s}",sub("\\+0?","",oo[[i]][2])))
      }
      
  }

  return(out)
  

}

FitDegree <- function(DEG, title, Nsim, DATAleg, WIDTH, HEIGHT ){

    #WIDTH=480
    #HEIGHT=480

    #tmp = max(DEG)
    #Max = 4*tmp
    
     DEG <- DEG[DEG > 0]
     
     data <- DEG     
     
     m_pl = displ$new(data)

     est = estimate_xmin(m_pl)

     m_pl$setXmin(est)

     gof <- bootstrap_p(m_pl, no_of_sims = Nsim, threads=4)
    
     x_lab="k"    ##degree
     y_lab="P(k)" ## the CDFs P(k) for the PPI network data
     leg_x = max(data)
     leg_y = 1
     
     png(filename=sprintf("PLOTS/%s_cdf.png",title), width=WIDTH, height=HEIGHT, units="px")

     d = plot(m_pl,draw=F)

     Xmax <- max(d$x) - max(d$x)*0.5

     ##build y-axis labels
     yTICKS  = round(lseqBy(min(d$y),1,0.5),4)
     yLABELS = changeSciNot(yTICKS)
         
     plot(m_pl, xlab=sprintf("%s",x_lab), ylab=y_lab, 
     panel.first=grid(col="grey60"), 
     pch=22, bg='black', axes = F, cex.lab = 1.5, yaxt='n' )
     box(col='black')

     axis(1, cex.axis = 1.5, font = 1.5, family = 'arial')
     axis(2, cex.axis = 1.5, font = 1.5, family = 'arial', at=yTICKS, labels=yLABELS)

     lines(m_pl, col=2, lwd=3)

    S1 = round(m_pl$xmin,2)
    S2 = round(m_pl$pars,2)
    S3 = round(gof$p,2)

    sdS1 = round(sd(gof$bootstraps$xmin),0)
    sdS2 = round(sd(gof$bootstraps$pars),2)

    errS1 = str_sub(as.character(sdS1),-1,-1)
    errS2 = str_sub(as.character(sdS2),-1,-1)
    
    
    fitl <- TeX(sprintf("Power-law $\\alpha = %.2f(%s), $k_{min} = %.0f(%s)",S2,errS2,S1,errS1))
    
    legend("bottomleft",c(DATAleg,fitl),lty=c(1,1),lwd=c(4,4),col=c('black',2),merge=TRUE, cex = 1.5)
    
    #legend("topright",c(DATAleg,expression(paste(P(k), " = ",frac(k^alpha,sigma1(alpha,k[min]))))),lty=c(1,1),lwd=c(4,4),col=c('black',2),merge=TRUE, cex = 1.5)
    #text(x=Xmax,y=0.2,substitute(paste(k[min], " = ", s1, ", ", alpha, " = ", s2), list(s1=S1,s2=S2)),cex=0.9)
    
     dev.off()


     return(new("law",
                fit=m_pl,
                p=as.numeric(gof$p),
                alpha=as.numeric(est$pars),
                SDxmin=as.numeric(sd(gof$bootstraps$xmin)),
                SDalpha=as.numeric(sd(gof$bootstraps$pars))
                )
            )
    
     
 }


##---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- DIRS[grepl("Clustering",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]


##No: of bootstrap iterations
bootStrap    <- vector(length=3)
bootStrap[1] <- 100
bootStrap[2] <- 1000
bootStrap[3] <- 5000

b=2 #we'll stick with 1000 iteteractions   

##Legend Titles
Legend <- vector(length=2)
Legend[1] <- "Presynaptic PPI"
Legend[2] <- "PSP PPI"

l=2;
if( grepl("Pre",subDIR[S]) ){
    l=1;
}

##---Check or create plot dir
if( !file_test("-d","PLOTS") ){
    dir.create("PLOTS")
}

##---load corresponding graph which was used to build the consensus matrices from 
gg <- igraph::read.graph(sprintf("%s/%s/%s.gml",OUT[3],subDIR[S],subDIR[S]),format="gml")
    
pFit <- FitDegree( as.vector(igraph::degree(graph=gg)), subDIR[S], bootStrap[b], Legend[l], WIDTH, HEIGHT )

cat(sprintf("%s, bootStrap: %.1f pvalue: %.3f alpha: %.3f \n", subDIR[S], bootStrap[b], pFit@p, pFit@alpha))

