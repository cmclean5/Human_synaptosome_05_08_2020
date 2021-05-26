source('../setUp.R')

require("minpack.lm")
#require("qpcR")
require(cowplot)

### Sigmoid function ### create a function to generate sigmoid pattern
sigmoid <- function(pars, xx){#, lower_asymptote, carrying_capacity, growth_rate, time_max) {

    a = as.numeric(pars[1])
    b = as.numeric(pars[2])
    c = as.numeric(pars[3])
    d = as.numeric(pars[4])

    return( a + ((b-a)/(1+exp(-c*(xx-d)))) )
}

residFun <- function(pars, observed, xx) observed - sigmoid(pars,xx)


plotSigmoid <- function( x, rates, model, alg="", pv=0 ){

    conf     <- NULL
    try(conf <- confint(model), FALSE)

    if( !is.null(conf) ){
        # adding the 95% confidence interval around the fitted coefficient
        lower = list(a=conf[1, 1], b=conf[2, 1], c=conf[3, 1], d=conf[4, 1])
        upper = list(a=conf[1, 2], b=conf[2, 2], c=conf[3, 2], d=conf[4, 2])
    }
        
    #fitted values
    y    <- model$m$lhs()
    yhat <- as.vector(fitted(model))

    if( !is.null(conf) ){
        ylower <- sigmoid(pars=lower, xx=x)
        yupper <- sigmoid(pars=upper, xx=x)
    } else {
        ylower <- rep(0, length(y))
        yupper <- rep(0, length(y))
    }
        
    #plot(y ~ x, ylim=c(ylow,yup), xlim=c(0,max(x)),
    #           main=alg, xlab="log2(Fe)", ylab="Fraction of Enriched Communities")
    #lines(x, yi,     lty = 1, lwd = 2, col = "black")   
    #lines(x, yhat,   lty = 2, lwd = 2, col = "red")   

    #if( !is.null(conf) ){
    #    lines(x, ylower, lty = 2, lwd = 1, col = "blue")
    #    lines(x, yupper, lty = 2, lwd = 1, col = "blue")
    #}

    df <- cbind( rep(alg,length(x)), x, y, yhat, ylower, yupper )
    
    R = length(rates)
    Rsize = rep(1,R)
    Rcol  = rep("grey",R)
    indx  = which(rates==-2)
    if( length(indx) != 0 ){
        Rsize[indx[1]] = 2
        Rcol[indx[1]]  = "black"
    }
    for( r in 1:R ){
        pp = list(a=0, b=1, c=rates[r], d=round(median(x)) )
        yi <-  sigmoid(pars=pp, xx=x)
        df <- cbind(df,yi)
    }

    #df <- cbind( rep(alg,length(x)), x, y, yhat, yi, ylower, yupper )
    colnames(df) <- c("alg", "x", "y", "yhat", "ylower", "yupper", sprintf("yiR%f", seq(1,R,1)))
    df <- as.data.frame(df)

    #--- labels
    qq   = as.vector(quantile(as.numeric(df$x)))
    xval = qq
    xlab = as.character(qq)

    ylow <- ifelse( min(y) < 0, min(y), 0)
    yup  <- ifelse( max(y) > 1, max(y), 1)
    #---

    plotCI=TRUE
    if( is.null(conf) ){ plotCI = FALSE }
    
    pv = as.numeric(pv)
    if( is.na(pv) ){ pv = 0 }
    
    gplot <- ggplot(df, aes(as.numeric(df$x)))+
        geom_point(aes(y=as.numeric(as.vector(df$y))),   shape=1, size=2.5)+
        geom_line(aes(y=as.numeric(as.vector(df$yhat))), linetype="dashed", color="red", size=2)+
        geom_line(aes(y=as.numeric(as.vector(df$yiR1))), linetype="solid", color=Rcol[1], size=Rsize[1])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR2))), linetype="solid", color=Rcol[2], size=Rsize[2])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR3))), linetype="solid", color=Rcol[3], size=Rsize[3])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR4))), linetype="solid", color=Rcol[4], size=Rsize[4])+
        geom_line(aes(y=as.numeric(as.vector(df$yiR4))), linetype="solid", color=Rcol[5], size=Rsize[5])+
        {if(plotCI)geom_line(aes(y=as.numeric(as.vector(df$ylower))), linetype="dashed", color="blue", size=2)}+
        {if(plotCI)geom_line(aes(y=as.numeric(as.vector(df$yupper))), linetype="dashed", color="blue", size=2)}+
        labs(x="log2(Fe)",y="Fraction of Enriched Communities",title=sprintf("%s, KS p.value = %3.e", alg, pv))+
        theme(axis.title.x=element_text(face="bold",size=rel(1.5)),
              axis.title.y=element_text(face="bold",size=rel(1.5)),
              legend.text=element_text(face="bold",size=rel(1.5)),
              plot.title=element_text(face="bold",size=rel(1.5)),
              legend.position="bottom")+
        coord_cartesian(xlim = c(min(x), max(x)), ylim=c(ylow, yup))+
        #scale_y_continuous(expand=c(0,0),limits=c(ylow,yup))+
        #scale_x_discrete(expand=c(0,0), limit=xval, labels=xlab)+
        theme(panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_line(colour="grey40",size=0.1),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(linetype="solid",fill=NA))+
        guides(color = FALSE,
               alpha = FALSE,
               size  = FALSE)
    
    return(list(gplot=gplot, df=df))
    
}

errorY <- function( x, N=100, SD=0.05){
    x = x+rnorm(N, 0, SD)
    return(sd(x)/sqrt(length(x)))
    
}

sigmaSqY <- function( x, N=100, SD=0.05){

    x = x+rnorm(N, 0, SD)

    N  = length(x)
    mu = mean(x)

    return( (1/N) * sum( (x-mu)^2 ) )
    
}

addNoise <- function( Y, MN=0, SD=0.05 ){
    return( Y+rnorm(length(Y),mean=MN, sd=SD) )
}

#googness of fit tests
gofs <- function(x, rate, model, sigma2=NULL, countDATA=TRUE ){

    y    = model$m$lhs()
    yhat = fitted(model)

    R  = length(rate)
    KS = list()
    
    for( r in 1:R ){
    
        pp = list(a=0, b=1, c=rate[r], d=round(median(x)) )
        yi = sigmoid(pars=pp, xx=x)

        KS[[r]]      = ks.test(yhat, yi)
        names(KS)[r] = sprintf("rate_%f", rate[r]) 
    }

    return(KS)
    
    #N = length(y)
    #m = length(models[[1]]$m$getPars())
    #v = (N - m)
    
    #if( is.null(sigma2 ) ){
    #    if( countDATA ){
    #        sigma2 = abs(yhat) #use Poission error for count data        
    #    } else {
    #        sigma2 = model$weights
    #    }
    #}

    #KS = ks.test(yhat, yi)
    #return(list(P=KS$p.value))
    
    #chi2 = sum( (y-yhat)^2/sigma2 )

    #aa = (mean(abs(y-yhat))/mean(y)) * 100
    
    #aa = sum( (y-yhat)^2/yi )
    #aa = ( sum((y-yhat)^2) / sum(y^2) )^0.5
    
    #P    = chi2^(0.5*(v-2)) * exp( -0.5*chi2 )

    #if( P > 1 ) { P = 1 }
    #if( P < 0 ) { P = 0 }
    
    #chi2 = sum( (y-yhat)^2/sigma2 )

    #chi2x=(1/(N-m)) *var(model$m$resid())/var(y)
    
    #return( list(Chi2=chi2, RChi2=(1/(N-m))*chi2, Chi2x=chi2x) )
    #return(list(Chi2=chi2, v=v, P=P))
    #return(list(aa=aa))
}

calIC <- function(model, type="AICc"){

    val = AICc( model )
    
    if( type == "AIC" ){
        val = AIC( model )
    }

    if( type == "BIC" ){
        val = BIC( model )
    }

    return(val)
    
}

#---OUT Dir
OUT    <- vector(length=4)
OUT[1] <- DIRS[grepl("EnrichmentPackage",DIRS)]
OUT[2] <- DIRS[grepl("Graphs",DIRS)]
OUT[3] <- DIRS[grepl("parameterFiles",DIRS)]
OUT[4] <- DIRS[grepl("Annotations",DIRS)]

plotDIR <- sprintf("%s/PLOTS",OUT[1])
if( !file_test("-d",plotDIR) ){
    dir.create(plotDIR)
}
#---

tt <- read.delim("test.data.csv",sep="\t", header=T, check.names=F)

N = length(colnames(tt))
x = as.numeric(colnames(tt)[3:N])

SDv   = c(0, 0.01, 0.05, 0.1)
SDlab = c("0","0.01","0.05","0.1")
s     = 4

for( s in 1:length(SDv) ){

    models <- list()
    
    for( i in 1:length(tt[,1]) ){

        y = as.numeric(tt[i,3:N])/as.numeric(tt[i,2])

        y = addNoise(y, SD=SDv[s])
    
   #vary = unlist(lapply(y, 2, FUN=sigmaSqY)) #variance_y or sigma_squared_y
    
    #pp = list(a=0, b=1, c=-1, d=round(median(x)) )
        pp = list(a=0, b=round(max(y)), c=-2, d=round(median(x)) )
    
    #df = data.frame(x,y)
    
    #fits 
    #m.s <- nlsLM(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = pp, trace = FALSE, weights = minpack.lm::wfct(1/fitted^2), data = df)
    #m.s <- nlsLM(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = pp, trace = FALSE, weights = 1/vary, data = df)
        m.s <- nlsLM(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = pp, trace = FALSE)
    
    #fit adding noise
    ##simDNoisy <- sigmoid(pars=pp,xx=x) + rnorm(length(x), 0, 0.05)
    #simDNoisy <- y + rnorm(length(x), 0, 0.05)
    #m.s       <- nlsLM(simDNoisy ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = pp, trace = FALSE)
        
    
    ### perform fit
    ##nls.out <- nls.lm(par=pp, fn = residFun, observed = simDNoisy,
    ##              xx = x, control = nls.lm.control(nprint=1))

    ##plot(x,simDNoisy, main="data")
    ##lines(x,sigmoid(as.list(coef(nls.out)), x), col=2, lwd=2)
    
    models[[i]]      <- m.s
    names(models)[i] <- as.character(tt[i,1]) 

    rm(m.s)
    
}


#gofp = list(a=0, b=1, c=-2, d=round(median(x)) )
rates = c(-4, -3, -2, -1, -0.5)
    
GPLOTS <- list()
gof    <- list()

for( i in 1:length(names(models)) ){

    ks  <- gofs(x, rates, models[[i]])
    tmp <- plotSigmoid( x=x, rates=rates, model=models[[i]], alg=names(models)[i], pv=ks[[3]]$p.value)

    gof[[i]]    <- ks
    GPLOTS[[i]] <- tmp$gplot

    names(gof)[i]    <- names(models)[i]
    names(GPLOTS)[i] <- names(models)[i]
    
}

p <- plot_grid(plotlist=GPLOTS, labels = "AUTO", label_size=50, label_fontface="bold")
#ggsave("Plots.png", p, width=20, height=20, device="png")
ggsave(sprintf("%s/Fitting_%s.png",plotDIR,SDlab[s]), p, width=20, height=20, device="png")

rm(models, GPLOTS, gof)

}
