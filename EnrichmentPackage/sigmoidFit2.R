source('../setUp.R')

require("minpack.lm")
require("qpcR")
#require("ggplot2")

### Sigmoid function ### create a function to generate sigmoid pattern
sigmoid <- function(pars, xx){#, lower_asymptote, carrying_capacity, growth_rate, time_max) {

    a = as.numeric(pars[1])
    b = as.numeric(pars[2])
    c = as.numeric(pars[3])
    d = as.numeric(pars[4])

    return( a + ((b-a)/(1+exp(-c*(xx-d)))) )
}

residFun <- function(pars, observed, xx) observed - sigmoid(pars,xx)


plotSigmoid <- function( x, model, alg="" ){

    conf <- confint(model)

    # adding the 95% confidence interval around the fitted coefficient
    lower = list(a=conf[1, 1], b=conf[2, 1], c=conf[3, 1], d=conf[4, 1])
    upper = list(a=conf[1, 2], b=conf[2, 2], c=conf[3, 2], d=conf[4, 2])

    #fitted values
    y = model$m$lhs()

    yhat   <- as.vector(fitted(model))
    ylower <- sigmoid(pars=lower, xx=x)
    yupper <- sigmoid(pars=upper, xx=x)
    
    plot(y ~ x)
    lines(x, yhat,   lty = 2, lwd = 2, col = "red")   
    lines(x, ylower, lty = 2, lwd = 1, col = "blue")
    lines(x, yupper, lty = 2, lwd = 1, col = "blue")

    df <- cbind( rep(alg,length(x)), x, y, yhat, ylower, yupper )
    colnames(df) <- c("alg", "x", "y", "yhat", "ylower", "yupper")

    return(df)
    
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


tt <- read.delim("test.data.csv",sep="\t", header=T, check.names=F)

N = length(colnames(tt))
x = as.numeric(colnames(tt)[3:N])

models <- list()
plots  <- list()

for( i in 1:length(tt[,1]) ){

    y = as.numeric(tt[i,3:N])/as.numeric(tt[i,2])

    pp = list(a=0, b=1, c=-1, d=round(median(x)) )
    #pp = list(a=0, b=round(max(y)), c=-1, d=round(median(x)) )

    df = data.frame(x,y)
    
    #fit without noise
    m.s <- nlsLM(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = pp, trace = FALSE )#, weights = minpack.lm::wfct(1/fitted^2), data = df)

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

        
    # goodness of fit 
    #aic            <- AICc(m.s)
    #aics[[i]]      <- aic
    #names(aics)[i] <- as.character(tt[i,1]) 

    rm(m.s)
    
}


DF <- data.frame()
for( i in 1:length(tt[,1]) ){
    if(  models[[i]]$convInfo$stopCode != -1 ){
        tmp <- plotSigmoid( x=x, model=models[[i]], alg=names(models)[i] )
        DF <- rbind(DF, tmp)
    }
}




# goodness of fit
ic   <- list()
n    <- length(names(models))
for( i in 1:n ){
    ic[[i]] <- calIC( models[[i]] )
    names(ic)[i] <- names(models)[i]
}

oo1 <- matrix(NA, nrow=(n), ncol=(n))
colnames(oo1) <- names(ic)
rownames(oo1) <- names(ic)

for( i in 1:n ){
    for( j in i:n ){
        oo1[i,j] = evidence(ic[[i]],ic[[j]])
        oo1[j,i] = oo1[i,j]
    }
}

#--- plotting
#size=2
size=3
pt1 = plotMetric(oo1, "PLOTS", "test", "AIC", "AIC", valSIZE=size)

