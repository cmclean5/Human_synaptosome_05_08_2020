rm(list=ls())

require("minpack.lm")
require("qpcR")


### Sigmoid function ### create a function to generate sigmoid pattern
sigmoid <- function(x, lower_asymptote, carrying_capacity, growth_rate, time_max) {
    return(lower_asymptote + ((carrying_capacity - lower_asymptote)/(1 + exp(-growth_rate * 
        (x - time_max)))))
}


tt <- read.delim("test.data.csv",sep="\t", header=T, check.names=F)
N = length(tt)
x = as.numeric(colnames(tt)[3:N])

models <- list()
aics   <- list()
plots  <- list()

for( i in 1:length(tt[,1]) ){

    y = as.numeric(tt[i,3:N])/as.numeric(tt[i,2])



#x = seq(0,10,0.1)
#y <- sigmoid(x, 0, 1, -1, 5) + rnorm(length(x), 0, 0.05)


    m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), start = list(a = 0, 
               b = 1, c = -1, d = round(median(x))), trace = FALSE)#TRUE)


    m.s <- nls(y ~ a + ((b - a)/(1 + exp(-c * (x - d)))), par = list(a = 0, 
               b = 1, c = -1, d = round(median(x))), trace = FALSE)#TRUE)

    plot(y ~ x)
    lines(x, fitted(m.s), lty = 2, lwd = 2, col = "red")
    # adding the 95% confidence interval around the fitted coefficient
    conf <- confint(m.s)

    lines(x, sigmoid(x, conf[1, 1], conf[2, 1], conf[3, 1], conf[4, 1]), lty = 2, 
          lwd = 1, col = "blue")
    lines(x, sigmoid(x, conf[1, 2], conf[2, 2], conf[3, 2], conf[4, 2]), lty = 2, 
          lwd = 1, col = "blue")


    models[[i]]   <- m.s
    names(models) <- as.character(tt[i,1]) 

        
    # goodness of fit 
    aic         <- AICc(m.s)
    aics[[i]]   <- aics
    names(aics) <- as.character(tt[i,1]) 
    
}
