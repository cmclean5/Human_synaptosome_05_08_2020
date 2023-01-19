rm(list=ls())

#---WIDTH and HEIGHT for plots (in px)
WIDTH=480
HEIGHT=480

#---WIDTH and HEIGHT for plots (in inches)
WIDTHin=6.4
HEIGHTin=6.4

#---WIDTH and HEIGHT for plots (in cm)
WIDTHcm=12.7
HEIGHTcm=12.7

#---WIDTH and HEIGHT for plots (in mm)
WIDTHmm=127
HEIGHTmm=127


#set required R libraries
#source("http://www.bioconductor.org/biocLite.R")
library(igraph);
library(lattice);
library(ggplot2);
library(stringr);
library(gtable);
library(grid);
library(ggrepel);
library(scales);
library(plyr);
library(biomaRt);
library(knitr);

##---For Bonferroni correction
alpha <- vector(length=3)
alpha[1] <- 0.05
alpha[2] <- 0.01
alpha[3] <- 0.001

##---For Bonferroni correction
stars    <- vector(length=3)
stars[1] <- "*"
stars[2] <- "**"
stars[3] <- "***"

##--- adding vertex properties
COLLAPSE <- vector(length=2)
COLLAPSE[1] <- ";"
COLLAPSE[2] <- "&"
c=1
