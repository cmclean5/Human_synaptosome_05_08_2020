source('../setUp.R')
library(corrplot)
library(RColorBrewer)

dd1 = readRDS("../POWERlawFIT/originalNtwrkMeasures.rds")

dd2 = read.delim("../EntropyRate/PPI_PSP/SignalEntropyRate.csv",sep="\t",header=T)

dd3=read.delim("../Consensus/PPI_PSP/PPI_PSP_Measures.csv",sep="\t",header=T)

indx = which(names(dd1)==subDIR[S])

dd1 = dd1[[indx]]


df=cbind(dd2[,c(4,5)],dd3[,c(3,4,5,6,8,9,10,12)])
#colnames(df) = c('SR_UP','SR_DOWN','DEGREE','CLOSENESS','BETWEENNESS',
#                 'CLUSTERING_COEFFICIENT','SEMI_LOCAL','SHORTEST_PATH',
#                 'PAGE_RANK','BRIDGENESS')
colnames(df) = c('SR_UP','SR_DOWN','DEG','CLN','BET','CC','CL','SP','PR','BR')


cor_matrix <- cor(df[,1:ncol(df)], use = 'complete.obs')

#pt = corrplot(cor_matrix, method="circle",type="upper", order="hclust")


png(filename = "PLOTS/corrplot.png", width=WIDTH, height=HEIGHT, unit="px");
corrplot(cor_matrix, method="circle",type="upper", order="AOE", tl.col = "black", tl.srt = 45, col = brewer.pal(n = 8, name = "PuOr"))
dev.off();

#pt = corrplot.mixed(cor_matrix, lower = "circle", upper = "number", tl.pos = "lt", diag = "u")
