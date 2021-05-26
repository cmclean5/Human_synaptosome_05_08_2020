source('setUp.R')

# magenta50per = t_col("FF00FF",percent=50)
t_col <- function(color, percent = 50, name = NULL) {
   ## https://www.dataanalytics.org.uk/make-transparent-colors-in-r/
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}
## END

#---OUT Dir
OUT    <- vector(length=3)
OUT[1] <- DIRS[grepl("GeneSets",DIRS)]
OUT[2] <- DIRS[grepl("Clustering",DIRS)]
OUT[3] <- DIRS[grepl("Graphs",DIRS)]

grdir <- sprintf("%s/%s",OUT[3],subDIR[S])
if( !file_test("-d",grdir) ){
    dir.create(grdir)
}
#---


#---READ IN GRAPH 
gg  <- igraph::read.graph(sprintf("%s/%s.gml",grdir,subDIR[S]),format="gml")
ids <- V(gg)$name;

LOWcol="#FFFFCC"
MIDcol="#00EEEE"
HIGHcol="#436EEE"
CorePSD95col="#FF00FF"

color_range <- colorRampPalette(c(LOWcol,MIDcol,HIGHcol))

Cprob=V(gg)$probSpectral

CprobTB=table(Cprob)

TBnames=names(CprobTB)

TBnames=TBnames[order(TBnames)]

col_table = cbind(TBnames,color_range(length(TBnames)))

val = col_table[match(V(gg)$probSpectral,col_table[,1]),2]

val = ifelse(V(gg)$CorePSD95==1,CorePSD95col,val)

gg = set.vertex.attribute(gg,"color",V(gg),val)


##---Write .gml graph to file
igraph::write.graph(gg, sprintf("%s/%s.gml",grdir,subDIR[S]), "gml")
##---Write .graphml graph to file
igraph::write.graph(gg, sprintf("%s/%s.graphml",grdir,subDIR[S]), "graphml")

}
