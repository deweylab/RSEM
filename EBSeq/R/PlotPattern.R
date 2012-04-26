PlotPattern<-function(Patterns){
	par(oma=c(3,3,3,3))
	PatternCol=rainbow(ncol(Patterns))
	heatmap(Patterns,col=PatternCol,Colv=NA,Rowv=NA,scale="none")

}

