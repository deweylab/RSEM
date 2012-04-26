CheckNg<-function(NewMean, NewVar,nterm, xlim, ylim){
	Ng=1=PolyFit_ENAR(NewMean[[1]],NewVar[[1]],nterm,"Mean","Variance","Ng=1",xlim, ylim)
	sortNg1=order(NewMean[[1]])
	Ng=2=PolyFit_ENAR(unlist(NewMean[c(2,4,6,8)]),unlist(NewVar[c(2,4,6,8)]),nterm,"Mean","Variance","Ng=2",xlim, ylim)
	sortNg2=order(unlist(NewMean[c(2,4,6,8)]))
	Ng=3=PolyFit_ENAR(unlist(NewMean[c(3,5,7,9)]),unlist(NewVar[c(3,5,7,9)]),nterm,"Mean","Variance","Ng=3",xlim, ylim)
	sortNg3=order(unlist(NewMean[c(3,5,7,9)]))

	ALL=PolyFit_ENAR(unlist(NewMean),unlist(NewVar),nterm,"Mean","Variance","",xlim, ylim)
	lines(log10(unlist(NewMean[c(2,4,6,8)]))[sortNg2],Ng=2$fit[sortNg2],col="green",lwd=2)
	lines(log10(unlist(NewMean[c(3,5,7,9)]))[sortNg3],Ng=3$fit[sortNg3],col="orange",lwd=2)
	lines(log10(unlist(NewMean[1]))[sortNg1],Ng=1$fit[sortNg1],col="pink",lwd=2)
	legend("topleft",col=c("red","pink","green","orange"),c("all","Ng=1","Ng=2","Ng=3"),lwd=2)
}







