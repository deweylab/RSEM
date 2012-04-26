PolyFitPlot <-
function(X , Y , nterms , xname="Estimated Mean", yname="Estimated Var", pdfname="", xlim=c(-1,5), ylim=c(-1,7), ChangeXY=F,col="red"){
	
	b=rep(NA,nterms)
	logX=matrix(rep(X, nterms),ncol=nterms, byrow=T)
	for (i in 1:nterms)
		logX[,i]=(log10(X))^i
	colnames(logX)=paste("logmu^",c(1:nterms))
	rownames(logX)=names(X)
	NotUse=c(names(X)[X==0],names(Y)[Y==0],names(X)[rowMeans(logX)==-Inf],names(X)[rowMeans(logX)==Inf])
	Use=names(X[!names(X)%in%NotUse])
	Lm=lm(log10(Y[Use])~logX[Use,1:nterms])
	b=summary(Lm)$coefficients[2:(nterms+1),1]
	d=summary(Lm)$coefficients[1,1]
	bvec=matrix(rep(b,length(X)),ncol=nterms,byrow=T)
	fit=rowSums(logX*bvec)+d
	main2=NULL
	if (ChangeXY==T){
		X.plot=log10(Y)
		Y.plot=log10(X)
		fit.X.plot=fit
		fit.Y.plot=log10(X)
	}
	else{
        X.plot=log10(X)
        Y.plot=log10(Y)
	    fit.X.plot=log10(X)
		fit.Y.plot=fit
				   }

	for (i in 1:nterms)
		main2=paste(main2,round(b[i],2),"*log(",xname,")^",i,"+")
	main=pdfname
	
	smoothScatter(X.plot, Y.plot ,main=main,xlim=xlim,ylim=ylim,xlab=xname,ylab=yname,axes=F)
	axis(1,at=seq(xlim[1],xlim[2],by=1), 10^seq(xlim[1],xlim[2],by=1))
	axis(2,at=seq(ylim[1],ylim[2],by=2), 10^seq(ylim[1],ylim[2],by=2))
	Sortit=order(fit.X.plot)
	lines(fit.X.plot[Sortit],fit.Y.plot[Sortit],col=col,lwd=3)
	output=list(b=b,d=d,lm=Lm,fit=fit,sort=Sortit)
	names(output$b)=paste(xname,"^",c(1:length(output$b)))
	output
}

