DenNHistTable <-
function(QList,Alpha,Beta,AList="F")
{	
	par(mfrow=c(3,4))
	plot(1, type="n", axes=F, xlab="", ylab="", main="No 3' end  No 5' end",cex.main=1)
	plot(1, type="n", axes=F, xlab="", ylab="",main="With 3' end No 5' end",cex.main=1)
	plot(1, type="n", axes=F, xlab="", ylab="",main="With 5' end No 3' end",cex.main=1)
	for (i in c(1,2,4,6,8)){
		 alpha.use=Alpha
	hist(QList[[i]][QList[[i]]<.98&QList[[i]]>0],prob=T,col="blue",breaks=100,main=ifelse(i==1,"With 5' end With 3' end",""),cex.main=1, xlim=c(0,1),xlab=paste("Q alpha=",round(alpha.use,2)," beta=",round(Beta[i],2),sep=""))
	if(i==1)mtext("Ng=1",side=4, cex=1)
	if(i==8)mtext("Ng=2", side=4,cex=1)
	tmpSize=length(QList[[i]][QList[[i]]<.98])

        tmpseq=seq(0.001,1,length=1000)
	ll=tmpseq
		 lines(ll,dbeta(ll,alpha.use,Beta[i]),col="green",lwd=2)
	legend("topright",c("Data","Fitted density"),col=c("blue","green"),lwd=2,cex=.5)
}
	
	for (i in c(3,5,7,9)){
		 alpha.use=Alpha
	hist(QList[[i]][QList[[i]]<.98&QList[[i]]>0],prob=T,col="blue",breaks=100,main=ifelse(i==1,"With 5' end With 3' end exons",""),xlim=c(0,1),xlab=paste("Q alpha=",round(alpha.use,2)," beta=",round(Beta[i],2),sep=""))
	if(i==9)mtext("Ng=3", side=4,cex=1)

	tmpSize=length(QList[[i]][QList[[i]]<.98])

        tmpseq=seq(0.001,1,length=1000)
	ll=tmpseq
		 lines(ll,dbeta(ll,alpha.use,Beta[i]),col="green",lwd=2)
	legend("topright",c("Data","Fitted density"),col=c("blue","green"),cex=.5, lwd=2)
}




	}

