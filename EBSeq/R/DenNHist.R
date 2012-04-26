DenNHist <-
function(QList,Alpha,Beta,name,AList="F",GroupName)
{
    if(!is.list(QList)) QList=list(QList) 	
	for (i in 1:length(QList)){
		if (AList=="F") alpha.use=Alpha
			if(AList=="T")  alpha.use=Alpha[i]
	hist(QList[[i]][QList[[i]]<.98&QList[[i]]>0],prob=T,col="blue",breaks=100,main=paste(GroupName[i],name,sep=" "),xlim=c(0,1),xlab=paste("Q alpha=",round(alpha.use,2)," beta=",round(Beta[i],2),sep=""))
	tmpSize=length(QList[[i]][QList[[i]]<.98])
        tmpseq=seq(0.001,1,length=1000)
        #tmpdensity=dbeta(tmpseq,AlphaResult,BetaResult[i])
        #points(tmpseq,tmpdensity, type="l",col="green")
	#ll=dbeta(tmpseq,Alpha,Beta[i])
	ll=tmpseq
		 lines(ll,dbeta(ll,alpha.use,Beta[i]),col="green",lwd=2)
	legend("topright",c("Data","Fitted density"),col=c("blue","green"),lwd=2)
}
	
	}

