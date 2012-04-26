QQP <-
function(QList,AlphaResult,BetaResult,name,AList="F",GroupName){
	
		    for (i in 1:length(BetaResult)){
				tmpSize=length(QList[[i]][QList[[i]]<1 & !is.na(QList[[i]])])
			if (AList=="F") rdpts=rbeta(tmpSize,AlphaResult,BetaResult[i])
				else rdpts=rbeta(tmpSize,AlphaResult[i],BetaResult[i])
	qqplot(QList[[i]][QList[[i]]<1], rdpts,xlab="estimated q's", ylab="simulated q's from fitted beta",main=paste(name,GroupName[i],sep=" "),xlim=c(0,1),ylim=c(0,1))
	fit=lm(sort(rdpts)~sort(QList[[i]][QList[[i]]<1  & !is.na(QList[[i]])]))
	abline(fit,col="red")
	
			}
}

