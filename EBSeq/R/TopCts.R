TopCts <-
function(pvalue, PP=NULL, TrueNames, TopNum){
	NumOfMethods=ncol(pvalue)
	puse=pvalue
	if(1%in%PP)puse[,PP==1]=1-pvalue[,PP==1]
	#puse.list=data.frame(puse)
	FD=matrix(rep(0,NumOfMethods*TopNum),ncol=NumOfMethods)
#	Rank=apply(puse,2,rank)
#	for(i in 1:TopNum)
#		FD[i,]=sapply(1:NumOfMethods, function(j)sum(!rownames(Rank)[Rank[,j]<=i]%in%TrueNames))	
#	FD=sapply(1:TopNum, function(i)sapply(1:NumOfMethods, function(j)sum(!rownames(Rank)[Rank[,j]<=i]%in%TrueNames)))
	for (s in 1:NumOfMethods){
		tmp=puse[,s]
		names(tmp)=rownames(puse)
		sorttmp=sort(tmp)
		for( c in 2:TopNum)
			FD[c, s]=FD[(c-1),s]+as.numeric(!names(sorttmp)[c]%in%TrueNames)
	}
	FD
	#matplot(TopNum,FD,type="l",ylim=c(0,1),xlab="Top DE selected", ylab="FDR")
	#legend("rightbottom",col=1:TopNum, lty=1:TopNum, names)
	}

