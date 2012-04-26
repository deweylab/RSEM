TPFDRplot <-
function(DESeqP, EBZ, TrueDE, main, FDR=NULL){
	Seq=seq(0.001,0.5,by=0.001)
	DETPR=rep(0,length(Seq))
	EBTPR=rep(0,length(Seq))
	DEFDR=rep(0,length(Seq))
	EBFDR=rep(0,length(Seq))
	DETPNum=rep(0,length(Seq))
    EBTPNum=rep(0,length(Seq))
    DEFDNum=rep(0,length(Seq))
    EBFDNum=rep(0,length(Seq))
	for (i in 1:length(Seq)){
		DESeqOnes=names(DESeqP)[DESeqP<=Seq[i]]
		if (length(FDR)==0) EBOnes=names(EBZ)[EBZ>=crit.fun(1-EBZ, Seq[i])]
		else if (FDR=="H") EBOnes=names(EBZ)[EBZ>=(1-Seq[i])]
			else EBOnes=names(EBZ)[EBZ>=FDR[i]]

		DETPNum[i]=sum(DESeqOnes%in%TrueDE)
		EBTPNum[i]=sum(EBOnes%in%TrueDE)
		DEFDNum[i]=sum(!DESeqOnes%in%TrueDE)
		EBFDNum[i]=sum(!EBOnes%in%TrueDE)
		
		DETPR[i]=DETPNum[i]/length(TrueDE)
		EBTPR[i]=EBTPNum[i]/length(TrueDE)
		DEFDR[i]=DEFDNum[i]/length(TrueDE)
		EBFDR[i]=EBFDNum[i]/length(TrueDE)
	}
	plot(Seq,DETPR,ylim=c(0,1),xlim=c(0,.5),type="l",col="red", main=paste(main, "TPR"),xlab="controled FDR level", ylab="TPR",lwd=2)
	lines(Seq,EBTPR,col="blue",lwd=2)
	legend("bottomright",lwd=2, col=c("red","blue"), c("DESeq","EBSeq"))

	plot(Seq,DEFDR,ylim=c(0,1),xlim=c(0,.5),type="l",col="red", main=paste(main, "FDR"),xlab="controled FDR level", ylab="TPR",lwd=2)
	lines(Seq,EBFDR,col="blue",lwd=2)
	legend("topleft", lwd=2, col=c("red","blue"), c("DESeq","EBSeq"))


	output=cbind( DETPR,EBTPR, DEFDR,EBFDR,DETPNum,EBTPNum,DEFDNum,EBFDNum)
}

