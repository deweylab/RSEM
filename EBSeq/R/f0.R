f0 <-
function(Input, AlphaIn, BetaIn, EmpiricalR, NumOfGroups, log)
{	
		 
		BetaVect=do.call(c,sapply(1:length(BetaIn),function(i)rep(BetaIn[i],NumOfGroups[i]),simplify=F))
		SampleNum=dim(Input)[2]
		#Product part
		ChooseParam1=round(Input+EmpiricalR-1)
		roundInput=round(Input)
		EachChoose=sapply(1:SampleNum, function(i)lchoose(ChooseParam1[,i], roundInput[,i]))
		
		SumEachIso=rowSums(Input)
		param1=AlphaIn + rowSums(EmpiricalR)
		param2=BetaVect + SumEachIso
		LogConst=rowSums(EachChoose)+lbeta(param1, param2)-lbeta(AlphaIn, BetaVect)


		if (log==F) FinalResult=exp(LogConst)
		if (log==T) FinalResult=LogConst
    FinalResult
}

