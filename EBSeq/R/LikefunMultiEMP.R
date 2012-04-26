LikefunMulti <-
function(ParamPool, InputPool)
{

NoneZeroLength=InputPool[[4]]
AlphaIn=ParamPool[1]
BetaIn=ParamPool[2:(1+NoneZeroLength)]
PIn=ParamPool[(2+NoneZeroLength):length(ParamPool)]
PInAll=c(1-sum(PIn),PIn)
ZIn=InputPool[[3]]
Input=InputPool[[2]]
InputSP=InputPool[[1]]
RIn=InputPool[[5]]
RInSP=InputPool[[6]]
NumIn=InputPool[[7]]
AllParti=InputPool[[8]]
PInMat=matrix(rep(1,nrow(Input)),ncol=1)%*%matrix(PInAll,nrow=1)
##Function here
FList=sapply(1:nrow(AllParti),function(i)sapply(1:nlevels(as.factor(AllParti[i,])),
			                        function(j)f0(do.call(cbind,InputSP[AllParti[i,]==j]),AlphaIn, BetaIn, 
	                                do.call(cbind,RInSP[AllParti[i,]==j]), NumIn, log=T)),
			                        simplify=F) 
FPartiLog=sapply(FList,rowSums)
#FMat=exp(FPartiLog)
FMat=FPartiLog
-sum(ZIn*(FMat+log(PInMat)))
}

