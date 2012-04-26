Likefun <-
function(ParamPool, InputPool)
{

NoneZeroLength=InputPool[[5]]
AlphaIn=ParamPool[1]
BetaIn=ParamPool[2:(1+NoneZeroLength)]
PIn=ParamPool[2+NoneZeroLength]
ZIn=InputPool[[4]]
Input=InputPool[[3]]
Input1=matrix(InputPool[[1]],nrow=nrow(Input))
Input2=matrix(InputPool[[2]],nrow=nrow(Input))
RIn=InputPool[[6]]
RInSP1=matrix(InputPool[[7]],nrow=nrow(Input))
RInSP2=matrix(InputPool[[8]],nrow=nrow(Input))
NumIn=InputPool[[9]]
##Function here
#LikelihoodFunction<- function(NoneZeroLength){
	F0=f0(Input, AlphaIn, BetaIn, RIn, NumIn, log=T)
	F1=f1(Input1, Input2, AlphaIn, BetaIn, RInSP1,RInSP2, NumIn, log=T)
		F0[F0==Inf]=min(!is.na(F0[F0!=Inf]))
		F1[F1==Inf]=min(!is.na(F1[F1!=Inf]))

	-sum((1-ZIn)*F0+ (1-ZIn)* log(1-PIn) + ZIn*F1 + ZIn*log(PIn))
}

