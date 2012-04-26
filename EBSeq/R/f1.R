f1 <-
function(Input1, Input2, AlphaIn, BetaIn, EmpiricalRSP1,EmpiricalRSP2,NumOfGroup, log){
	F0.1=f0(Input1, AlphaIn, BetaIn, EmpiricalRSP1, NumOfGroup, log=log)
	F0.2=f0(Input2, AlphaIn, BetaIn, EmpiricalRSP2, NumOfGroup, log=log)
	
	if (log==F) Result=F0.1*F0.2
	if (log==T) Result=F0.1+F0.2
	Result
}

