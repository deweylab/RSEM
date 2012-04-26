GetMultiPP <- function(EBout){
	PP=EBout$PPDE	
	MAP=colnames(EBout$f)[apply(EBout$f,1,which.max)]
	AllParti=EBout$AllParti
	out=list(PP=PP, MAP=MAP,Patterns=AllParti)
}
