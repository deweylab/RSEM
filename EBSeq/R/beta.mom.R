beta.mom <-
function(qs.in){
	xbar<-mean(qs.in)
	s2<-var(qs.in)
	term<-(xbar*(1-xbar))/s2
	alpha.hat<-xbar*(term-1)
	beta.hat<-(1-xbar)*(term-1)
	return(c(alpha.hat,beta.hat))
}

