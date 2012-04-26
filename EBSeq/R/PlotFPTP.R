PlotFPTP <-
function(TopNum, FPR, TPR,names)
{
	 
	  matplot(FPR, TPR,xlim=c(0,.1), ylim=c(0,1) ,type="l",lwd=2, xlab="FPR", ylab="TPR")
	      legend("bottomright",col=1:TopNum,lwd=2, lty=1:TopNum, names)


}

