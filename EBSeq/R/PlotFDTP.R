PlotFDTP <-
function(TopNum, FDR, TPR,names)
{
  
  matplot(FDR, TPR, xlim=c(0,.5), ylim=c(0,1) ,type="l",lwd=2,xlab="FDR", ylab="TPR")
    legend("bottomright",col=1:TopNum, lwd=2, lty=1:TopNum, names)


}

