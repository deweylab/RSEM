PlotTopCts <-
function(TopNum, FD, names)
{
    matplot(c(1:TopNum) , FD,type="l",xlab="Top DE selected", lwd=2, log="y", ylab="FD")
    legend("topleft",col=1:TopNum, lwd=2, lty=1:TopNum, names)

}

