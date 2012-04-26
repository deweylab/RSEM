crit_fun<-function (PPEE, thre) 
{
    y <- cumsum(sort(PPEE))/(1:length(PPEE))
    mm <- y < thre
    index <- sum(mm)
    if (index > 0) {
        out <- 1 - sort(PPEE)[index]
	    }		
    if (index == 0) {
		        out <- 1
				    }
    names(out) <- NULL
    return(out)
}

