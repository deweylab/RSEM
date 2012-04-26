
QuantileNorm=function(Data, Quantile){
	#SortData=apply(Data, 2, sort)
	QtilePt=apply(Data, 2, function(i)quantile(i, Quantile))
	Size= QtilePt * prod(QtilePt) ^ (-1/ncol(Data))
	Size
	}

