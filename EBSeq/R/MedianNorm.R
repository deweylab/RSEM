MedianNorm=function(Data){

    geomeans <- exp(rowMeans(log(Data)))
	apply(Data, 2, function(cnts) median((cnts/geomeans)[geomeans >  0]))
}
