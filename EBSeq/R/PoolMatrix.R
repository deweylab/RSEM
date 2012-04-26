PoolMatrix <-
function(Data,reads,type)
{
poolnames=names(Data)
poolM=NULL
for (po in 1:8)
	poolM=cbind(poolM,Data[[po]][,1])
rownames(poolM)=rownames(Data[[1]])
colnames(poolM)=poolnames

#poolValue=poolM*reads
poolValue=poolM
for (col in 1:8)
	poolValue[,col]=poolM[,col]*reads[col]
poolValue=round(poolValue)
if (type=="G")
	{
		poolM=cbind(Data[[1]][,2],poolM)
		poolValue=cbind(Data[[1]][,2],poolValue)
		colnames(poolM)=c("Groups",poolnames)
		colnames(poolValue)=c("Groups",poolnames)
	}
poolOutput=list(poolM=poolM,poolValue=poolValue)
}

