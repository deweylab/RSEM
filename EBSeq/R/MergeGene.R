MergeGene <-
function(GeneSIMout, Num, Path="./"){
NumSample=ncol(GeneSIMout[[i]]$generateData)

NumGene=rep(0,Num)
for (i in 1:Num)NumGene[i]=nrow(GeneSIMout[[i]]$generateData)

MinNumGene=min(NumGene)
AproxNumDE=length(GeneSIMout[[1]]$TrueDE)
	
GeneMergeTable=matrix(rep(0,12),nrow=6)
	for(i in 1:Num)GeneMergeTable=GeneMergeTable+GeneSIMout[[i]][[1]]
	GeneMergeTable=GeneMergeTable/Num
	GeneMergeTable=round(GeneMergeTable,2)
	          
	GeneMergeDVD=rep(0,2)
	  for(i in 1:Num)GeneMergeDVD=GeneMergeDVD+GeneSIMout[[i]][[3]]
		  GeneMergeDVD=round(GeneMergeDVD/Num,2) 
				          
	  GeneMergePhi=matrix(rep(0,2),nrow=2)
		  for(i in 1:Num)GeneMergePhi=GeneMergePhi+GeneSIMout[[i]][[4]]
			  GeneMergePhi=round(GeneMergePhi/Num,2)
## Write
TXTname=paste(paste(Path,paste("Gene","DVD",GeneMergeDVD[1], GeneMergeDVD[2],"Phi",GeneMergePhi[1], GeneMergePhi[2],"Sample",NumSample,sep="_"),sep=""),".txt",sep="")
write.table(GeneMergeTable, file=TXTname)


####### Note everytime # DE genes and # total genes may different. (since NA issue)
  GeneMergeFD=matrix(rep(0,5*MinNumGene),ncol=5)
  GeneMergeFD.p=matrix(rep(0,5*MinNumGene),ncol=5)
  GeneMergeTP.p=matrix(rep(0,5*MinNumGene),ncol=5)
  GeneMergeFN.p=matrix(rep(0,5*MinNumGene),ncol=5)
  GeneMergeTN.p=matrix(rep(0,5*MinNumGene),ncol=5)

  GeneMergeFDR=matrix(rep(0,5*MinNumGene),ncol=5)
  GeneMergeTPR=matrix(rep(0,5*MinNumGene),ncol=5)
  GeneMergeFPR=matrix(rep(0,5*MinNumGene),ncol=5)


  for(i in 1:Num){
	# Make sure names in the same order
	# Get FD number for each number of genes found
    TotalNum=nrow(GeneSIMout[[i]]$generateData)
	NumDE=length(GeneSIMout[[i]]$TrueDE)
	EBSeqNames=names(GeneSIMout[[i]]$EBSeqPP)
    tmpMatrix=cbind(GeneSIMout[[i]]$DESeqP[EBSeqNames],GeneSIMout[[i]]$edgeRP[EBSeqNames], exp(GeneSIMout[[i]]$BaySeqPP[EBSeqNames,2]),GeneSIMout[[i]]$BBSeqP[EBSeqNames],GeneSIMout[[i]]$EBSeqPP)
	# Bayseq and EBseq are PP. Others are p value 
    tmpFD=TopCts(tmpMatrix, c(0,0,1,0,1), GeneSIMout[[i]]$TrueDE[GeneSIMout[[i]]$TrueDE%in%EBSeqNames], MinNumGene)
    # Get percentage for FP, TP, TN, FN!
	tmpFD.p=tmpFD/TotalNum
	# TP = Find - FD
	tmpTP.p=(c(1:MinNumGene)-tmpFD)/TotalNum
	# FN = TrueDE - TP
	tmpFN.p=NumDE/TotalNum - tmpTP.p
	# TN = TrueEE - FD
	tmpTN.p=(TotalNum-NumDE)/TotalNum - tmpFD.p
	
	
	tmpFDR=tmpFD.p/(tmpFD.p+tmpTP.p)
	tmpFPR=tmpFD.p/(tmpFD.p+tmpTN.p)
	tmpTPR=tmpTP.p/(tmpFN.p+tmpTP.p)
	GeneMergeFDR=GeneMergeFDR+tmpFDR
	GeneMergeTPR=GeneMergeTPR+tmpTPR
	GeneMergeFPR=GeneMergeFPR+tmpFPR

    GeneMergeFD.p=GeneMergeFD.p+tmpFD.p
	GeneMergeTP.p=GeneMergeTP.p+tmpTP.p
	GeneMergeFN.p=GeneMergeFN.p+tmpFN.p
	GeneMergeTN.p=GeneMergeTN.p+tmpTN.p

	GeneMergeFD=GeneMergeFD+tmpFD
 }   
  GeneMergeFD=GeneMergeFD/Num
  GeneMergeFD.p=GeneMergeFD.p/Num
  GeneMergeTP.p=GeneMergeTP.p/Num
  GeneMergeFN.p=GeneMergeFN.p/Num
  GeneMergeTN.p=GeneMergeTN.p/Num

  GeneMergeFDR=GeneMergeFDR/Num
  GeneMergeTPR=GeneMergeTPR/Num
  GeneMergeFPR=GeneMergeFPR/Num


PlotTopName=paste(paste(Path,paste("Top","Gene","DVD",GeneMergeDVD[1], GeneMergeDVD[2],"Phi",GeneMergePhi[1], GeneMergePhi[2],"Sample",NumSample, sep="_"),sep=""),".pdf",sep="")

TrueDELength=length(GeneSIMout[[i]]$TrueDE[GeneSIMout[[i]]$TrueDE%in%EBSeqNames])
pdf(PlotTopName)
  PlotTopCts(TrueDELength,GeneMergeFD[1:TrueDELength,],c("DESeq","edgeR","BaySeq","BBSeq","EBSeq"))
dev.off()


PlotFDName=paste(paste(Path,paste("FDTP","Gene","DVD",GeneMergeDVD[1], GeneMergeDVD[2],"Phi",GeneMergePhi[1], GeneMergePhi[2],"Sample",NumSample,sep="_"),sep=""),".pdf",sep="")
pdf(PlotFDName)
  PlotFDTP(MinNumGene,GeneMergeFDR, GeneMergeTPR, c("DESeq","edgeR","BaySeq","BBSeq","EBSeq"))
dev.off()

PlotFPName=paste(paste(Path,paste("FPRTP","Gene","DVD",GeneMergeDVD[1], GeneMergeDVD[2],"Phi",GeneMergePhi[1], GeneMergePhi[2],"Sample",NumSample,sep="_"),sep=""),".pdf",sep="")
pdf(PlotFPName)
  PlotFPTP(MinNumGene,GeneMergeFPR, GeneMergeTPR, c("DESeq","edgeR","BaySeq","BBSeq","EBSeq"))
  dev.off()


out=list(GeneMergeTable=GeneMergeTable, GeneMergeDVD=GeneMergeDVD, GeneMergePhi=GeneMergePhi, GeneMergeFD=GeneMergeFD)


}

