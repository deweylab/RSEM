MergeIso <-
function(IsoSIMout, Num, Path="./"){
NumSample=ncol(do.call(rbind, IsoSIMout[[i]]$generateData))

NumIso=rep(0,Num)
for (i in 1:Num)NumIso[i]=nrow(do.call(rbind, IsoSIMout[[i]]$generateData))

MinNumIso=min(NumIso)
AproxNumDE=length(unlist(IsoSIMout[[1]]$TrueDE))
	
IsoMergeTable=matrix(rep(0,60),nrow=10)
	for(i in 1:Num)IsoMergeTable=IsoMergeTable+IsoSIMout[[i]][[1]]
	IsoMergeTable=IsoMergeTable/Num
	IsoMergeTable=round(IsoMergeTable,2)
	          
	IsoMergeDVD=rep(0,2)
	  for(i in 1:Num)IsoMergeDVD=IsoMergeDVD+IsoSIMout[[i]][[3]]
		  IsoMergeDVD=round(IsoMergeDVD/Num,2) 
				          
	  IsoMergePhi=matrix(rep(0,18),nrow=2)
		  for(i in 1:Num)IsoMergePhi=IsoMergePhi+IsoSIMout[[i]][[4]]
			  IsoMergePhi=round(IsoMergePhi/Num,2)
## Write
TXTname=paste(paste("../IsoOutput/",paste("Iso","DVD",IsoMergeDVD[1], IsoMergeDVD[2],"Sample",NumSample,sep="_"),sep=""),".txt",sep="")
write.table(IsoMergeTable, file=TXTname)


####### Note everytime # DE genes and # total genes may different. (since NA issue)
  IsoMergeFD=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeFD.p=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeTP.p=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeFN.p=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeTN.p=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeFDR=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeTPR=matrix(rep(0,5*MinNumIso),ncol=5)
  IsoMergeFPR=matrix(rep(0,5*MinNumIso),ncol=5)

  for(i in 1:Num){
	# Make sure names in the same order
	# Get FD number for each number of genes found
	# columns are samples 
    TotalNum=nrow(do.call(rbind, IsoSIMout[[i]]$generateData))
	NumDE=length(unlist(IsoSIMout[[i]]$TrueDE))
	EBSeqNames=names(IsoSIMout[[i]]$EBSeqPP)
    tmpMatrix=cbind(IsoSIMout[[i]]$DESeqP[EBSeqNames],IsoSIMout[[i]]$edgeRP[EBSeqNames], exp(IsoSIMout[[i]]$BaySeqPP[EBSeqNames,2]),IsoSIMout[[i]]$BBSeqP[EBSeqNames],IsoSIMout[[i]]$EBSeqPP)
	# Bayseq and EBseq are PP. Others are p value 
    tmpFD=TopCts(tmpMatrix, c(0,0,1,0,1), unlist(IsoSIMout[[i]]$TrueDE)[unlist(IsoSIMout[[i]]$TrueDE)%in%EBSeqNames], MinNumIso)
    # Get percentage for FP, TP, TN, FN!
	tmpFD.p=tmpFD/TotalNum
	# TP = Find - FD
	tmpTP.p=(outer(c(1:MinNumIso),rep(1,5))-tmpFD)/TotalNum
	# FN = TrueDE - TP
	tmpFN.p=NumDE/TotalNum - tmpTP.p
	# TN = TrueEE - FD
	tmpTN.p=(TotalNum-NumDE)/TotalNum - tmpFD.p
	
	tmpFDR=tmpFD.p/(tmpFD.p+tmpTP.p)
	tmpFPR=tmpFD.p/(tmpFD.p+tmpTN.p)
	tmpTPR=tmpTP.p/(tmpFN.p+tmpTP.p)
	IsoMergeFDR=IsoMergeFDR+tmpFDR
	IsoMergeTPR=IsoMergeTPR+tmpTPR
	IsoMergeFPR=IsoMergeFPR+tmpFPR

    IsoMergeFD.p=IsoMergeFD.p+tmpFD.p
	IsoMergeTP.p=IsoMergeTP.p+tmpTP.p
	IsoMergeFN.p=IsoMergeFN.p+tmpFN.p
	IsoMergeTN.p=IsoMergeTN.p+tmpTN.p

	IsoMergeFD=IsoMergeFD+tmpFD
 }   
  IsoMergeFD=IsoMergeFD/Num
  IsoMergeFD.p=IsoMergeFD.p/Num
  IsoMergeTP.p=IsoMergeTP.p/Num
  IsoMergeFN.p=IsoMergeFN.p/Num
  IsoMergeTN.p=IsoMergeTN.p/Num
  IsoMergeFDR=IsoMergeFDR/Num
  IsoMergeTPR=IsoMergeTPR/Num
  IsoMergeFPR=IsoMergeFPR/Num

PlotTopName=paste(paste(Path,paste("Top","Iso","DVD",IsoMergeDVD[1], IsoMergeDVD[2],"Sample",NumSample, sep="_"),sep=""),".pdf",sep="")

TrueDELength=length(unlist(IsoSIMout[[i]]$TrueDE)[unlist(IsoSIMout[[i]]$TrueDE)%in%EBSeqNames])
pdf(PlotTopName)
  PlotTopCts(TrueDELength,IsoMergeFD[1:TrueDELength,],c("DESeq","edgeR","BaySeq","BBSeq","EBSeq"))
dev.off()


PlotFDName=paste(paste(Path,paste("FDTP","Iso","DVD",IsoMergeDVD[1], IsoMergeDVD[2],"Sample",NumSample,sep="_"),sep=""),".pdf",sep="")
pdf(PlotFDName)
  PlotFDTP(MinNumIso,IsoMergeFDR, IsoMergeTPR, c("DESeq","edgeR","BaySeq","BBSeq","EBSeq"))
dev.off()

PlotFPName=paste(paste(Path,paste("FPRTP","Iso","DVD",IsoMergeDVD[1], IsoMergeDVD[2],"Sample",NumSample,sep="_"),sep=""),".pdf",sep="")
pdf(PlotFPName)
  PlotFPTP(MinNumIso,IsoMergeFPR, IsoMergeTPR, c("DESeq","edgeR","BaySeq","BBSeq","EBSeq"))
  dev.off()


out=list(IsoMergeTable=IsoMergeTable, IsoMergeDVD=IsoMergeDVD, IsoMergePhi=IsoMergePhi, IsoMergeFD=IsoMergeFD)


}

