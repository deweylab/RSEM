EBTest <-
function(Data,NgVector=NULL,Vect5End=NULL,Vect3End=NULL,Conditions, sizeFactors, maxround, tau=NULL,CI=NULL,CIthre=NULL, Pool=F, NumBin=1000,ApproxVal=10^-10)
{
	Dataraw=Data
	AllZeroNames=which(rowMeans(Data)==0)
	NotAllZeroNames=which(rowMeans(Data)>0)
	if(length(AllZeroNames)>0) print("Remove transcripts with all zero")
	Data=Data[NotAllZeroNames,]
	if(!is.null(NgVector))NgVector=NgVector[NotAllZeroNames]
	if(!length(sizeFactors)==ncol(Data))sizeFactors=sizeFactors[NotAllZeroNames,]

	if(is.null(NgVector))NgVector=rep(1,nrow(Data))

	#Rename Them
	IsoNamesIn=rownames(Data)
	Names=paste("I",c(1:dim(Data)[1]),sep="")
	names(IsoNamesIn)=Names
	rownames(Data)=paste("I",c(1:dim(Data)[1]),sep="")
	names(NgVector)=paste("I",c(1:dim(Data)[1]),sep="")
	

	if(!length(sizeFactors)==ncol(Data)){
		rownames(sizeFactors)=rownames(Data)
		colnames(sizeFactors)=Conditions
	}
	
	NumOfNg=nlevels(as.factor(NgVector))
	NameList=sapply(1:NumOfNg,function(i)Names[NgVector==i],simplify=F)
	names(NameList)=paste("Ng",c(1:NumOfNg),sep="")
	NotNone=NULL
	for (i in 1:NumOfNg) {
		if (length(NameList[[i]])!=0) 
			NotNone=c(NotNone,names(NameList)[i])
		}
	NameList=NameList[NotNone]
		
	NoneZeroLength=length(NameList)
	DataList=vector("list",NoneZeroLength)
	DataList=sapply(1:NoneZeroLength , function(i) Data[NameList[[i]],],simplify=F)
	names(DataList)=names(NameList)
    
	NumEachGroup=sapply(1:NoneZeroLength , function(i)dim(DataList)[i])
	# Unlist 
	DataList.unlist=do.call(rbind, DataList)

	# Divide by SampleSize factor
	
	if(length(sizeFactors)==ncol(Data))
	DataList.unlist.dvd=t(t( DataList.unlist)/sizeFactors)
	
	if(length(sizeFactors)!=ncol(Data))
	DataList.unlist.dvd=DataList.unlist/sizeFactors
	
	# Get FC and VarPool for pooling - Only works on 2 conditions
	if(ncol(Data)==2){
	DataforPoolSP.dvd1=matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[1]],nrow=dim(DataList.unlist)[1])	
	DataforPoolSP.dvd2=matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[2]],nrow=dim(DataList.unlist)[1])
	MeanforPoolSP.dvd1=rowMeans(DataforPoolSP.dvd1)
	MeanforPoolSP.dvd2=rowMeans(DataforPoolSP.dvd2)
	FCforPool=MeanforPoolSP.dvd1/MeanforPoolSP.dvd2
	names(FCforPool)=rownames(Data)
	FC_Use=which(FCforPool>=quantile(FCforPool[!is.na(FCforPool)],.25) & 
								  FCforPool<=quantile(FCforPool[!is.na(FCforPool)],.75))
	
	Var_FC_Use=apply( DataList.unlist.dvd[FC_Use,],1,var )
	Mean_FC_Use=(MeanforPoolSP.dvd1[FC_Use]+MeanforPoolSP.dvd2[FC_Use])/2
	MeanforPool=(MeanforPoolSP.dvd1+MeanforPoolSP.dvd2)/2
	FC_Use2=which(Var_FC_Use>=Mean_FC_Use)
	Var_FC_Use2=Var_FC_Use[FC_Use2]
	Mean_FC_Use2=Mean_FC_Use[FC_Use2]
	Phi=mean((Var_FC_Use2-Mean_FC_Use2)/Mean_FC_Use2^2)
	VarEst=	MeanforPool*(1+MeanforPool*Phi)
	print(Phi)
	}

	#DataListSP Here also unlist.. Only two lists
	DataListSP=vector("list",nlevels(Conditions))
	DataListSP.dvd=vector("list",nlevels(Conditions))
	SizeFSP=DataListSP
	MeanSP=DataListSP
	VarSP=DataListSP
	GetPSP=DataListSP
	RSP=DataListSP
	CISP=DataListSP
	tauSP=DataListSP
	NumSampleEachCon=rep(NULL,nlevels(Conditions))

	for (lv in 1:nlevels(Conditions)){
		DataListSP[[lv]]= matrix(DataList.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist)[1])
		rownames(DataListSP[[lv]])=rownames(DataList.unlist)
		DataListSP.dvd[[lv]]= matrix(DataList.unlist.dvd[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
		NumSampleEachCon[lv]=ncol(DataListSP[[lv]])

	if(ncol(DataListSP[[lv]])==1 & !is.null(CI)){
		CISP[[lv]]=matrix(CI[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
		tauSP[[lv]]=matrix(tau[,Conditions==levels(Conditions)[lv]],nrow=dim(DataList.unlist.dvd)[1])
	}
	# no matter sizeFactors is a vector or a matrix. Matrix should be columns are the normalization factors
	# may input one for each 
	if(length(sizeFactors)==ncol(Data))SizeFSP[[lv]]=sizeFactors[Conditions==levels(Conditions)[lv]]
	if(length(sizeFactors)!=ncol(Data))SizeFSP[[lv]]=sizeFactors[,Conditions==levels(Conditions)[lv]]
	
	
	MeanSP[[lv]]=rowMeans(DataListSP.dvd[[lv]])
	
	if(length(sizeFactors)==ncol(Data))PrePareVar=sapply(1:ncol( DataListSP[[lv]]),function(i)( DataListSP[[lv]][,i]- SizeFSP[[lv]][i]*MeanSP[[lv]])^2 /SizeFSP[[lv]][i])
	if(length(sizeFactors)!=ncol(Data))PrePareVar=sapply(1:ncol( DataListSP[[lv]]),function(i)( DataListSP[[lv]][,i]- SizeFSP[[lv]][,i]*MeanSP[[lv]])^2 /SizeFSP[[lv]][,i])

	if(ncol(DataListSP[[lv]])==1 & !is.null(CI))
		VarSP[[lv]]=as.vector(((DataListSP[[lv]]/tauSP[[lv]]) * CISP[[lv]]/(CIthre*2))^2)
	if(ncol(DataListSP[[lv]])!=1){
		VarSP[[lv]]=rowSums(PrePareVar)/ncol( DataListSP[[lv]])
		names(MeanSP[[lv]])=rownames(DataList.unlist)
		names(VarSP[[lv]])=rownames(DataList.unlist)
		GetPSP[[lv]]=MeanSP[[lv]]/VarSP[[lv]]
		RSP[[lv]]=MeanSP[[lv]]*GetPSP[[lv]]/(1-GetPSP[[lv]])
	}
}
	
	
	MeanList=rowMeans(DataList.unlist.dvd)
	VarList=apply(DataList.unlist.dvd, 1, var)
	if(ncol(Data)==2)PoolVar=VarEst
	if(!ncol(Data)==2){
		CondWithRep=which(NumSampleEachCon>1)
		VarCondWithRep=do.call(cbind,VarSP[CondWithRep])
		PoolVar=rowMeans(VarCondWithRep)
	}
	GetP=MeanList/PoolVar
	
    EmpiricalRList=MeanList*GetP/(1-GetP) 
	EmpiricalRList[EmpiricalRList==Inf]	=max(EmpiricalRList[EmpiricalRList!=Inf])
	
	if(ncol(Data)!=2){
	Varcbind=do.call(cbind,VarSP)
	VarrowMin=apply(Varcbind,1,min)
	}

	if(ncol(Data)==2){
		Varcbind=VarEst
		VarrowMin=VarEst
	}
	# 
	# 
	GoodData=names(MeanList)[EmpiricalRList>0 &  VarrowMin!=0 & EmpiricalRList!=Inf & !is.na(VarrowMin) & !is.na(EmpiricalRList)]
	NotIn=names(MeanList)[EmpiricalRList<=0 | VarrowMin==0 | EmpiricalRList==Inf |  is.na(VarrowMin) | is.na(EmpiricalRList)]
	#print(paste("ZeroVar",sum(VarrowMin==0), "InfR", length(which(EmpiricalRList==Inf)), "Poi", length(which(EmpiricalRList<0)), ""))
	EmpiricalRList.NotIn=EmpiricalRList[NotIn]
	EmpiricalRList.Good=EmpiricalRList[GoodData]
	EmpiricalRList.Good[EmpiricalRList.Good<1]=1+EmpiricalRList.Good[EmpiricalRList.Good<1]
	if(length(sizeFactors)==ncol(Data))
	EmpiricalRList.Good.mat= outer(EmpiricalRList.Good, sizeFactors)	
	if(!length(sizeFactors)==ncol(Data))
	EmpiricalRList.Good.mat=EmpiricalRList.Good* sizeFactors[GoodData,]


	# Only Use Data has Good q's
	DataList.In=sapply(1:NoneZeroLength, function(i)DataList[[i]][GoodData[GoodData%in%rownames(DataList[[i]])],],simplify=F)
	DataList.NotIn=sapply(1:NoneZeroLength, function(i)DataList[[i]][NotIn[NotIn%in%rownames(DataList[[i]])],],simplify=F)
	DataListIn.unlist=do.call(rbind, DataList.In)
	DataListNotIn.unlist=do.call(rbind, DataList.NotIn)
	
	DataListSPIn=vector("list",nlevels(Conditions))
	DataListSPNotIn=vector("list",nlevels(Conditions))
	EmpiricalRList.Good.mat.SP=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)){
		DataListSPIn[[lv]]= matrix(DataListIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListIn.unlist)[1])
	if(length(NotIn)>0){	DataListSPNotIn[[lv]]= matrix(DataListNotIn.unlist[,Conditions==levels(Conditions)[lv]],nrow=dim(DataListNotIn.unlist)[1])
	rownames(DataListSPNotIn[[lv]])=rownames(DataListNotIn.unlist)
	}
	rownames(DataListSPIn[[lv]])=rownames(DataListIn.unlist)
	EmpiricalRList.Good.mat.SP[[lv]]=matrix(EmpiricalRList.Good.mat[,Conditions==levels(Conditions)[lv]],nrow=dim(EmpiricalRList.Good.mat)[1])
}	

	NumOfEachGroupIn=sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.In[[i]])[1]))
	NumOfEachGroupNotIn=sapply(1:NoneZeroLength, function(i)max(0,dim(DataList.NotIn[[i]])[1]))

	#Initialize SigIn & ...
	AlphaIn=0.5
	BetaIn=rep(0.5,NoneZeroLength)
	PIn=0.5

	####use while to make an infinity round?
	UpdateAlpha=NULL
	UpdateBeta=NULL
	UpdateP=NULL
	UpdatePFromZ=NULL
    Timeperround=NULL 
	for (times in 1:maxround){
    	temptime1=proc.time()
		UpdateOutput=suppressWarnings(LogN(DataListIn.unlist,DataListSPIn, EmpiricalRList.Good.mat ,EmpiricalRList.Good.mat.SP,  NumOfEachGroupIn, AlphaIn, BetaIn, PIn, NoneZeroLength))
    	print(paste("iteration", times, "done",sep=" "))
		AlphaIn=UpdateOutput$AlphaNew
    	BetaIn=UpdateOutput$BetaNew
    	PIn=UpdateOutput$PNew
		PFromZ=UpdateOutput$PFromZ
    	F0Out=UpdateOutput$F0Out
		F1Out=UpdateOutput$F1Out
		UpdateAlpha=rbind(UpdateAlpha,AlphaIn)
   		UpdateBeta=rbind(UpdateBeta,BetaIn)
    	UpdateP=rbind(UpdateP,PIn)
		UpdatePFromZ=rbind(UpdatePFromZ,PFromZ)
		temptime2=proc.time()
		Timeperround=c(Timeperround,temptime2[3]-temptime1[3])
		print(paste("time" ,Timeperround[times],sep=" "))
		Z.output=UpdateOutput$ZNew.list[!is.na(UpdateOutput$ZNew.list)]
   		Z.NA.Names=UpdateOutput$zNaNName
		}
		#Remove this } after testing!!
		 
#    	if (times!=1){  
#        	if((UpdateAlpha[times]-UpdateAlpha[times-1])^2+UpdateBeta[times]-UpdateBeta[times-1])^2+UpdateR[times]-UpdateR[times-1])^2+UpdateP[times]-UpdateP[times-1])^2<=10^(-6)){ 
#           		Result=list(Sig=SigIn, Miu=MiuIn, Tau=TauIn)
#           		break
#        }
#    }
#}

##########Change Names############
## Only z are for Good Ones
## Others are for ALL Data
GoodData=GoodData[!GoodData%in%Z.NA.Names]
IsoNamesIn.Good=IsoNamesIn[GoodData]
RealName.Z.output=Z.output
RealName.F0=F0Out
RealName.F1=F1Out
names(RealName.Z.output)=IsoNamesIn.Good
names(RealName.F0)=IsoNamesIn.Good
names(RealName.F1)=IsoNamesIn.Good


RealName.EmpiricalRList=sapply(1:NoneZeroLength,function(i)EmpiricalRList[names(EmpiricalRList)%in%NameList[[i]]], simplify=F)
RealName.MeanList=sapply(1:NoneZeroLength,function(i)MeanList[names(MeanList)%in%NameList[[i]]], simplify=F)
RealName.C1MeanList=sapply(1:NoneZeroLength,function(i)MeanSP[[1]][names(MeanSP[[1]])%in%NameList[[i]]], simplify=F)
RealName.C2MeanList=sapply(1:NoneZeroLength,function(i)MeanSP[[2]][names(MeanSP[[2]])%in%NameList[[i]]], simplify=F)
RealName.C1VarList=sapply(1:NoneZeroLength,function(i)VarSP[[1]][names(VarSP[[1]])%in%NameList[[i]]], simplify=F)
RealName.C2VarList=sapply(1:NoneZeroLength,function(i)VarSP[[2]][names(VarSP[[2]])%in%NameList[[i]]], simplify=F)
RealName.DataList=sapply(1:NoneZeroLength,function(i)DataList[[i]][rownames(DataList[[i]])%in%NameList[[i]],], simplify=F)



RealName.VarList=sapply(1:NoneZeroLength,function(i)VarList[names(VarList)%in%NameList[[i]]], simplify=F)
RealName.PoolVarList=sapply(1:NoneZeroLength,function(i)PoolVar[names(PoolVar)%in%NameList[[i]]], simplify=F)


RealName.QList1=sapply(1:NoneZeroLength,function(i)GetPSP[[1]][names(GetPSP[[1]])%in%NameList[[i]]], simplify=F)
RealName.QList2=sapply(1:NoneZeroLength,function(i)GetPSP[[2]][names(GetPSP[[2]])%in%NameList[[i]]], simplify=F)


for (i in 1:NoneZeroLength){
tmp=NameList[[i]]
names=IsoNamesIn[tmp]

RealName.MeanList[[i]]=RealName.MeanList[[i]][NameList[[i]]]
RealName.VarList[[i]]=RealName.VarList[[i]][NameList[[i]]]
RealName.QList1[[i]]=RealName.QList1[[i]][NameList[[i]]]
RealName.QList2[[i]]=RealName.QList2[[i]][NameList[[i]]]
RealName.EmpiricalRList[[i]]=RealName.EmpiricalRList[[i]][NameList[[i]]]
RealName.C1MeanList[[i]]=RealName.C1MeanList[[i]][NameList[[i]]]
RealName.C2MeanList[[i]]=RealName.C2MeanList[[i]][NameList[[i]]]
RealName.PoolVarList[[i]]=RealName.PoolVarList[[i]][NameList[[i]]]
RealName.C1VarList[[i]]=RealName.C1VarList[[i]][NameList[[i]]]
RealName.C2VarList[[i]]=RealName.C2VarList[[i]][NameList[[i]]]
RealName.DataList[[i]]=RealName.DataList[[i]][NameList[[i]],]

names(RealName.MeanList[[i]])=names
names(RealName.VarList[[i]])=names
if(ncol(DataListSP[[1]])!=1){
	names(RealName.QList1[[i]])=names
	names(RealName.C1VarList[[i]])=names
}
if(ncol(DataListSP[[2]])!=1){
	names(RealName.QList2[[i]])=names
	names(RealName.C2VarList[[i]])=names
}

names(RealName.EmpiricalRList[[i]])=names
names(RealName.C1MeanList[[i]])=names
names(RealName.C2MeanList[[i]])=names
names(RealName.PoolVarList[[i]])=names
rownames(RealName.DataList[[i]])=names


}


#########posterior part for other data set here later############
AllNA=unique(c(Z.NA.Names,NotIn))
z.list.NotIn=NULL
AllF0=c(RealName.F0)
AllF1=c(RealName.F1)
AllZ=RealName.Z.output

if (length(AllNA)>0){
	Ng.NA=NgVector[AllNA]
	AllNA.Ngorder=AllNA[order(Ng.NA)]
	NumOfEachGroupNA=rep(0,NoneZeroLength)
	NumOfEachGroupNA.tmp=tapply(Ng.NA,Ng.NA,length)
	names(NumOfEachGroupNA)=c(1:NoneZeroLength)
	NumOfEachGroupNA[names(NumOfEachGroupNA.tmp)]=NumOfEachGroupNA.tmp
	PNotIn=rep(1-ApproxVal,length(AllNA.Ngorder))
	MeanList.NotIn=MeanList[AllNA.Ngorder]
	R.NotIn.raw=MeanList.NotIn*PNotIn/(1-PNotIn) 
	if(length(sizeFactors)==ncol(Data))
	R.NotIn=outer(R.NotIn.raw,sizeFactors)
	if(!length(sizeFactors)==ncol(Data))
	R.NotIn=R.NotIn.raw*sizeFactors[NotIn,]
	R.NotIn1=matrix(R.NotIn[,Conditions==levels(Conditions)[1]],nrow=nrow(R.NotIn))
	R.NotIn2=matrix(R.NotIn[,Conditions==levels(Conditions)[2]],nrow=nrow(R.NotIn))
    
	DataListNotIn.unlistWithZ=DataList.unlist[AllNA.Ngorder,]
	DataListSPNotInWithZ=vector("list",nlevels(Conditions))
	for (lv in 1:nlevels(Conditions)) 
		DataListSPNotInWithZ[[lv]] = matrix(DataListSP[[lv]][AllNA.Ngorder,],nrow=length(AllNA.Ngorder))
		F0=f0(DataListNotIn.unlistWithZ,  AlphaIn, BetaIn, R.NotIn, NumOfEachGroupNA, log=F)
    	F1=f1(DataListSPNotInWithZ[[1]], DataListSPNotInWithZ[[2]], AlphaIn, BetaIn, R.NotIn1,R.NotIn2, NumOfEachGroupNA, log=F)
	z.list.NotIn=PIn*F1/(PIn*F1+(1-PIn)*F0)
#	names(z.list.NotIn)=IsoNamesIn.Good=IsoNamesIn[which(Names%in%NotIn)]
	names(z.list.NotIn)=IsoNamesIn[AllNA.Ngorder]

	AllZ=c(RealName.Z.output,z.list.NotIn)
	AllZ=AllZ[IsoNamesIn]
	AllZ[is.na(AllZ)]=0
	F0.NotIn=F0
	F1.NotIn=F1
	names(F0.NotIn)=IsoNamesIn[names(F0)]
    names(F1.NotIn)=IsoNamesIn[names(F1)]
	AllF0=c(RealName.F0,F0.NotIn)
	AllF1=c(RealName.F1,F1.NotIn)
	AllF0=AllF0[IsoNamesIn]
	AllF1=AllF1[IsoNamesIn]
	AllF0[is.na(AllF0)]=0
	AllF1[is.na(AllF1)]=0
}
#############Result############################
Result=list(Alpha=UpdateAlpha,Beta=UpdateBeta,P=UpdateP,PFromZ=UpdatePFromZ, Z=RealName.Z.output,PoissonZ=z.list.NotIn, RList=RealName.EmpiricalRList, MeanList=RealName.MeanList, VarList=RealName.VarList, QList1=RealName.QList1, QList2=RealName.QList2, C1Mean=RealName.C1MeanList, C2Mean=RealName.C2MeanList,C1EstVar=RealName.C1VarList, C2EstVar=RealName.C2VarList, PoolVar=RealName.PoolVarList , DataList=RealName.DataList,PPDE=AllZ,f0=AllF0, f1=AllF1,
			AllZeroIndex=AllZeroNames)
}

