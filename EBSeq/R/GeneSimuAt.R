GeneSimuAt<-function(DVDconstant=NULL, DVDqt1=NULL, DVDqt2=NULL, Conditions, NumofSample, NumofGene=NULL, DEGeneProp, Phiconstant=NULL, Phi.qt1=NULL, Phi.qt2=NULL, Meanconstant=NULL,NormFactor=NULL, OnlyData=T)
{
# 2012 feb 1 
# paired level simulation

data(GeneEBresultGouldBart2)
if(is.null(NormFactor)) NormFactor=rep(1,NumofSample)

#MeansC1=rowMeans(GeneV.norm1.NZ.b2[,1:4])
#MeansC2=rowMeans(GeneV.norm1.NZ.b2[,5:8])
MeansC1=GeneEBresultGouldBart2$C1Mean[[1]]
MeansC2=GeneEBresultGouldBart2$C2Mean[[1]]

MeanDVD=MeansC1/MeansC2

if(is.null(DVDconstant))DVDLibrary=MeanDVD[MeanDVD<quantile(MeanDVD[MeanDVD!=Inf],DVDqt2) & MeanDVD>quantile(MeanDVD[MeanDVD!=Inf],DVDqt1)]


# If DVD constant, use constant when generate
# If not, use DVDLibrary

MeanInputraw=GeneEBresultGouldBart2$MeanList[[1]]
#MeanInputraw=rowMeans(GeneV.norm1.NZ.b2)
#Var1=apply(GeneV.norm1.NZ.b2[,1:4],1,var)
#Var2=apply(GeneV.norm1.NZ.b2[,5:8],1,var)
#VarInput=(Var1 + Var2)/2
#If NumofGene.raw=NULL, empirical # of Gene
#If !=NULL , Input a 9-vector
NumofGene.raw=length(MeanInputraw)

# here phi denotes r -- which is 1/phi' in which sigma^2=mu(1+mu phi')
# In negative binomial 
# size is 1/phi'
# rnbinom(100,size=100,mu=10) 
# var(qq)
#[1] 10.93687 
# qq=rnbinom(100,size=10,mu=10)
# var(qq)
#[1] 24.01404

#PhiInput.raw=(MeanInputraw^2) / (VarInput - MeanInputraw)
PhiInput.raw=GeneEBresultGouldBart2$RList[[1]]
if (length(Phiconstant)==0){
	PhiLibrary=PhiInput.raw[1/(PhiInput.raw)<quantile(1/(PhiInput.raw),Phi.qt2) & 1/(PhiInput.raw)>quantile(1/(PhiInput.raw),Phi.qt1)]
    PhiInputNames=sample(names(PhiLibrary),NumofGene.raw,replace=T)
	PhiInput=PhiInput.raw[PhiInputNames]


}

if (length(Phiconstant)!=0)PhiInput=rep(Phiconstant,length(MeanInputraw))
if(length(Meanconstant)==0)MeanInput=GeneEBresultGouldBart2$MeanList[[1]][PhiInputNames]
if(length(Meanconstant)!=0)MeanInput=rep(Meanconstant,length(GeneEBresultGouldBart2$MeanList[[1]]))

# Wanna DENumbers be proportion to 2 
DEGeneNumbers=round(NumofGene.raw*DEGeneProp/2)*2
GeneNames=paste("G",c(1:NumofGene.raw),sep="_")
names(PhiInput)=GeneNames
names(MeanInput)=GeneNames

#########
# data
#########
EEList=sapply(1:NumofGene.raw, function(j) sapply(1:NumofSample, function(i)rnbinom(1,mu=NormFactor[i]*MeanInput[j], size=PhiInput[j])))




    generateDataraw=t(EEList)
	if(length(DVDconstant)==0){
		DVDSample=sample(DVDLibrary,DEGeneNumbers,replace=T)
		for(j in 1:NumofGene.raw){
	    	 if (j<=(DEGeneNumbers/2)) generateDataraw[j,((NumofSample/2)+1):NumofSample]=sapply(((NumofSample/2) +1):NumofSample, function(i)rnbinom(1, size=PhiInput[j], mu=DVDSample[j]*MeanInput[j]*NormFactor[i]),simplify=T)
	     	if (j>=((DEGeneNumbers/2)+1) & j <=DEGeneNumbers) generateDataraw[j,1:(NumofSample/2)]=sapply(1:(NumofSample/2),function(i)rnbinom(1, size=MeanInput[j], mu= DVDSample[j]*MeanInput[j]*NormFactor[i]),simplify=T)
}
	 }
	if(length(DVDconstant)!=0){
        for(j in 1:NumofGene.raw){
             if (j<=(DEGeneNumbers/2)) generateDataraw[j,((NumofSample/2)+1):NumofSample]=sapply((NumofSample/2+1):NumofSample, function(i)rnbinom(1, size=MeanInput[j],mu=DVDconstant*MeanInput[j]*NormFactor[i]),simplify=T)
             if (j>=((DEGeneNumbers/2)+1) & j <=DEGeneNumbers) generateDataraw[j,1:(NumofSample/2)]=sapply(1:(NumofSample/2),function(i)rnbinom(1, size=MeanInput[j],mu=DVDconstant*MeanInput[j]*NormFactor[i]),simplify=T)
		}
	}
rownames(generateDataraw)=GeneNames
MeanVector=rowMeans(generateDataraw)
VarVector=apply(generateDataraw,1,var)
MOV.post=MeanVector/VarVector



### Remove MOV=NA
generateData=generateDataraw
generateData=generateData[!is.na(MOV.post)& MeanVector>2 & MeanVector<10000 ,] 
print(paste("NA MOV's",sum(is.na(MOV.post)),sum( MeanVector<2), sum(MeanVector>10000)))
## DE
NumDENow=sum(rownames(generateData)%in%rownames(generateDataraw)[1:DEGeneNumbers])

if(length(NumofGene)!=0)
    generateData=generateData[c(sample(1:NumDENow,round(NumofGene*DEGeneProp),replace=F),round( (dim(generateData)[1]+1-NumofGene*(1-DEGeneProp)):dim(generateData)[1])),]


UseName=rownames(generateData)

TrueDE=UseName[UseName%in%rownames(generateDataraw)[1:DEGeneNumbers]]
phiuse=PhiInput[rownames(generateData)]
meanuse=MeanInput[rownames(generateData)]

#ArtiNames=rownames(generateData)[(DEGeneNumbers+1):(2*DEGeneNumbers)]
#Noise=sample(c(1,ncol(generateData)),DEGeneNumbers,replace=T)
TrueDELength=length(TrueDE)
AtLoc=sample(c(1:length(Conditions)), TrueDELength, replace=T)
AtFold=sample(c(4,6,8,10),TrueDELength, replace=T)

AtNames_Level=vector("list",4)
names(AtNames_Level)=c(4,6,8,10)
for(i in 1:TrueDELength){
generateData[(TrueDELength+i),AtLoc[i]]=generateData[(TrueDELength+i),AtLoc[i]]*AtFold[i]
AtNames_Level[[as.character(AtFold[i])]]=c(AtNames_Level[[as.character(AtFold[i])]],rownames(generateData)[TrueDELength+i])
}


if(OnlyData==T){
	OutName=paste("Gene",c(1:nrow(generateData)),sep="_")
	names(OutName)=rownames(generateData)
    OutData=generateData
    rownames(OutData)=as.vector(OutName)
	OutAt=as.vector(OutName[AtNames_Level])
	OutTrueDE=as.vector(OutName[TrueDE])
    output=list(data=OutData, TrueDE=OutTrueDE,Outliers=OutAt)
	return(output)
	}
## DESeq

cds=newCountDataSet(round(generateData),Conditions)
cds=estimateSizeFactors(cds)
Sizes=sizeFactors(cds)
if(dim(generateData)[2]>4)cds=estimateVarianceFunctions(cds)
else  cds=estimateVarianceFunctions(cds, method="blind")

res=nbinomTest(cds, "1", "2")
ResAdj=res$padj
names(ResAdj)=res$id
SmallPValueName=names(ResAdj)[which(ResAdj<=.05)]
print(paste("DESEq found",length(SmallPValueName)))
print(paste("In True DE",sum(SmallPValueName%in%TrueDE)))

print("DESeq Size factors")
print(Sizes)

## DESeq each group
## Ours
NewData=generateData


#source("/z/Comp/kendziorskigroup/ningleng/RNASEQ/CODE/FinalV/NBBetaBiasUniqueP_PoolVar_SpeedUp_MDFPoi_NoNormVar.R")
#source("/z/Comp/kendziorskigroup/ningleng/RNASEQ/CODE/FinalV/NBBetaBiasUniqueP_PoolVar_SpeedUp_MDFPoi_NoNormPoolR.R")

EBresult=EBTest(NewData,rep(1,dim(NewData)[1]), rep(1,dim(NewData)[1]), rep(1,dim(NewData)[1]),Conditions,sizeFactors=Sizes,5)

#EBres2=NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi_NoNormPoolR(NewData,rep(1,dim(NewData)[1]), rep(1,dim(NewData)[1]), rep(1,dim(NewData)[1]),Conditions,sizeFactors=Sizes,5)


zlist.unlist=EBresult[[5]]
fdr=max(.5,crit_fun(1-zlist.unlist,.05))
EBDE=names(zlist.unlist)[which(zlist.unlist>fdr)]
EBDE.Poi=names(EBresult[[6]])[which(EBresult[[6]]>fdr)]
zlist.unlist.whole=c(EBresult[[5]],EBresult[[6]])
print(paste("Soft EB Poi",length(EBDE.Poi)))
EBDE=c(EBDE, EBDE.Poi)
print(paste("Soft EB found",length(EBDE)))
print(paste("In True DE",sum(EBDE%in%TrueDE)))

EBDE95=names(zlist.unlist)[which(zlist.unlist>.95)]
EBDE95.Poi=names(EBresult[[6]])[which(EBresult[[6]]>.95)]
print(paste("Hard Poi found",length(EBDE95.Poi)))
EBDE95=c(EBDE95, EBDE95.Poi)
print(paste("Hard EB found" ,length(EBDE95)))
print(paste("In True DE",sum(EBDE95%in%TrueDE)))

### edgeR
library(edgeR,lib.loc="~/RCODE")
edgeRList.b2=DGEList(NewData,group=Conditions)
if(length(Phiconstant)==1){
	edgeRList.b2=estimateCommonDisp(edgeRList.b2)
	edgeRRes.b2=exactTest(edgeRList.b2)
}
if(length(Phiconstant)==0){
	edgeRList.b2=estimateCommonDisp(edgeRList.b2)	
	edgeRList.b2=estimateTagwiseDisp(edgeRList.b2)
	edgeRRes.b2=exactTest(edgeRList.b2, common.disp = FALSE)
}
edgeRPvalue.b2.raw=edgeRRes.b2[[1]][[3]]
edgeRPvalue.b2=p.adjust(edgeRPvalue.b2.raw, method="BH")
names(edgeRPvalue.b2)=rownames(NewData)
edgeRSmallpvalue=names(which(edgeRPvalue.b2<.05))
print(paste("edgeR found",length(edgeRSmallpvalue)))
print(paste("In True DE",sum(edgeRSmallpvalue%in%TrueDE)))

### Bayseq
library(baySeq, lib.loc="~/RCODE")
library(snow, lib.loc="~/RCODE")
cl <- makeCluster(4, "SOCK")
groups <- list(NDE = rep(1,NumofSample), DE = rep(c(1,2),each=NumofSample/2))
CD <- new("countData", data = NewData, replicates = Conditions, libsizes = as.integer(colSums(NewData)), groups = groups)
CDP.NBML <- getPriors.NB(CD, samplesize = dim(NewData)[1], estimation = "QL", cl = cl)
CDPost.NBML <- getLikelihoods.NB(CDP.NBML, pET = "BIC", cl = cl)
bayseqPost=CDPost.NBML@posteriors
rownames(bayseqPost)=rownames(NewData)
bayseqDE=rownames(NewData)[bayseqPost[,2]>log(.95)]
print(paste("bayseq found",length(bayseqDE)))
print(paste("In True DE",sum(bayseqDE%in%TrueDE)))


### BBSeq
library("BBSeq",lib.loc="~/RCODE")
CondM=cbind(rep(1,NumofSample),rep(c(0,1),each=NumofSample/2))
output=free.estimate(NewData,CondM)
beta.free = output$betahat.free
p.free = output$p.free
psi.free = output$psi.free
names(p.free)=rownames(NewData)
p.free.adj=p.adjust(p.free,method="BH")
# Top p free?
#out.model=constrained.estimate(NewData,CondM, gn=3, beta.free ,psi.free)
#p.constrained = out.model$p.model
BBDE=names(p.free.adj)[which(p.free.adj<.05)]
print(paste("BBSeq found",length(BBDE)))
print(paste("In True DE",sum(BBDE%in%TrueDE)))


#########################
# Generate table
Table=matrix(rep(0,12),ncol=2)
colnames(Table)=c("Power","FDR")
rownames(Table)=c("DESeq","edgeR","BaySeq","BBSeq","EBSeq_ModifiedSoft","EBSeq_Hard")

	Length=length(TrueDE)
	Table[1,1]=sum(SmallPValueName%in%TrueDE)/Length
	Table[2,1]=sum(edgeRSmallpvalue%in%TrueDE)/Length
	Table[3,1]=sum(bayseqDE%in%TrueDE)/Length
	Table[4,1]=sum(BBDE%in%TrueDE)/Length
	Table[5,1]=sum(EBDE%in%TrueDE)/Length
	Table[6,1]=sum(EBDE95%in%TrueDE)/Length
	Table[1,2]=sum(!SmallPValueName%in%TrueDE)/length(SmallPValueName)
	Table[2,2]=sum(!edgeRSmallpvalue%in%TrueDE)/length(edgeRSmallpvalue)
	Table[3,2]=sum(!bayseqDE%in%TrueDE)/length(bayseqDE)
	Table[4,2]=sum(!BBDE%in%TrueDE)/length(BBDE)
	Table[5,2]=sum(!EBDE%in%TrueDE)/length(EBDE)
	Table[6,2]=sum(!EBDE95%in%TrueDE)/length(EBDE95)
	Table=round(Table,2)

ValueTable=matrix(rep(0,12),ncol=2)
colnames(ValueTable)=c("Power","FDR")
rownames(ValueTable)=c("DESeq","edgeR","BaySeq","BBSeq","EBSeq_ModifiedSoft","EBSeq_Hard")
 	ValueTable[1,1]=sum(SmallPValueName%in%TrueDE)
	ValueTable[2,1]=sum(edgeRSmallpvalue%in%TrueDE)
	ValueTable[3,1]=sum(bayseqDE%in%TrueDE)
	ValueTable[4,1]=sum(BBDE%in%TrueDE)
	ValueTable[5,1]=sum(EBDE%in%TrueDE)
	ValueTable[6,1]=sum(EBDE95%in%TrueDE)
	ValueTable[1,2]=sum(!SmallPValueName%in%TrueDE)
	ValueTable[2,2]=sum(!edgeRSmallpvalue%in%TrueDE)
	ValueTable[3,2]=sum(!bayseqDE%in%TrueDE)
	ValueTable[4,2]=sum(!BBDE%in%TrueDE)
	ValueTable[5,2]=sum(!EBDE%in%TrueDE)
	ValueTable[6,2]=sum(!EBDE95%in%TrueDE)


AtFoundTable=matrix(rep(0,24),ncol=4)
colnames(AtFoundTable)=paste("Level",c(1:4),sep="_")
rownames(Table)=c("DESeq","edgeR","BaySeq","BBSeq","EBSeq_ModifiedSoft","EBSeq_Hard")
for(i in 1:4){
 	AtFoundTable[1,i]=sum(SmallPValueName%in%AtNames_Level[[i]])
	AtFoundTable[2,i]=sum(edgeRSmallpvalue%in%AtNames_Level[[i]])
	AtFoundTable[3,i]=sum(bayseqDE%in%AtNames_Level[[i]])
	AtFoundTable[4,i]=sum(BBDE%in%AtNames_Level[[i]])
	AtFoundTable[5,i]=sum(EBDE%in%AtNames_Level[[i]])
	AtFoundTable[6,i]=sum(EBDE95%in%AtNames_Level[[i]])	
	}

	
if(length(DVDconstant)==0)DVD=c(quantile(MeanDVD[MeanDVD!=Inf],DVDqt1), quantile(MeanDVD[MeanDVD!=Inf],DVDqt2))
if(length(DVDconstant)!=0) DVD=DVDconstant
if(length(Phiconstant)==0)Phi=c(quantile(PhiInput.raw,Phi.qt1), quantile(PhiInput.raw,Phi.qt2))
if(length(Phiconstant)!=0) Phi=Phiconstant
OUT=list(Table=Table, ValueTable=ValueTable, DVD=DVD, Phi=Phi, generateData=NewData, TrueDE=TrueDE,phi.vector=phiuse,mean.vector=meanuse,NormFactor=NormFactor, DESeqP=ResAdj, edgeRP=edgeRPvalue.b2, EBSeqPP=zlist.unlist.whole, BaySeqPP=bayseqPost,BBSeqP=p.free.adj,EBoutput=EBresult,  AtFoundTable= AtFoundTable,Outliers=AtNames_Level)



}


