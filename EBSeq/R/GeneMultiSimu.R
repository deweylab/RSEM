GeneMultiSimu<-
function(DVDconstant=NULL, DVDqt1=NULL, DVDqt2=NULL, Conditions,AllParti, NumofSample, NumofGene=NULL, DEGeneProp, Phiconstant=NULL, Phi.qt1=NULL, Phi.qt2=NULL, Meanconstant=NULL,NormFactor=NULL, OnlyData=T)
{
# 2012 feb 1 paired simulation
if(is.null(NormFactor)) NormFactor=rep(1,NumofSample)
data(GeneEBresultGouldBart2)
MeansC1=GeneEBresultGouldBart2$C1Mean[[1]]
MeansC2=GeneEBresultGouldBart2$C2Mean[[1]]

MeanDVD=MeansC1/MeansC2

if(is.null(DVDconstant))DVDLibrary=MeanDVD[MeanDVD<quantile(MeanDVD[MeanDVD!=Inf],DVDqt2) & MeanDVD>quantile(MeanDVD[MeanDVD!=Inf],DVDqt1)]
if(!is.null(DVDconstant))DVDLibrary=DVDconstant

# If DVD constant, use constant when generate
# If not, use DVDLibrary

MeanInputraw=GeneEBresultGouldBart2$MeanList[[1]]

if(length(NumofGene)!=0)
NumofGene.raw=NumofGene*2

if(length(NumofGene)==0)
NumofGene.raw=length(MeanInputraw)


PhiInput.raw=GeneEBresultGouldBart2$RList[[1]]
if (length(Phiconstant)==0){
	PhiLibrary=PhiInput.raw[(1/PhiInput.raw)<quantile(1/PhiInput.raw,Phi.qt2) & 1/PhiInput.raw>quantile(1/PhiInput.raw,Phi.qt1)]
	PhiInputNames=sample(names(PhiLibrary),NumofGene.raw,replace=T)
	PhiInput=PhiInput.raw[PhiInputNames]
}

if (length(Phiconstant)!=0)PhiInput=rep(Phiconstant,length(MeanInputraw))
if(length(Meanconstant)==0)MeanInput=GeneEBresultGouldBart2$MeanList[[1]][PhiInputNames]
if(length(Meanconstant)!=0)MeanInput=rep(Meanconstant,length(GeneEBresultGouldBart2$MeanList[[1]]))

# length(DEGeneNumbers) should be num of patterns -1. the others EE
PatternGeneNumbers=round(NumofGene.raw*DEGeneProp/2)*2
names(PatternGeneNumbers)=rownames(AllParti)
EEWhich=which(rowSums(AllParti)==ncol(AllParti))
DEGeneNumbers=PatternGeneNumbers[-EEWhich]


OutGeneNumbers=round(NumofGene*DEGeneProp/2)*2
names(OutGeneNumbers)=rownames(AllParti)
OutDEGeneNumbers=OutGeneNumbers[-EEWhich]
OutEEGeneNumbers=OutGeneNumbers[EEWhich]
OutGenePatterns=c(unlist(sapply(1:length(OutDEGeneNumbers),
							  function(i)rep(names(OutDEGeneNumbers)[i],OutDEGeneNumbers[i]),simplify=F)),
				  rep(names(OutEEGeneNumbers),OutEEGeneNumbers))

GeneNames=paste("G",c(1:NumofGene.raw),sep="_")
names(PhiInput)=GeneNames
names(MeanInput)=GeneNames
#########
# data
#########
EEList=sapply(1:NumofGene.raw, function(j) sapply(1:NumofSample, function(i)suppressWarnings(rnbinom(1,mu=NormFactor[i]*MeanInput[j], size=PhiInput[j]))))

generateDataraw=t(EEList)
DVDSample=sample(DVDLibrary,sum(DEGeneNumbers),replace=T)

DErawNames=vector("list",length(DEGeneNumbers))
st=1
for(i in 1:length(DEGeneNumbers)){
	for(j in st:(st+DEGeneNumbers[i]-1)){
		NumGroup=max(AllParti[names(DEGeneNumbers)[i],])
		SampleGroup=sample(NumGroup,NumGroup)
		DVDSampleEach=c(1,DVDSample[j]^c(1:(NumGroup-1)))
		for(k in 1:NumGroup){
		CondWhich=which(AllParti[names(DEGeneNumbers)[i],]==SampleGroup[k])
		SampleChoose=which(Conditions%in%colnames(AllParti)[CondWhich])
		generateDataraw[j,SampleChoose]=sapply(1:length(SampleChoose), function(i)suppressWarnings(rnbinom(1, size=PhiInput[j], mu=DVDSampleEach[k]*MeanInput[j]*NormFactor[i])),simplify=T)
		}}
		DErawNames[[i]]=GeneNames[st:(st+DEGeneNumbers[i]-1)]
		st=st+DEGeneNumbers[i]
}

rownames(generateDataraw)=GeneNames
MeanVector=rowMeans(generateDataraw)
VarVector=apply(generateDataraw,1,var)
MOV.post=MeanVector/VarVector
EErawNames=GeneNames[!GeneNames%in%unlist(DErawNames)]


### Remove MOV=NA
generateData=generateDataraw
generateData=generateData[!is.na(MOV.post)& MeanVector>2 & MeanVector<10000 ,] 
InName=rownames(generateData)
#print(paste("NA MOV's",sum(is.na(MOV.post)),sum( MeanVector<2), sum(MeanVector>10000)))
## DE
##################################
FinalDEInName=sapply(1:length(DEGeneNumbers),function(i)InName[InName%in%DErawNames[[i]]][1:OutDEGeneNumbers[i]],simplify=F)
FinalEEInName=InName[InName%in%EErawNames][1:OutEEGeneNumbers]
FinalNames=c(unlist(FinalDEInName),FinalEEInName)

generateData=generateData[FinalNames,]
########################################

UseName=rownames(generateData)
phiuse=PhiInput[rownames(generateData)]
meanuse=MeanInput[rownames(generateData)]

OutName=paste("Gene",c(1:nrow(generateData)),sep="_")
names(OutName)=rownames(generateData)
OutData=generateData
rownames(OutData)=as.vector(OutName)
names(OutGenePatterns)=as.vector(OutName)
output=list(data=OutData, Patterns=OutGenePatterns)
}
