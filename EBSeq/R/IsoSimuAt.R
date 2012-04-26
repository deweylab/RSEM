IsoSimuAt<-function(DVDconstant=NULL, DVDqt1=NULL, DVDqt2=NULL, Conditions, NumofSample, NumofIso=NULL, DEIsoProp, Phiconstant=NULL, Phi.qt1=NULL, Phi.qt2=NULL,NormFactor=NULL, OnlyData=T)
{
#Ng paired 2012 feb 1
if(is.null(NormFactor)) NormFactor=rep(1,NumofSample)
data(IsoEBresultGouldBart2)

MeansC1=IsoEBresultGouldBart2$C1Mean
MeansC2=IsoEBresultGouldBart2$C2Mean
MeanDVD=sapply(1:9,function(i) MeansC1[[i]]/MeansC2[[i]])
if (length(DVDconstant)==0) DVDLibrary= unlist(MeanDVD)[unlist(MeanDVD)<quantile(unlist(MeanDVD)[unlist(MeanDVD)!=Inf],DVDqt2) & unlist(MeanDVD)>quantile(unlist(MeanDVD)[unlist(MeanDVD)!=Inf],DVDqt1)]




VarInput=IsoEBresultGouldBart2$VarList
VarInputNg=list(VarInput[[1]],unlist(VarInput[c(2,4,6,8)]),unlist(VarInput[c(3,5,7,9)]))

if(length(NumofIso)==0) NumofIso=sapply(1:3,function(i)length(VarInputNg[[i]]))
PhiInput.raw=IsoEBresultGouldBart2$RList
PhiInput.raw.Ng=list(PhiInput.raw[[1]],unlist(PhiInput.raw[c(2,4,6,8)]),unlist(PhiInput.raw[c(3,5,7,9)]))


if (length(Phiconstant)==0){
	PhiLibrary=sapply(1:3,function(i)PhiInput.raw.Ng[[i]][1/PhiInput.raw.Ng[[i]]<quantile(1/PhiInput.raw.Ng[[i]],Phi.qt2) & 1/PhiInput.raw.Ng[[i]]>quantile(1/PhiInput.raw.Ng[[i]],Phi.qt1)],simplify=F)
	PhiIndex=sapply(1:3, function(i)sample(names(PhiLibrary[[i]]),NumofIso[[i]],replace=T),simplify=F)
	PhiInputNg=sapply(1:3, function(i)PhiLibrary[[i]][PhiIndex[[i]]])
}
if (length(Phiconstant)!=0)PhiInputNg=sapply(1:3,function(i)rep(Phiconstant,NumofIso[[i]]),simplify=F)

# Wanna DENumbers be proportion to 2 
DEIsoNumbers=round(NumofIso*DEIsoProp/2)*2
IsoNames=sapply(1:3,function(i)paste("I",i,c(1:NumofIso[i]),sep="_"),simplify=F)
MeanNg=list(IsoEBresultGouldBart2$MeanList[[1]],unlist(IsoEBresultGouldBart2$MeanList[c(2,4,6,8)]),
unlist(IsoEBresultGouldBart2$MeanList[c(3,5,7,9)]))
MeanInputNg=sapply(1:3, function(i)MeanNg[[i]][PhiIndex[[i]]])

for(i in 1:3){
	names(MeanInputNg[[i]])=IsoNames[[i]]
	names(PhiInputNg[[i]])=IsoNames[[i]]
	}

#########
# data
#########
EEList=sapply(1:3,function(i) sapply(1:NumofIso[[i]], function(j)sapply(1:NumofSample,function(h) rnbinom(1,mu=MeanInputNg[[i]][j]*NormFactor[h], size=PhiInputNg[[i]][j]))),simplify=F)


generateDataraw=vector("list",3)
MeanVector=vector("list",3)
VarVector=vector("list",3)
MOV.post=vector("list",3)


for(g in 1:3){
    generateDataraw[[g]]=t(EEList[[g]][,1:NumofIso[g]])
	if(length(DVDconstant)==0){
		for(j in 1:NumofIso[g]){
	    	 if (j<=(DEIsoNumbers[g]/2)) generateDataraw[[g]][j,((NumofSample/2)+1):NumofSample]=sapply((NumofSample/2+1):NumofSample, function(h)rnbinom(1, size=PhiInputNg[[g]][j], mu=sample(DVDLibrary,1)*MeanInputNg[[g]][j]*NormFactor[h]), simplify=T)
	     	if (j>=((DEIsoNumbers[g]/2)+1) & j <=DEIsoNumbers[g]) generateDataraw[[g]][j,1:(NumofSample/2)]=sapply(1:(NumofSample/2),function(h) rnbinom(1, size=MeanInputNg[[g]][j], mu= sample(DVDLibrary,1)*MeanInputNg[[g]][j]*NormFactor[h]),simplify=T)
}
	 }
	if(length(DVDconstant)!=0){
        for(j in 1:NumofIso[g]){
             if (j<=(DEIsoNumbers[g]/2)) generateDataraw[[g]][j,((NumofSample/2)+1):NumofSample]=sapply((NumofSample/2+1):NumofSample, function(h)rnbinom(1, DVDconstant*MeanInputNg[[g]][j]*NormFactor[h]),simplify=T)
             if (j>=((DEIsoNumbers[g]/2)+1) & j <=DEIsoNumbers[g]) generateDataraw[[g]][j,1:(NumofSample/2)]=sapply(1:(NumofSample/2),function(h) rnbinom(1, DVDconstant*MeanInputNg[[g]][j]*NormFactor[h]),simplify=T)
		}
	}
rownames(generateDataraw[[g]])=IsoNames[[g]][1:NumofIso[g]]
MeanVector[[g]]=rowMeans(generateDataraw[[g]])
VarVector[[g]]=apply(generateDataraw[[g]],1,var)
MOV.post[[g]]=MeanVector[[g]]/VarVector[[g]]
}


### Remove MOV=NA
generateData=generateDataraw
for (i in 1:3) generateData[[i]]=generateData[[i]][!is.na(MOV.post[[i]]),] 
print(paste("NA MOV's",sum(is.na(unlist(MOV.post)))))
#tmpmean=sapply(1:9,function(i)rowMeans(generateData[[i]]))
#tmpvar=sapply(1:9,function(i)apply(generateData[[i]],1,var))
#source("plot_functions.R")
#CheckSimuNg(tmpmean,tmpvar,c(-1,5),c(-1,7))




## DE
UseName=sapply(1:3, function(i)rownames(generateData[[i]]),simplify=F)
TrueDE=sapply(1:3, function(i)UseName[[i]][UseName[[i]] %in% rownames(generateData[[i]])[1:DEIsoNumbers[i]]],simplify=F)
TrueDE.unlist=do.call(c,TrueDE)

TrueDELength=sapply(TrueDE,length)

AtNames_Level=vector("list",4)
AtLoc=vector("list",3)
AtFold=vector("list",3)
names(AtNames_Level)=c(4,6,8,10)


for(j in 1:3){
AtLoc[[j]]=sample(c(1:length(Conditions)), TrueDELength[j], replace=T)
AtFold[[j]]=sample(c(4,6,8,10),TrueDELength[j], replace=T)

for(i in 1:TrueDELength[j]){

generateData[[j]][(TrueDELength[j]+i),AtLoc[[j]][i]]=generateData[[j]][(TrueDELength[j]+i),AtLoc[[j]][i]]*AtFold[[j]][i]
AtNames_Level[[as.character(AtFold[[j]][i])]]=c(AtNames_Level[[as.character(AtFold[[j]][i])]],rownames(generateData[[j]])[TrueDELength[j]+i])
}
}
phiuse=sapply(1:3,function(i)PhiInputNg[[i]][UseName[[i]]])
meanuse=sapply(1:3,function(i)MeanInputNg[[i]][UseName[[i]]])

#generateDataNg=list(generateData[[1]], do.call(rbind,generateData[c(2,4,6,8)]), do.call(rbind,generateData[c(3,5,7,9)]))
generateDataNg=generateData

#if(OnlyData==T){

OutName=sapply(1:3,function(i)paste("Iso",i,c(1:nrow(generateDataNg[[i]])),sep="_"))
for(i in 1:3)names(OutName[[i]])=rownames(generateDataNg[[i]])
OutData=generateDataNg
for(i in 1:3)rownames(OutData[[i]])=as.vector(OutName[[i]])
OutTrueDE=as.vector(unlist(OutName)[TrueDE.unlist])
OutAt=as.vector(unlist(OutName)[AtNames <- Level])

output=list(data=OutData, TrueDE=OutTrueDE, Outliers=OutAt)
#	return(output)
#    }
	}
