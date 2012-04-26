IsoSimu=function(DVDconstant=NULL, DVDqt1=NULL, DVDqt2=NULL, Conditions, NumofSample, NumofIso=NULL, DEIsoProp, Phiconstant=NULL, Phi.qt1=NULL, Phi.qt2=NULL,NormFactor=NULL, OnlyData=T)
{
# 2012 feb 1 
# paired simulation
data(IsoEBresultGouldBart2)
if(is.null(NormFactor)) NormFactor=rep(1,NumofSample)

MeansC1=IsoEBresultGouldBart2$C1Mean
MeansC2=IsoEBresultGouldBart2$C2Mean
MeanDVD=sapply(1:9,function(i) MeansC1[[i]]/MeansC2[[i]])
# DVD library with each group here
if (length(DVDconstant)==0) DVDLibrary= unlist(MeanDVD)[unlist(MeanDVD)<quantile(unlist(MeanDVD)[unlist(MeanDVD)!=Inf],DVDqt2) & unlist(MeanDVD)>quantile(unlist(MeanDVD)[unlist(MeanDVD)!=Inf],DVDqt1)]



# If DVD constant, use constant when generate
# If not, use DVDLibrary

VarInput=IsoEBresultGouldBart2$VarList
VarInputNg=list(VarInput[[1]],unlist(VarInput[c(2,4,6,8)]),unlist(VarInput[c(3,5,7,9)]))
#If NumofIso=NULL, empirical # of Iso
#If !=NULL , Input a 9-vector
if(length(NumofIso)==0) NumofIso.raw=sapply(1:3,function(i)length(VarInputNg[[i]]))
if(length(NumofIso)!=0) NumofIso.raw=NumofIso*2

PhiInput.raw=IsoEBresultGouldBart2$RList
PhiInput.raw.Ng=list(PhiInput.raw[[1]],unlist(PhiInput.raw[c(2,4,6,8)]),unlist(PhiInput.raw[c(3,5,7,9)]))


if (length(Phiconstant)==0){
	PhiLibrary=sapply(1:3,function(i)PhiInput.raw.Ng[[i]][1/PhiInput.raw.Ng[[i]]<quantile(1/PhiInput.raw.Ng[[i]],Phi.qt2) & 1/PhiInput.raw.Ng[[i]]>quantile(1/PhiInput.raw.Ng[[i]],Phi.qt1)],simplify=F)
	PhiIndex=sapply(1:3, function(i)sample(names(PhiLibrary[[i]]),NumofIso.raw[[i]],replace=T),simplify=F)
	PhiInputNg=sapply(1:3, function(i)PhiLibrary[[i]][PhiIndex[[i]]])
}
if (length(Phiconstant)!=0)PhiInputNg=sapply(1:3,function(i)rep(Phiconstant,NumofIso.raw[[i]]),simplify=F)

# Wanna DENumbers be proportion to 2 
DEIsoNumbers=round(NumofIso.raw*DEIsoProp/2)*2
IsoNames=sapply(1:3,function(i)paste("I",i,c(1:NumofIso.raw[i]),sep="_"),simplify=F)
MeanNg=list(IsoEBresultGouldBart2$MeanList[[1]],unlist(IsoEBresultGouldBart2$MeanList[c(2,4,6,8)]),
unlist(IsoEBresultGouldBart2$MeanList[c(3,5,7,9)]))
MeanInputNg=sapply(1:3, function(i)MeanNg[[i]][PhiIndex[[i]]])

for(i in 1:3){
	names(MeanInputNg[[i]])=IsoNames[[i]]
	names(PhiInputNg[[i]])=IsoNames[[i]]
	}

##############################
# Get Ng version to every one
##############################


#########
# data
#########
EEList=sapply(1:3,function(i) sapply(1:NumofIso.raw[[i]], function(j)sapply(1:NumofSample,function(h) rnbinom(1,mu=MeanInputNg[[i]][j]*NormFactor[h], size=PhiInputNg[[i]][j]))),simplify=F)


generateDataraw=vector("list",3)
MeanVector=vector("list",3)
VarVector=vector("list",3)
MOV.post=vector("list",3)


for(g in 1:3){
    generateDataraw[[g]]=t(EEList[[g]][,1:NumofIso.raw[g]])
	if(length(DVDconstant)==0){
		for(j in 1:NumofIso.raw[g]){
	    	 if (j<=(DEIsoNumbers[g]/2)) generateDataraw[[g]][j,((NumofSample/2)+1):NumofSample]=sapply((NumofSample/2+1):NumofSample, function(h)suppressWarnings(rnbinom(1, size=PhiInputNg[[g]][j], mu=sample(DVDLibrary,1)*MeanInputNg[[g]][j]*NormFactor[h])), simplify=T)
	     	if (j>=((DEIsoNumbers[g]/2)+1) & j <=DEIsoNumbers[g]) generateDataraw[[g]][j,1:(NumofSample/2)]=sapply(1:(NumofSample/2),function(h) suppressWarnings(rnbinom(1, size=MeanInputNg[[g]][j], mu= sample(DVDLibrary,1)*MeanInputNg[[g]][j]*NormFactor[h])),simplify=T)
}
	 }
	if(length(DVDconstant)!=0){
        for(j in 1:NumofIso.raw[g]){
             if (j<=(DEIsoNumbers[g]/2)) generateDataraw[[g]][j,((NumofSample/2)+1):NumofSample]=sapply((NumofSample/2+1):NumofSample, function(h)suppressWarnings(rnbinom(1, DVDconstant*MeanInputNg[[g]][j]*NormFactor[h])),simplify=T)
             if (j>=((DEIsoNumbers[g]/2)+1) & j <=DEIsoNumbers[g]) generateDataraw[[g]][j,1:(NumofSample/2)]=sapply(1:(NumofSample/2),function(h) wuppressWarnings(rnbinom(1, DVDconstant*MeanInputNg[[g]][j]*NormFactor[h])),simplify=T)
		}
	}
rownames(generateDataraw[[g]])=IsoNames[[g]][1:NumofIso.raw[g]]
MeanVector[[g]]=rowMeans(generateDataraw[[g]])
VarVector[[g]]=apply(generateDataraw[[g]],1,var)
MOV.post[[g]]=MeanVector[[g]]/VarVector[[g]]
}


### Remove MOV=NA
generateData=generateDataraw
for (i in 1:3) generateData[[i]]=generateData[[i]][!is.na(MOV.post[[i]]),] 
#print(paste("NA MOV's",sum(is.na(unlist(MOV.post)))))
NumDENow=sapply(1:3, function(i)sum(rownames(generateData[[i]])%in%rownames(generateDataraw[[i]])[1:DEIsoNumbers[i]]))

if(length(NumofIso)!=0){
	    for(i in 1:3)
		generateData[[i]]=generateData[[i]][c(sample(1:NumDENow[i],round(NumofIso[i]*DEIsoProp),replace=F),round( (dim(generateData[[i]])[1]+1-NumofIso[i]*(1-DEIsoProp)):dim(generateData[[i]])[1])),]
}
generateDataNg=generateData

## DE
UseName=sapply(1:3, function(i)rownames(generateData[[i]]),simplify=F)
TrueDE=sapply(1:3, function(i)UseName[[i]][UseName[[i]] %in% rownames(generateDataraw[[i]])[1:DEIsoNumbers[i]]],simplify=F)
TrueDE.unlist=do.call(c,TrueDE)

phiuse=sapply(1:3,function(i)PhiInputNg[[i]][UseName[[i]]])
meanuse=sapply(1:3,function(i)MeanInputNg[[i]][UseName[[i]]])

#if(OnlyData==T){
    
OutName=sapply(1:3,function(i)paste("Iso",i,c(1:nrow(generateDataNg[[i]])),sep="_"))
for(i in 1:3)names(OutName[[i]])=rownames(generateDataNg[[i]])
OutData=generateDataNg
for(i in 1:3)rownames(OutData[[i]])=as.vector(OutName[[i]])
OutTrueDE=as.vector(unlist(OutName)[TrueDE.unlist])
output=list(data=OutData, TrueDE=OutTrueDE)


#output=list(data=generateDataNg, TrueDE=TrueDE.unlist)
return(output)
#    }
# Now only OnlyData=T version
}

