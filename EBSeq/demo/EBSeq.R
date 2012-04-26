library(EBSeq)
set.seed(13)

# Section 3.1

GeneGenerate=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL,
  Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000,
  DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.1, Phi.qt2=.9,
  Meanconstant=NULL, OnlyData=T)
GeneData=GeneGenerate$data
GeneTrueDENames=GeneGenerate$TrueDE
str(GeneData)
str(GeneTrueDENames)

Sizes=MedianNorm(GeneData)

EBres=EBTest(Data=GeneData, 
  Conditions=as.factor(rep(c(1,2),each=5)),sizeFactors=Sizes, maxround=5)

PP=GetPP(EBres)
str(PP)
DEfound=names(PP)[which(PP>=.95)]
str(DEfound)
sum(DEfound%in%GeneTrueDENames)

QQP(QList=EBres$QList1, AlphaResult=EBres[[1]][5,1], 
  BetaResult=EBres[[2]][5,1], name="Gene Simulation", AList="F", GroupName=NULL)
DenNHist(QList=EBres$QList1, Alpha=EBres[[1]][5,1], Beta=EBres[[2]][5,1], 
  name="Gene Simulation", AList="F", GroupName=NULL)

# Section 3.2

IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, 
  Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, 
  NumofIso=c(1000,2000,3000), DEIsoProp=.1, Phiconstant=NULL, 
  Phi.qt1=.25, Phi.qt2=.75, OnlyData=T )
str(IsoGenerate)

IsoMat=do.call(rbind,IsoGenerate$data)
str(IsoMat)

IsoSizes=MedianNorm(IsoMat)

IsoNames=rownames(IsoMat)
str(IsoNames)
GeneNames=paste("Gene",c(1:3000),sep="_")
IsosGeneNames=c(GeneNames[1:1000],rep(GeneNames[1001:2000],each=2),
  rep(GeneNames[2001:3000],each=3))
NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,1001:1003,3001:3003)]

IsoEBres=EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
  Conditions=as.factor(rep(c(1,2),each=5)),sizeFactors=IsoSizes, maxround=5)
IsoPP=GetPP(IsoEBres)
str(IsoPP)
IsoDE=IsoPP[which(IsoPP>=.95)]
str(IsoDE)
sum(names(IsoDE)%in%IsoGenerate$TrueDE)

par(mfrow=c(2,2))
PolyFitValue=vector("list",3)
for(i in 1:3)
  PolyFitValue[[i]]=PolyFitPlot(IsoEBres$C1Mean[[i]], 
    IsoEBres$C1EstVar[[i]],5)

PolyAll=PolyFitPlot(unlist(IsoEBres$C1Mean), unlist(IsoEBres$C1EstVar),5)
lines(log10(IsoEBres$C1Mean[[1]][PolyFitValue[[1]]$sort]), 
  PolyFitValue[[1]]$fit[PolyFitValue[[1]]$sort],col="yellow")
lines(log10(IsoEBres$C1Mean[[2]][PolyFitValue[[2]]$sort]), 
  PolyFitValue[[2]]$fit[PolyFitValue[[2]]$sort],col="pink")
lines(log10(IsoEBres$C1Mean[[3]][PolyFitValue[[3]]$sort]), 
  PolyFitValue[[3]]$fit[PolyFitValue[[3]]$sort],col="green")
legend("topleft",c("All Isoforms","Ng = 1","Ng = 2","Ng = 3"),
  col=c("red","yellow","pink","green"),lty=1,lwd=3,box.lwd=2)

par(mfrow=c(2,2))
QQP(QList=IsoEBres$QList1, AlphaResult=IsoEBres[[1]][5,],
 BetaResult=IsoEBres[[2]][5,], 
 name="Isoforms", AList="F", GroupName=paste("Ng = ",c(1:3),sep=""))

DenNHist(QList=IsoEBres$QList1, Alpha=IsoEBres[[1]][5,], 
  Beta=IsoEBres[[2]][5,], 
  name="Isoforms", AList="F", GroupName=paste("Ng = ",c(1:3),sep=""))

# Section 3.3

Conditions=c("C1","C1","C2","C2","C3","C3")
PosParti=GetPatterns(Conditions)
PosParti

Parti=PosParti[-3,]
Parti

MultiData=GeneMultiSimu(Conditions=Conditions,AllParti=Parti,
          NumofSample=6,NumofGene=1000,DEGeneProp=c(.7,.1,.1,.1),
          DVDqt1=.98,DVDqt2=.99,Phi.qt1=.25,Phi.qt2=.75)
str(MultiData)

MultiSize=MedianNorm(MultiData$data)
MultiRes=EBMultiTest(MultiData$data,NgVector=NULL,Conditions=Conditions,
           AllParti=Parti, sizeFactors=MultiSize, maxround=5)
MultiPP=GetMultiPP(MultiRes)
names(MultiPP)
MultiPP$PP[1:10,]
MultiPP$MAP[1:10]
MultiPP$Patterns
sum(MultiPP$MAP==MultiData$Patterns)

# EOF