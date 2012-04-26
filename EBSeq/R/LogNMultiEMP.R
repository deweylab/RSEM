LogNMulti <-
function(Input, InputSP, EmpiricalR, EmpiricalRSP, NumOfEachGroup, AlphaIn, BetaIn,  PIn, NoneZeroLength, AllParti, Conditions)
{

        #For each gene (m rows of Input---m genes)
        #Save each gene's F0, F1 for further likelihood calculation. 
 		FList=sapply(1:nrow(AllParti),function(i)sapply(1:nlevels(as.factor(AllParti[i,])),
			     	   function(j)f0(do.call(cbind,InputSP[AllParti[i,]==j]),AlphaIn, BetaIn, 
									 do.call(cbind,EmpiricalRSP[AllParti[i,]==j]), NumOfEachGroup, log=T)),
					  simplify=F) 
		FPartiLog=sapply(FList,rowSums)
		FMat=exp(FPartiLog)
		rownames(FMat)=rownames(Input)
        #Get z
		#Use data.list in logfunction
        PInMat=matrix(rep(1,nrow(Input)),ncol=1)%*%matrix(PIn,nrow=1)
		FmultiP=FMat*PInMat
		Denom=rowSums(FmultiP)
		ZEach=apply(FmultiP,2,function(i)i/Denom)
		zNaNName1=names(Denom)[is.na(Denom)]
		# other NAs in LikeFun
		LF=ZEach*(log(FmultiP))
		zNaNMore=rownames(LF)[which(is.na(rowSums(LF)))]
		zNaNName=unique(c(zNaNName1,zNaNMore))
		zGood=which(!rownames(LF)%in%zNaNName)
		ZEachGood=ZEach[zGood,]
		###Update P
        PFromZ=colSums(ZEach[zGood,])/length(zGood)
        FGood=FMat[zGood,]
		### MLE Part ####
        # Since we dont wanna update p and Z in this step
        # Each Ng for one row
		
		NumGroupVector=rep(c(1:NoneZeroLength),NumOfEachGroup)
		
		NumGroupVector.zGood=NumGroupVector[zGood]
		NumOfEachGroup.zGood=tapply(NumGroupVector.zGood,NumGroupVector.zGood,length)

        StartValue=c(AlphaIn, BetaIn,PIn[-1])
		InputSPGood=sapply(1:length(InputSP),function(i)InputSP[[i]][zGood,],simplify=F)
        EmpiricalRSPGood=sapply(1:length(EmpiricalRSP),function(i)EmpiricalRSP[[i]][zGood,],simplify=F)

		Result<-optim(StartValue,LikefunMulti,InputPool=list(InputSPGood,Input[zGood,],ZEach[zGood,], 
					 NoneZeroLength,EmpiricalR[zGood, ],EmpiricalRSPGood, NumOfEachGroup.zGood, AllParti))
		AlphaNew= Result$par[1]
		BetaNew=Result$par[2:(1+NoneZeroLength)]
        PNewNo1=Result$par[(2+NoneZeroLength):length(Result$par)]
		PNew=c(1-sum(PNewNo1),PNewNo1)
		##
        Output=list(AlphaNew=AlphaNew,BetaNew=BetaNew,PNew=PNew,ZEachNew=ZEach, ZEachGood=ZEachGood, 
					PFromZ=PFromZ, zGood=zGood, zNaNName=zNaNName,FGood=FGood)
        Output
    }

