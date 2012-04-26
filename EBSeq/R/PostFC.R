PostFC=function(EBoutput) {
	GeneRealMeanC1=unlist(EBoutput$C1Mean)
	GeneRealMeanC2=unlist(EBoutput$C2Mean)
	GeneRealMean=(GeneRealMeanC1+GeneRealMeanC2)/2

	GeneRealFC=GeneRealMeanC1/GeneRealMeanC2

	GeneR=unlist(EBoutput$RList)
	GeneR[GeneR<=0 | is.na(GeneR)]=GeneRealMean[GeneR<=0 | is.na(GeneR)]*.99/.01

	GeneAlpha=EBoutput[[1]][nrow(EBoutput[[1]]),]
	GeneBeta=unlist(sapply(1:length(EBoutput$C1Mean),function(i)rep(EBoutput[[2]][nrow(EBoutput[[1]]),i],length(EBoutput$C1Mean[[i]]))))
	GeneBeta=as.vector(GeneBeta)
	# Post alpha = alpha + r_C1 * 3
	# Post beta = beta + Mean_C1 * 3
	# Post Mean of q in C1 P_q_C1= P_a/ (P_a + P_b)
	# Post FC = (1-p_q_c1)/p_q_c1 /( (1-p_q_c2)/p_q_c2)

	GenePostAlpha=GeneAlpha+3*GeneR
	GenePostBetaC1=GeneBeta+3*GeneRealMeanC1
	GenePostBetaC2=GeneBeta+3*GeneRealMeanC2
	GenePostQC1=GenePostAlpha/(GenePostAlpha+GenePostBetaC1)
	GenePostQC2=GenePostAlpha/(GenePostAlpha+GenePostBetaC2)

	GenePostFC=((1-GenePostQC1)/(1-GenePostQC2))*(GenePostQC2/GenePostQC1)
	Out=list(GenePostFC=GenePostFC, GeneRealFC=GeneRealFC)

}
