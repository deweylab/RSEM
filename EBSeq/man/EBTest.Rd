\name{EBTest}
\alias{EBTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Using EM algorithm to calculate the posterior probabilities of being DE
}
\description{
Base on the assumption of NB-Beta Empirical Bayes model, the EM algorithm is used to get the posterior probability of being DE.
}
\usage{
EBTest(Data, NgVector=NULL, Vect5End=NULL, Vect3End=NULL, Conditions, sizeFactors, maxround,tau=NULL,CI=NULL,CIthre=NULL, Pool=F, NumBin=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{Data}{
A data matrix contains expression values for each transcript .(Gene level or Isoform level.). In which rows should be transcripts and columns should be samples.
}
  \item{NgVector}{
A vector contains the Ng value of each isoform. If the isoform is in a gene with 2 isoforms, Ng should be 2. Ng could be only 1, 2 or 3. If it's gene level data, Ngvector should all be 1. The vector length should be the same as the number of rows in Data.
}
  \item{Vect5End}{
A vector contains the 5' end information of each isoform. It should be 1 if the isoform contains 5' end and otherwise should be 0. If it's gene level data, Vect5End should all be 1. The vector length should be the same as the number of rows in Data.
(Not recommended)
}
  \item{Vect3End}{
A vector contains the 3' end information of each isoform. It should be 1 if the
isoform contains 3' end and otherwise should be 0. If it's gene level data, Vect3End should all be 1. The vector length should be the same as the number of rows in Data.
(Not recommended)
}
  \item{Conditions}{
A vector indicates the condition each sample belongs to. 
}


  \item{sizeFactors}{
The normalization factors. 
The normalization factors could be a vector with lane specitic numbers.
Or it could be a matrix with lane and transcript specific numbers.
}
  \item{maxround}{
Number of iterations. The suggested value is 5.
}

\item{tau}{
The tau value from RSEM output. If the data has no replicates within condition, 
EBSeq will use the CI of tau to capture the variation from mapping
uncertainty and estimate the variance.
	}
\item{CI}{
The CI of each tau from RSEM output	
	}
\item{CIthre}{
The threshold of CI RSEM used.
	}
\item{Pool, NumBin}{
Working without replicates, we should define the Pool=T in the
 EBTest function to enable pooling.
By defining NumBin = 1000, EBSeq will group the genes with similar means
together into 1,000 bins.
With the assumption that no more than 50\% genes are DE in the data set,
We take genes whose FC are in the 25\% - 75\% quantile of the FC's  as the
candidate genes.
For each bin, the bin-wise variance estimation would be the median of the
cross condition variance estimations of the candidate genes within that bin.
We use the cross condition variance estimations for the candidate genes
and the bin-wise variance estimations of the host bin for the non-candidate genes.
}

}

\details{
For each transcript gi within condition, the model assumes:
X_gis|mu_gi ~ NB (r_gi0 * l_s, q_gi)
q_gi|alpha, beta^N_g,b_gi ~ Beta (alpha, neta^N_g,b_gi)
In which the l_s is the sizeFactors of sample s.

The function will test:
H0: q_giC1 = q_giC2
H1: q_giC1 != q_giC2


}
\value{
\item{Alpha }{Fitted parameter alpha of the prior beta distribution. Rows are the values for each iteration.}
\item{Beta }{Fitted parameter beta of the prior beta distribution. Rows are the values for each iteration.}
\item{P, PFromZ }{ The bayes estimator of being DE.Rows are the values for each iteration.}
\item{Z, PoissonZ}{ The Posterior Probability of being DE for each transcript. (Maybe not in the same order of input)}
\item{RList}{ The fitted values of r for each transcript.}
\item{MeanList}{The mean of each transcript. (Cross conditions)}
\item{VarList}{The variance of each transcript. (Cross conditions, using the expression values devided by it's sizeFactors)}
\item{QListi1}{The fitted q values of each transcript within condition 1.}
\item{QListi2}{The fitted q values of each transcript within condition 2.}
\item{C1Mean}{The mean of each transcript within Condition 1}
\item{C2Mean}{The mean of each transcript within Condition 2}
\item{C1EstVar}{The estimated variance of each transcript within Condition 1}
\item{C2EstVar}{The estimated variance of each transcript within Condition 2}
\item{PoolVar}{The variance of each transcript. (The pooled value of within condition EstVar)}
\item{DataList}{A List of data that grouped with Ng and bias.}
\item{PPDE}{The Posterior Probability of being DE for each transcript. (The same order of input)}


}
\references{
}
\author{
Ning Leng
}
\note{
}


\seealso{
}
\examples{
#Simulate Gene level data
GeneGenerate=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, OnlyData="Y")
GeneData=GeneGenerate$data

# Run EBSeq
# sizeFactors could be obtained by MedianNorm, QuantileNorm or RankNorm
EBres=EBTest(Data=GeneData, NgVector=rep(1,10^4), Vect5End=rep(1,10^4), Vect3End=rep(1,10^4), Conditions=as.factor(rep(c(1,2),each=5)), sizeFactors=rep(1,10),maxround=5)

# Isoform Level
IsoGenerate=IsoSimu(DVDconstant=NULL, DVDqt1=.97, DVDqt2=.98, Conditions=as.factor(rep(c(1,2),each=5)), NumofSample=10, NumofIso=c(1000,2000,3000), DEIsoProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, OnlyData=T )

IsoMat=do.call(rbind,IsoGenerate$data)
IsoNames=rownames(IsoMat)

Ngvector=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=Ngvector$IsoformNgTrun

IsoEBres=EBTest(Data=IsoMat, NgVector=IsoNgTrun, Conditions=as.factor(rep(c(1,2),each=5)),sizeFactors=rep(1,10), maxround=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
