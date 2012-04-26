\name{PolyFitPlot}
\alias{PolyFitPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fit the mean-var relationship using polynomial regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PolyFitPlot(X, Y, nterms, xname = "Estimated Mean", yname = "Estimated Var", pdfname = "", xlim =  c(-1,5), ylim = c(-1,7), ChangeXY = F, col = "red")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The first group of values want to be fitted by the polynomial regression. ( e.g Mean of the data. )
}
  \item{Y}{
The second group of values want to be fitted by the polynomial regression. ( e.g. variance of the data.) The length of Y should be the same as the length of X.
}
  \item{nterms}{
How many polynomial terms want to be used.
}
  \item{xname}{
Name of the x axis.
}
  \item{yname}{
Name of the y axis.
}
  \item{pdfname}{
Name of the plot.
}
  \item{xlim}{
The x limits of the plot. 
}
  \item{ylim}{
The y limits of the plot.

}
  \item{ChangeXY}{
If ChangeXY is setted to be TRUE, X will be treated as the dependent variable and Y will be treated as the independent one. Default is FALSE.
}
  \item{col}{
Color of the fitted line.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	The PolyFitPlot function provides a smooth scatter plot of two variables and their best fitting line of polynomial regression.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ning Leng
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
GeneData=GeneSimu(DVDconstant=4, DVDqt1=NULL, DVDqt2=NULL, Conditions=rep(c(1,2),each=5), NumofSample=10, NumofGene=10000, DEGeneProp=.1, Phiconstant=NULL, Phi.qt1=.25, Phi.qt2=.75, Meanconstant=NULL, OnlyData="Y")

EBres=NBBetaEB.bias.uniqueP_PoolVarSpeedUp_MDFPoi(Data=GeneData,NgVector=rep(1,10^4), Vect5End=rep(1,10^4), Vect3End=rep(1,10^4), Conditions=as.factor(rep(c(1,2),each=5)), maxround=5)

poly=PolyFitPlot(X=EBres$MeanList[[1]], Y=EBres$PoolVar[[1]], nterms=5, xname = "mean", yname = "var", pdfname=NULL, xlim = c(0,4.5),ylim = c(-2,8), ChangeXY = F, col = "red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
