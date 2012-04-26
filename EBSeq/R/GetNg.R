GetNg<- function(IsoformName, GeneName){
	GeneNg = tapply(IsoformName, GeneName, length)
	IsoformNg = GeneNg[GeneName]
	names(IsoformNg) = IsoformName
	GeneNgTrun=GeneNg
	GeneNgTrun[GeneNgTrun>3]=3
	IsoformNgTrun=IsoformNg
	IsoformNgTrun[IsoformNgTrun>3]=3
	out=list( GeneNg=GeneNg, GeneNgTrun=GeneNgTrun, IsoformNg=IsoformNg, IsoformNgTrun=IsoformNgTrun)
	}
