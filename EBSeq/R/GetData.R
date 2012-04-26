GetData <-
function(path,Name1,Name2,type)
{
Data=vector("list",8)
Filenames=NULL
Tablenames=NULL
for (name in 1:4)
	{
		if (type=="I")
			Filenames=c(Filenames,paste(path,Name1,name,"_isoform_nus.tab",sep=""))  
		if (type=="G")	
			Filenames=c(Filenames,paste(path,Name1,name,"_gene_nus.tab",sep=""))  
		Tablenames=c(Tablenames,paste(Name1,name,sep=""))
	}
for (name in 1:4)
	{
		if (type=="I")
			Filenames=c(Filenames,paste(path,Name2,name,"_isoform_nus.tab",sep=""))
		if (type=="G")
			Filenames=c(Filenames,paste(path,Name2,name,"_gene_nus.tab",sep=""))
		Tablenames=c(Tablenames,paste(Name2,name,sep=""))
	}


names(Data)=Tablenames
for (file in 1:8)
	{
		temp=read.table(Filenames[file],header=T)
		temp2=as.matrix(temp[-1])
		rownames(temp2)=as.vector(as.matrix(temp[1]))
		Data[[file]]=temp2
	}
	Data
}

