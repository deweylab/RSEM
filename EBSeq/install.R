#!/usr/bin/env Rscript

result <- suppressWarnings(tryCatch({
       	  library("EBSeq", lib.loc = ".")
       }, error = function(err) {
            tryCatch({
		source("http://www.bioconductor.org/biocLite.R") 
    		biocLite("EBSeq", lib = ".")
    		library("EBSeq", lib.loc = ".")
		}, error = function(err) {
		     tryCatch({
		         cat("Failed to install the latest version of EBSeq from Bioconductor! Try to install EBSeq v1.1.5 locally instead.\n")
      			 install.packages(c("blockmodeling_0.1.8.tar.gz", "EBSeq_1.1.5.tar.gz"), lib = ".", repos = NULL)
      			 library("EBSeq", lib.loc = ".") 
			 }, error = function(err) { cat("Failed to install EBSeq v1.1.5 locally!\n") })
		    })
	}))


