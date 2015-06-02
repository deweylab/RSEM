# run_spp.R
# =============
# Author: Anshul Kundaje, Computer Science Dept., MIT
# Email: anshul@kundaje.net
# Last updated: Feb 12, 2012
# =============
# MANDATORY ARGUMENTS
# -c=<ChIP_tagAlign/BAMFile>, full path and name of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
# MANDATORY ARGUMENT FOR PEAK CALLING
# -i=<Input_tagAlign/BAMFile>, full path and name of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)
# OPTIONAL ARGUMENTS
# -s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-500:5:1500
# -speak=<strPeak>, user-defined cross-correlation peak strandshift
# -x=<min>:<max>, strand shifts to exclude (This is mainly to avoid phantom peaks) default=10:(readlen+10)
# -p=<nodes> , number of parallel processing nodes, default=NULL
# -fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling
# -npeak=<numPeaks>, threshold on number of peaks to call
# -tmpdir=<tempdir> , Temporary directory (if not specified R function tempdir() is used)
# -filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name
# OUTPUT PARAMETERS
# -odir=<outputDirectory> name of output directory (If not set same as ChIP file directory is used)
# -savn=<narrowpeakfilename> OR -savn NarrowPeak file name
# -savr=<regionpeakfilename> OR -savr RegionPeak file name
# -savd=<rdatafile> OR -savd , save Rdata file
# -savp=<plotdatafile> OR -savp , save cross-correlation plot
# -out=<resultfile>, append peakshift result to a file 
#      format:Filename<tab>numReads<tab>estFragLen<tab>corr_estFragLen<tab>PhantomPeak<tab>corr_phantomPeak<tab>argmin_corr<tab>min_corr<tab>Normalized SCC (NSC)<tab>Relative SCC (RSC)<tab>QualityTag
# -rf , if plot or rdata or narrowPeak file exists replace it. If not used then the run is aborted if the plot or Rdata or narrowPeak file exists
# -clean, if present will remove the original chip and control files after reading them in. CAUTION: Use only if the script calling run_spp.R is creating temporary files

args <- commandArgs(trailingOnly=TRUE); # Read Arguments from command line
nargs = length(args); # number of arguments

# ###########################################################################
# AUXILIARY FUNCTIONS
# ###########################################################################

print.usage <- function() {
# ===================================
# Function will print function usage
# ===================================	
	cat('Usage: Rscript run_spp.R <options>\n',file=stderr())
	cat('MANDATORY ARGUMENTS\n',file=stderr())
	cat('-c=<ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped)(FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz) \n',file=stderr())
	cat('MANDATORY ARGUMENTS FOR PEAK CALLING\n',file=stderr())
	cat('-i=<Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz) \n',file=stderr())
	cat('OPTIONAL ARGUMENTS\n',file=stderr())
	cat('-s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-500:5:1500\n',file=stderr())
	cat('-speak=<strPeak>, user-defined cross-correlation peak strandshift\n',file=stderr())
	cat('-x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10)\n',file=stderr())
	cat('-p=<nodes> , number of parallel processing nodes, default=0\n',file=stderr())
	cat('-fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling\n',file=stderr())
	cat('-npeak=<numPeaks>, threshold on number of peaks to call\n',file=stderr())    
	cat('-tmpdir=<tempdir> , Temporary directory (if not specified R function tempdir() is used)\n',file=stderr())
	cat('-filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name\n',file=stderr())
	cat('OUTPUT ARGUMENTS\n',file=stderr())
	cat('-odir=<outputDirectory> name of output directory (If not set same as ChIP file directory is used)\n',file=stderr())
	cat('-savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)\n',file=stderr())
	cat('-savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width peaks with regions of enrichment)\n',file=stderr())
	cat('-savd=<rdatafile> OR -savd, save Rdata file\n',file=stderr())
	cat('-savp=<plotdatafile> OR -savp, save cross-correlation plot\n',file=stderr())
	cat('-out=<resultfile>, append peakshift/phantomPeak results to a file\n',file=stderr())
	cat('     format:Filename<tab>numReads<tab>estFragLen<tab>corr_estFragLen<tab>PhantomPeak<tab>corr_phantomPeak<tab>argmin_corr<tab>min_corr<tab>Normalized SCC (NSC)<tab>Relative SCC (RSC)<tab>QualityTag)\n',file=stderr())
	cat('-rf, if plot or rdata or narrowPeak file exists replace it. If not used then the run is aborted if the plot or Rdata or narrowPeak file exists\n',file=stderr())
	cat('-clean, if present will remove the original chip and control files after reading them in. CAUTION: Use only if the script calling run_spp.R is creating temporary files\n',file=stderr())
} # end: print.usage()

get.file.parts <- function(file.fullpath) {
# ===================================
# Function will take a file name with path and split the file name into 
# path, fullname, name and ext
# ===================================	
	if (! is.character(file.fullpath)) {
		stop('File name must be a string')
	}
	
	file.parts <- strsplit(as.character(file.fullpath), .Platform$file.sep, fixed=TRUE)[[1]] # split on file separator
	
	if (length(file.parts) == 0) { # if empty file name
		return(list(path='',
						fullname='',
						name='',
						ext='')
		)
	} else {
		if (length(file.parts) == 1) { # if no path then just the file name itself
			file.path <- '.'
			file.fullname <- file.parts
		} else {
			file.path <- paste(file.parts[1:(length(file.parts)-1)], collapse=.Platform$file.sep) # 1:last-1 token is path
			file.fullname <- file.parts[length(file.parts)] # last token is filename
		}        
		file.fullname.parts <- strsplit(file.fullname,'.',fixed=TRUE)[[1]] # split on .
		if (length(file.fullname.parts) == 1) { # if no extension
			file.ext <- ''
			file.name <- file.fullname.parts
		} else {
			file.ext <- paste('.', file.fullname.parts[length(file.fullname.parts)], sep="") # add the . to the last token
			file.name <- paste(file.fullname.parts[1:(length(file.fullname.parts)-1)], collapse=".")
		}
		return(list(path=file.path,
						fullname=file.fullname,
						name=file.name,
						ext=file.ext))
	}         	
} # end: get.file.parts()

parse.arguments <- function(args) {
# ===================================
# Function will parse arguments
# ===================================	
	# Set arguments to default values
	chip.file <- NA  # main ChIP tagAlign/BAM file name
	isurl.chip.file <- FALSE # flag indicating whether ChIP file is a URL
	control.file <- NA # control tagAlign/BAM file name
	isurl.control.file <- FALSE # flag indicating whether control file is a URL   
	sep.min <- -500  # min strand shift
	sep.max <- 1500  # max strand shift
	sep.bin <- 5    # increment for strand shift
	sep.peak <- NA # user-defined peak shift
	exclude.min <- 10 # lowerbound of strand shift exclusion region
	exclude.max <- NaN # upperbound of strand shift exclusion region
	n.nodes <- NA # number of parallel processing nodes
	fdr <- 0.01 # false discovery rate threshold for peak calling
	npeak <- NA # threshold on number of peaks to call
	temp.dir <- tempdir() # temporary directory
	chrname.rm.pattern <- NA # chromosome name pattern used to remove tags
	output.odir <- NA # Output directory name
	output.npeak.file <- NA # Output narrowPeak file name
	output.rpeak.file <- NA # Output regionPeak file name
	output.rdata.file <- NA # Rdata file
	output.plot.file <- NA  # cross correlation plot file
	output.result.file <- NA # result file
	replace.flag <- FALSE # replace file flag
	clean.files.flag <- FALSE # file deletion flag
	
	# Parse arguments   
	for (each.arg in args) {
		
		if (grepl('^-c=',each.arg)) { #-c=<chip.file>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				chip.file <- arg.split[2] # second part is chip.file
			} else {
				stop('No tagAlign/BAM file name provided for parameter -c=')
			}
			
		} else if (grepl('^-i=',each.arg)) { #-i=<control.file>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				control.file <- arg.split[2] # second part is control.file
			} else {
				stop('No tagAlign/BAM file name provided for parameter -i=')
			}
			
		} else if (grepl('^-s=',each.arg)) { #-s=<sep.min>:<sep.bin>:<sep.max>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				sep.vals <- arg.split[2] # second part is sepmin:sepbin:sepmax
				sep.vals.split <- strsplit(sep.vals,':',fixed=TRUE)[[1]] # split on :                
				if (length(sep.vals.split) != 3) { # must have 3 parts
					stop('Strand shift limits must be specified as -s=sepmin:sepbin:sepmax')                    
				} else {
					if (any(is.na(as.numeric(sep.vals.split)))) { # check that sep vals are numeric
						stop('Strand shift limits must be numeric values')
					}
					sep.min <- round(as.numeric(sep.vals.split[1]))
					sep.bin <- round(as.numeric(sep.vals.split[2]))
					sep.max <- round(as.numeric(sep.vals.split[3]))
					if ((sep.min > sep.max) || (sep.bin > (sep.max - sep.min)) || (sep.bin < 0)) {
						stop('Illegal separation values -s=sepmin:sepbin:sepmax')
					}
				}                                    
			} else {
				stop('Strand shift limits must be specified as -s=sepmin:sepbin:sepmax')
			}
			
		} else if (grepl('^-speak=',each.arg)) { #-speak=<sep.peak> , user-defined cross-correlation peak strandshift
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				sep.peak <- arg.split[2] # second part is <sep.peak>
				if (is.na(as.numeric(sep.peak))) { # check that sep.peak is numeric
					stop('-speak=<sep.peak>: User defined peak shift must be numeric')
				}
				sep.peak <- as.numeric(sep.peak)
			} else {
				stop('User defined peak shift must be provided as -speak=<sep.peak>')
			}
			
		} else if (grepl('^-x=',each.arg)) { #-x=<exclude.min>:<exclude.max>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				exclude.vals <- arg.split[2] # second part is excludemin:excludemax
				exclude.vals.split <- strsplit(exclude.vals,':',fixed=TRUE)[[1]] # split on :                
				if (length(exclude.vals.split) != 2) { # must have 2 parts
					stop('Exclusion limits must be specified as -x=excludemin:excludemax')
				} else {
					if (any(is.na(as.numeric(exclude.vals.split)))) { # check that exclude vals are numeric
						stop('Exclusion limits must be numeric values')
					}
					exclude.min <- round(as.numeric(exclude.vals.split[1]))                    
					exclude.max <- round(as.numeric(exclude.vals.split[2]))
					if (exclude.min > exclude.max) {
						stop('Illegal exclusion limits -x=excludemin:excludemax')
					}                    
				}                                    
			} else {
				stop('Exclusion limits must be specified as -x=excludemin:excludemax')
			}
			
		} else if (grepl('^-p=',each.arg)) { #-p=<n.nodes> , number of parallel processing nodes, default=NULL            
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				n.nodes <- arg.split[2] # second part is numnodes
				if (is.na(as.numeric(n.nodes))) { # check that n.nodes is numeric
					stop('-p=<numnodes>: numnodes must be numeric')
				}
				n.nodes <- round(as.numeric(n.nodes))
			} else {
				stop('Number of parallel nodes must be provided as -p=<numnodes>')
			}
			
		} else if (grepl('^-fdr=',each.arg)) { #-fdr=<fdr> , false discovery rate, default=0.01            
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				fdr <- arg.split[2] # second part is fdr
				if (is.na(as.numeric(fdr))) { # check that fdr is numeric
					stop('-fdr=<falseDiscoveryRate>: false discovery rate must be numeric')
				}
				fdr <- as.numeric(fdr)
			} else {
				stop('False discovery rate must be provided as -fdr=<fdr>')
			}
			
		} else if (grepl('^-npeak=',each.arg)) { #-npeak=<numPeaks> , number of peaks threshold, default=NA            
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				npeak <- arg.split[2] # second part is npeak
				if (is.na(as.numeric(npeak))) { # check that npeak is numeric
					stop('-npeak=<numPeaks>: threshold on number of peaks must be numeric')
				}
				npeak <- round(as.numeric(npeak))
			} else {
				stop('Threshold on number of peaks must be provided as -npeak=<numPeaks>')
			}
			
		} else if (grepl('^-tmpdir=',each.arg)) { #-tmpdir=<temp.dir>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				temp.dir <- arg.split[2] # second part is temp.dir
			} else {
				stop('No temporary directory provided for parameter -tmpdir=')
			}

		} else if (grepl('^-filtchr=',each.arg)) { #-filtchr=<chrname.rm.pattern>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				chrname.rm.pattern <- arg.split[2] # second part is chrname.rm.pattern
			} else {
				stop('No pattern provided for parameter -filtchr=')
			}
			
		} else if (grepl('^-odir=',each.arg)) { #-odir=<output.odir>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				output.odir <- arg.split[2] # second part is output.odir
			} else {
				stop('No output directory provided for parameter -odir=')
			}
			
		} else if (grepl('^-savn',each.arg)) { # -savn=<output.npeak.file> OR -savn , save narrowpeak
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2])) {                
				output.npeak.file <- arg.split[2] #-savn=
			} else if (each.arg=='-savn') {
				output.npeak.file <- NULL # NULL indicates get the name from the main file name
			} else {
				stop('Argument for saving narrowPeak file must be -savn or -savn=<filename>')
			}
			
		} else if (grepl('^-savr',each.arg)) { # -savr=<output.rpeak.file> OR -savr , save regionpeak
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2])) {                
				output.rpeak.file <- arg.split[2] #-savr=
			} else if (each.arg=='-savr') {
				output.rpeak.file <- NULL # NULL indicates get the name from the main file name
			} else {
				stop('Argument for saving regionPeak file must be -savr or -savr=<filename>')
			}
			
		} else if (grepl('^-savd',each.arg)) { # -savd=<output.rdata.file> OR -savd , save Rdata file
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2])) {                
				output.rdata.file <- arg.split[2] #-savd=
			} else if (each.arg=='-savd') {
				output.rdata.file <- NULL # NULL indicates get the name from the main file name
			} else {
				stop('Argument for saving Rdata file must be -savd or -savd=<filename>')
			}
			
		} else if (grepl('^-savp',each.arg)) { # -savp=<output.plot.file> OR -savp , save cross-correlation plot                       
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2])) {                
				output.plot.file <- arg.split[2] #-savp=
			} else if (each.arg=='-savp') {
				output.plot.file <- NULL # NULL indicates get the name from the main file name
			} else {
				stop('Argument for saving Rdata file must be -savp or -savp=<filename>')
			}
			
		} else if (grepl('^-out=',each.arg)) { #-out=<output.result.file>
			
			arg.split <- strsplit(each.arg,'=',fixed=TRUE)[[1]] # split on =
			if (! is.na(arg.split[2]) ) {
				output.result.file <- arg.split[2] # second part is output.result.file
			} else {
				stop('No result file provided for parameter -out=')
			}
		
		} else if (each.arg == '-rf') {
			
			replace.flag <- TRUE
		
		} else if (each.arg == '-clean') {
			
			clean.files.flag <- TRUE
			
		} else {
			
			stop('Illegal argument ',each.arg)
		}        
	}
	# End: for loop
	
	# Check mandatory arguments
	if (is.na(chip.file)) {
		stop('-c=<tagAlign/BAMFileName> is a mandatory argument')
	}
	
	if (is.na(control.file) && ! is.na(output.npeak.file)) {
		stop('-i=<tagAlign/BAMFileName> is required for peak calling')
	}
	
	# Check if ChIP and control files are URLs
	if (grepl('^http://',chip.file)) {
		isurl.chip.file <- TRUE
	}
	if (grepl('^http://',control.file)) {
		isurl.control.file <- TRUE
	}
	
	# If ChIP file is a URL output.odir MUST be specified
	if (isurl.chip.file && is.na(output.odir)) {
		stop('If ChIP file is a URL, then output directory MUST be specified')
	}
	
	# Check that ChIP and control files exist
	if (isurl.chip.file) {
		if (system(paste('wget -q --spider',chip.file)) != 0) {
			stop('ChIP file URL not valid: ',chip.file)
		}
	} else if (!file.exists(chip.file)) {
		stop('ChIP File:',chip.file,' does not exist')
	}
	
	if (!is.na(control.file)) {
		if (isurl.control.file) {
			if (system(paste('wget -q --spider',control.file)) != 0) {
				stop('Control file URL not valid: ',control.file)
			}
		} else if (!file.exists(control.file)) {
			stop('Control File:',control.file,' does not exist')
		}   
	}
	
	# Correct other arguments
	if (is.na(output.odir)) { # Reconstruct output.odir if not provided
		output.odir <- get.file.parts(chip.file)$path
	}
	
	if (is.null(output.npeak.file)) { # Reconstruct output.npeak.file if NULL
		output.npeak.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '_VS_', get.file.parts(control.file)$name,'.narrowPeak', sep=""))
	}
	
	if (is.null(output.rpeak.file)) { # Reconstruct output.rpeak.file if NULL
		output.rpeak.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '_VS_', get.file.parts(control.file)$name,'.regionPeak', sep=""))
	}
	
	if (is.null(output.rdata.file)) { # Reconstruct output.rdata.file if NULL
		output.rdata.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '.Rdata', sep=""))
	}
	
	if (is.null(output.plot.file)) { # Reconstruct output.plot.file if NULL
		output.plot.file <- file.path(output.odir, paste(get.file.parts(chip.file)$name, '.pdf', sep=""))
	}
	
	return(list(chip.file=chip.file,
					isurl.chip.file=isurl.chip.file,
					control.file=control.file,
					isurl.control.file=isurl.control.file,
					sep.range=c(sep.min,sep.bin,sep.max),
					sep.peak=sep.peak,
					ex.range=c(exclude.min,exclude.max),
					n.nodes=n.nodes,
					fdr=fdr,
					npeak=npeak,
					temp.dir=temp.dir,
					chrname.rm.pattern=chrname.rm.pattern,
					output.odir=output.odir,
					output.npeak.file=output.npeak.file,
					output.rpeak.file=output.rpeak.file,
					output.rdata.file=output.rdata.file,
					output.plot.file=output.plot.file,
					output.result.file=output.result.file,
					replace.flag=replace.flag,
					clean.files.flag=clean.files.flag))          
} # end: parse.arguments()

read.align <- function(align.filename) {
# ===================================
# Function will read a tagAlign or BAM file
# ===================================	
	if (grepl('(\\.bam)?.*(\\.tagAlign)',align.filename)) { # if tagalign file
		chip.data <- read.tagalign.tags(align.filename)
		# get readlength info
		tmpDataRows <- read.table(align.filename,nrows=500)
		chip.data$read.length <- round(median(tmpDataRows$V3 - tmpDataRows$V2))
	} else if (grepl('(\\.tagAlign)?.*(\\.bam)',align.filename)) { # if bam file
		# create BAM file name
		bam2align.filename <- sub('\\.bam','.tagAlign',align.filename)
		# generate command to convert bam to tagalign
		command <- vector(length=2)
		command[1] <- sprintf("samtools view -F 0x0204 -o - %s",align.filename)
		command[2] <- paste("awk 'BEGIN{FS=" , '"\t"' , ";OFS=", '"\t"} {if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"}}', "' 1> ", bam2align.filename, sep="")
		# command[2] <- paste("awk 'BEGIN{OFS=", '"\t"} {if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"}}', "' 1> ", bam2align.filename, sep="")	
		command <- paste(command,collapse=" | ")
		# Run command
		status <- system(command,intern=FALSE,ignore.stderr=FALSE)
		if ((status != 0) || !file.exists(bam2align.filename)) {
			cat(sprintf("Error converting BAM to tagalign file: %s\n",align.filename),file=stderr())
			q(save="no",status=1)
		}
		# read converted BAM file
		chip.data <- read.tagalign.tags(bam2align.filename)
		# get readlength info
		tmpDataRows <- read.table(bam2align.filename,nrows=500)
		chip.data$read.length <- round(median(tmpDataRows$V3 - tmpDataRows$V2))
		# delete temporary tagalign file
		file.remove(bam2align.filename)	
	} else {
		cat(sprintf("Error:Unknown file format for file:%s\n",align.fname),file=stderr())
		q(save="no",status=1)	
	}
	return(chip.data)
} # end: read.align()

print.run.params <- function(params){
# ===================================
# Output run parameters
# ===================================		
	cat('################\n',file=stdout())
	cat(iparams$chip.file,
			iparams$control.file,
			iparams$sep.range,
			iparams$sep.peak,
			iparams$ex.range,
			iparams$n.nodes,
			iparams$fdr,
			iparams$npeak,
			iparams$output.odir,
			iparams$output.npeak.file,
			iparams$output.rpeak.file,
			iparams$output.rdata.file,
			iparams$output.plot.file,
			iparams$output.result.file,
			iparams$replace.flag,  
			labels=c('ChIP data:','Control data:', 'strandshift(min):','strandshift(step):','strandshift(max)','user-defined peak shift',
					'exclusion(min):','exclusion(max):','num parallel nodes:','FDR threshold:','NumPeaks Threshold:','Output Directory:',
					'narrowPeak output file name:', 'regionPeak output file name:', 'Rdata filename:', 
					'plot pdf filename:','result filename:','Overwrite files?:'),
			fill=18,
			file=stdout())	
	cat('\n',file=stdout())	
} # end: print.run.parameters()

check.replace.flag <- function(params){
# ===================================
# Check if files exist
# ===================================
# If replace.flag is NOT set, check if output files exist and abort if necessary
	if (! iparams$replace.flag) {
		if (! is.na(iparams$output.npeak.file)) {
			if (file.exists(iparams$output.npeak.file)) {
				cat('narrowPeak file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
				q(save="no",status=1)
			}
		}
		if (! is.na(iparams$output.rpeak.file)) {
			if (file.exists(iparams$output.rpeak.file)) {
				cat('regionPeak file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
				q(save="no",status=1)
			}
		}    
		if (! is.na(iparams$output.plot.file)) {
			if (file.exists(iparams$output.plot.file)) {
				cat('Plot file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
				q(save="no",status=1)
			}
		}
		if (! is.na(iparams$output.rdata.file)) {
			if (file.exists(iparams$output.rdata.file)) {
				cat('Rdata file already exists. Aborting Run. Use -rf if you want to overwrite\n',file=stderr())
				q(save="no",status=1)
			}
		}
	}	
}

# #############################################################################
# MAIN FUNCTION
# #############################################################################

# Check number of arguments
minargs = 1;
maxargs = 17;
if (nargs < minargs | nargs > maxargs) {
	print.usage()
	q(save="no",status=1)
}

# Parse arguments
# iparams$chip.file
# iparams$isurl.chip.file
# iparams$control.file
# iparams$isurl.control.file
# iparams$sep.range
# iparams$sep.peak
# iparams$ex.range
# iparams$n.nodes
# iparams$fdr
# iparams$npeak
# iparams$temp.dir
# iparams$output.odir
# iparams$output.npeak.file
# iparams$output.rpeak.file
# iparams$output.rdata.file
# iparams$output.plot.file
# iparams$output.result.file
# iparams$replace.flag
# iparams$clean.files.flag
iparams <- parse.arguments(args)

# Print run parameters
print.run.params(iparams)

# Check if output files exist 
check.replace.flag(iparams)

# curr.chip.file and curr.control.file always point to the original ChIP and control files on disk
# ta.chip.filename & ta.control.filename always point to the final but temporary versions of the ChIP and control files that will be passed to read.align

# Download ChIP and control files if necessary to temp.dir
if (iparams$isurl.chip.file) {
	curr.chip.file <- file.path(iparams$temp.dir, get.file.parts(iparams$chip.file)$fullname) # file is downloaded to temp.dir. Has same name as URL suffix
	cat('Downloading ChIP file:',iparams$chip.file,"\n",file=stdout())
	if (system(paste('wget -N -q -P',iparams$temp.dir,iparams$chip.file)) != 0) {
		stop('Error downloading ChIP file:',iparams$chip.file)
	}
} else {
	curr.chip.file <- iparams$chip.file # file is in original directory
}

if (iparams$isurl.control.file) {
	curr.control.file <- file.path(iparams$temp.dir, get.file.parts(iparams$control.file)$fullname) # file is downloaded to temp.dir. Has same name as URL suffix
	cat('Downloading control file:',iparams$control.file,"\n",file=stdout())
	if (system(paste('wget -N -q -P',iparams$temp.dir,iparams$control.file)) != 0) {
		stop('Error downloading Control file:',iparams$control.file)
	}
} else {
	curr.control.file <- iparams$control.file # file is in original directory
}

# unzip ChIP and input files if required AND copy to temp directory
if (get.file.parts(curr.chip.file)$ext == '.gz') {
	ta.chip.filename <- tempfile(get.file.parts(curr.chip.file)$name, tmpdir=iparams$temp.dir) # unzip file to temp.dir/[filename with .gz removed][randsuffix]
	cat('Decompressing ChIP file\n',file=stdout())
	if (system(paste("gunzip -c",curr.chip.file,">",ta.chip.filename)) != 0) {
		stop('Unable to decompress file:', iparams$chip.file)
	}
	if (iparams$clean.files.flag) { # Remove original file if clean.files.flag is set
		file.remove(curr.chip.file)
	}
} else {
	ta.chip.filename <- tempfile(get.file.parts(curr.chip.file)$fullname, tmpdir=iparams$temp.dir)
	if (iparams$clean.files.flag) {
		file.rename(curr.chip.file,ta.chip.filename) # move file to temp.dir/[filename][randsuffix]
	} else {
		file.copy(curr.chip.file,ta.chip.filename) # copy file to temp.dir/[filename][randsuffix]		
	}	
}

if (! is.na(iparams$control.file)) {
	if (get.file.parts(curr.control.file)$ext == '.gz') {
		ta.control.filename <- tempfile(get.file.parts(curr.control.file)$name, tmpdir=iparams$temp.dir) # unzip file to temp.dir/[filename with .gz removed][randsuffix]
		cat('Decompressing control file\n',file=stdout())        
		if (system(paste("gunzip -c",curr.control.file,">",ta.control.filename)) != 0) {
			stop('Unable to decompress file:', iparams$control.file)
		}
		if (iparams$clean.files.flag) { # Remove original file if clean.files.flag is set
			file.remove(curr.control.file)
		}				
	} else {
		ta.control.filename <- tempfile(get.file.parts(curr.control.file)$fullname, tmpdir=iparams$temp.dir) # copy file to temp.dir/[filename][randsuffix]
		
		if (iparams$clean.files.flag) {
			file.rename(curr.control.file,ta.control.filename) # move file to temp.dir/[filename][randsuffix]
		} else {
			file.copy(curr.control.file,ta.control.filename) # copy file to temp.dir/[filename][randsuffix]		
		}			
	}
}

# Remove downloaded files
if (iparams$isurl.chip.file & file.exists(curr.chip.file))  {
	file.remove(curr.chip.file)
}

if (! is.na(iparams$control.file)) {
	if (iparams$isurl.control.file & file.exists(curr.control.file)) {
		file.remove(curr.control.file)
	}
}

# Load SPP library
library(spp)

# Read ChIP tagAlign/BAM files
cat("Reading ChIP tagAlign/BAM file",iparams$chip.file,"\n",file=stdout())
chip.data <- read.align(ta.chip.filename)
cat("ChIP data read length",chip.data$read.length,"\n",file=stdout())
file.remove(ta.chip.filename) # Delete temporary file
if (length(chip.data$tags)==0) {
	stop('Error in ChIP file format:', iparams$chip.file)
}
# Remove illegal chromosome names
if (! is.na(iparams$chrname.rm.pattern)) {
	selectidx <- which(grepl(iparams$chrname.rm.pattern,names(chip.data$tags))==FALSE)
	chip.data$tags <- chip.data$tags[selectidx]
	chip.data$quality <- chip.data$quality[selectidx]
}
chip.data$num.tags <- sum(unlist(lapply(chip.data$tags,function(d) length(d))))

# Read Control tagAlign/BAM files
if (! is.na(iparams$control.file)) {
	cat("Reading Control tagAlign/BAM file",iparams$control.file,"\n",file=stdout())
	control.data <- read.align(ta.control.filename)
	file.remove(ta.control.filename) # Delete temporary file    
	if (length(control.data$tags)==0) {
		stop('Error in control file format:', iparams$chip.file)
	}    
	cat("Control data read length",control.data$read.length,"\n",file=stdout())
	# Remove illegal chromosome names
	if (! is.na(iparams$chrname.rm.pattern)) {
		selectidx <- which(grepl(iparams$chrname.rm.pattern,names(control.data$tags))==FALSE)
		control.data$tags <- control.data$tags[selectidx]
		control.data$quality <- control.data$quality[selectidx]
	}
	control.data$num.tags <- sum(unlist(lapply(control.data$tags,function(d) length(d))))
}

# Open multiple processes if required
if (is.na(iparams$n.nodes)) {
	cluster.nodes <- NULL
} else {
	library(snow)
	cluster.nodes <- makeCluster(iparams$n.nodes)
}

# #################################    
# Calculate cross-correlation for various strand shifts
# #################################    
cat("Calculating peak characteristics\n",file=stdout())
# crosscorr
# $cross.correlation : Cross-correlation profile as an $x/$y data.frame
# $peak : Position ($x) and height ($y) of automatically detected cross-correlation peak.
# $whs: Optimized window half-size for binding detection (based on the width of the cross-correlation peak) 
crosscorr <- get.binding.characteristics(chip.data,
		srange=iparams$sep.range[c(1,3)],
		bin=iparams$sep.range[2],
		accept.all.tags=T,
		cluster=cluster.nodes)
if (!is.na(iparams$n.nodes)) {
	stopCluster(cluster.nodes)
}

# Smooth the cross-correlation curve if required
cc <- crosscorr$cross.correlation
crosscorr$min.cc <- crosscorr$cross.correlation[ length(crosscorr$cross.correlation$y) , ] # minimum value and shift of cross-correlation
cat("Minimum cross-correlation value", crosscorr$min.cc$y,"\n",file=stdout())
cat("Minimum cross-correlation shift", crosscorr$min.cc$x,"\n",file=stdout())
sbw <- 2*floor(ceiling(5/iparams$sep.range[2]) / 2) + 1 # smoothing bandwidth
cc$y <- runmean(cc$y,sbw,alg="fast")

# Compute cross-correlation peak
bw <- ceiling(2/iparams$sep.range[2]) # crosscorr[i] is compared to crosscorr[i+/-bw] to find peaks
peakidx <- (diff(cc$y,bw)>=0) # cc[i] > cc[i-bw]
peakidx <- diff(peakidx,bw)
peakidx <- which(peakidx==-1) + bw        

# exclude peaks from the excluded region
if ( is.nan(iparams$ex.range[2]) ) {
	iparams$ex.range[2] <- chip.data$read.length+10
}
peakidx <- peakidx[(cc$x[peakidx] < iparams$ex.range[1]) | (cc$x[peakidx] > iparams$ex.range[2]) | (cc$x[peakidx] < 0) ]    
cc <- cc[peakidx,]

# Find max peak position and other peaks within 0.9*max_peakvalue that are further away from maxpeakposition   
maxpeakidx <- which.max(cc$y)
maxpeakshift <- cc$x[maxpeakidx]
maxpeakval <- cc$y[maxpeakidx]
peakidx <-which((cc$y >= 0.9*maxpeakval) & (cc$x >= maxpeakshift)) 
cc <- cc[peakidx,]

# sort the peaks and get the top 3
sortidx <- order(cc$y,decreasing=TRUE)
sortidx <- sortidx[c(1:min(3,length(sortidx)))]
cc.peak <- cc[sortidx,]

# Override peak shift if user supplies peak shift
if (! is.na(iparams$sep.peak)) {
	cc.peak <- approx(crosscorr$cross.correlation$x,crosscorr$cross.correlation$y,iparams$sep.peak,rule=2)
}
cat("Top 3 cross-correlation values", paste(cc.peak$y,collapse=","),"\n",file=stdout())
cat("Top 3 estimates for fragment length",paste(cc.peak$x,collapse=","),"\n",file=stdout())

# Reset values in crosscorr
crosscorr$peak$x <- cc.peak$x[1]
crosscorr$peak$y <- cc.peak$y[1]

# Compute window half size
whs.thresh <- crosscorr$min.cc$y + (crosscorr$peak$y - crosscorr$min.cc$y)/3
crosscorr$whs <- max(crosscorr$cross.correlation$x[crosscorr$cross.correlation$y >= whs.thresh])
cat("Window half size",crosscorr$whs,"\n",file=stdout())

# Compute phantom peak coefficient
ph.peakidx <- which( ( crosscorr$cross.correlation$x >= ( chip.data$read.length - round(2*iparams$sep.range[2]) ) ) & 
                     ( crosscorr$cross.correlation$x <= ( chip.data$read.length + round(2*iparams$sep.range[2]) ) ) )
ph.peakidx <- ph.peakidx[ which.max(crosscorr$cross.correlation$y[ph.peakidx]) ]
crosscorr$phantom.cc <- crosscorr$cross.correlation[ph.peakidx,]
cat("Phantom peak location",crosscorr$phantom.cc$x,"\n",file=stdout())
cat("Phantom peak Correlation",crosscorr$phantom.cc$y,"\n",file=stdout())
crosscorr$phantom.coeff <- crosscorr$peak$y / crosscorr$phantom.cc$y
crosscorr$phantom.coeff <- crosscorr$peak$y / crosscorr$min.cc$y
cat("Normalized Strand cross-correlation coefficient (NSC)",crosscorr$phantom.coeff,"\n",file=stdout())
crosscorr$rel.phantom.coeff <- (crosscorr$peak$y - crosscorr$min.cc$y) / (crosscorr$phantom.cc$y - crosscorr$min.cc$y)
cat("Relative Strand cross-correlation Coefficient (RSC)",crosscorr$rel.phantom.coeff,"\n",file=stdout())
crosscorr$phantom.quality.tag <- NA
if ( (crosscorr$rel.phantom.coeff >= 0) & (crosscorr$rel.phantom.coeff < 0.25) ) {
	crosscorr$phantom.quality.tag <- -2
} else if ( (crosscorr$rel.phantom.coeff >= 0.25) & (crosscorr$rel.phantom.coeff < 0.5) ) {
	crosscorr$phantom.quality.tag <- -1
} else if ( (crosscorr$rel.phantom.coeff >= 0.5) & (crosscorr$rel.phantom.coeff < 1) ) {
	crosscorr$phantom.quality.tag <- 0
} else if ( (crosscorr$rel.phantom.coeff >= 1) & (crosscorr$rel.phantom.coeff < 1.5) ) {
	crosscorr$phantom.quality.tag <- 1
} else if ( (crosscorr$rel.phantom.coeff >= 1.5) ) {
	crosscorr$phantom.quality.tag <- 2
}
cat("Phantom Peak Quality Tag",crosscorr$phantom.quality.tag,"\n",file=stdout())

# Output result to result file if required
#Filename\tnumReads\tPeak_shift\tPeak_Correlation\tRead_length\tPhantomPeak_Correlation\tMin_Correlation_Shift\tMin_Correlation\tNormalized_CrossCorrelation_Coefficient\tRelative_CrossCorrelation_Coefficient\tQualityTag)
if (! is.na(iparams$output.result.file)) {
	cat(get.file.parts(iparams$chip.file)$fullname,
			chip.data$num.tags,
			paste(cc.peak$x,collapse=","),
			paste(cc.peak$y,collapse=","),
			crosscorr$phantom.cc$x,
			crosscorr$phantom.cc$y,
			crosscorr$min.cc$x,
			crosscorr$min.cc$y,
			crosscorr$phantom.coeff,
			crosscorr$rel.phantom.coeff,
			crosscorr$phantom.quality.tag,
			sep="\t",
			file=iparams$output.result.file,
			append=TRUE)
	cat("\n",
			file=iparams$output.result.file,
			append=TRUE)    
}

# Save figure if required
if (! is.na(iparams$output.plot.file)) {
	pdf(file=iparams$output.plot.file,width=5,height=5)
	par(mar = c(4,3.5,2,0.5), mgp = c(1.5,0.5,0), cex = 0.8);
	plot(crosscorr$cross.correlation,
			type='l',
			xlab=sprintf("strand-shift (%s)",paste(cc.peak$x,collapse=",")),
			ylab="cross-correlation")
	abline(v=cc.peak$x,lty=2,col=2)
	abline(v=crosscorr$phantom.cc$x,lty=2,col=4)
	title(main=get.file.parts(iparams$chip.file)$fullname,
	      sub=sprintf("NSC=%g,RSC=%g,Qtag=%d",crosscorr$phantom.coeff,crosscorr$rel.phantom.coeff,crosscorr$phantom.quality.tag))
	dev.off();    
}

# Save RData file if required
if (! is.na(iparams$output.rdata.file)) {
	save(iparams,
			crosscorr,
			cc.peak,
			file=iparams$output.rdata.file);    
}

# #################################    
# Call peaks
# #################################

if ( !is.na(iparams$output.npeak.file) || !is.na(iparams$output.rpeak.file) ) {
	
	# Remove local tag anomalies
	cat('Removing read stacks\n',file=stdout())
	chip.data <- remove.local.tag.anomalies(chip.data$tags)
	control.data <- remove.local.tag.anomalies(control.data$tags)
	
	# Open multiple processes if required
	if (is.na(iparams$n.nodes)) {
		cluster.nodes <- NULL
	} else {
		cluster.nodes <- makeCluster(iparams$n.nodes)
	}
	
	# Find peaks
	cat('Finding peaks\n',file=stdout())
	if (!is.na(iparams$npeak)) {
		iparams$fdr <- 0.99
	}
	narrow.peaks <- find.binding.positions(signal.data=chip.data,control.data=control.data,fdr=iparams$fdr,method=tag.lwcc,whs=crosscorr$whs,cluster=cluster.nodes)
	if (!is.na(iparams$n.nodes)) {
		stopCluster(cluster.nodes)
	}
	cat(paste("Detected",sum(unlist(lapply(narrow.peaks$npl,function(d) length(d$x)))),"peaks"),"\n",file=stdout())
	
	# Write to narrowPeak file
	if (!is.na(iparams$output.npeak.file)) {
		write.narrowpeak.binding(narrow.peaks,iparams$output.npeak.file,margin=round(crosscorr$whs/2),npeaks=iparams$npeak)
		system(paste('gzip -f ',iparams$output.npeak.file))
	}
	
	# Compute and write regionPeak file
	if (!is.na(iparams$output.rpeak.file)) {
		region.peaks <- add.broad.peak.regions(chip.data,control.data,narrow.peaks,window.size=max(50,round(crosscorr$whs/4)),z.thr=10)
		write.narrowpeak.binding(region.peaks,iparams$output.rpeak.file,margin=round(crosscorr$whs/2),npeaks=iparams$npeak)
		system(paste('gzip -f ',iparams$output.rpeak.file))
	}
	
	# Save Rdata file    
	if (! is.na(iparams$output.rdata.file)) {
		save(iparams,
				crosscorr,
				cc.peak,
				narrow.peaks,
				region.peaks,
				file=iparams$output.rdata.file);    
	}
	
}


