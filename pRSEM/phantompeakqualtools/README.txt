===========================
Anshul Kundaje
Date: Feb 13 2012
Email: anshul@kundaje.net
Version: 2.0
===========================
This set of programs operate on mapped Illumina single-end read datasets in tagAlign or BAM format.
They can be used to 
(1) Compute the predominant fragment length based on strand cross-correlation peak
(2) Compute Data quality measures based on strand cross-correlation analysis and relative phantom peak
(3) Call Peaks and regions for punctate binding ChIP-seq datasets

===========================
CITATIONS:
===========================
If you are using the code or results in any formal publication please cite
[1] Anshul Kundaje, Computer Science Dept., MIT, ENCODE Consortium, http://code.google.com/p/phantompeakqualtools, Feb 2013
[2] Kharchenko PK, Tolstorukov MY, Park PJ, Design and analysis of ChIP-seq experiments for DNA-binding proteins Nat Biotechnol. 2008 Dec;26(12):1351-9

===========================
DEPENDENCIES:
===========================
unix,bash,R-2.10 and above,awk,samtools,boost C++ libraries
R packages: SPP, caTools, snow
NOTE: The current package does not run on a MAC or WINDOWS.

===========================
FILES:
===========================
(1) spp_1.10.1.tar.gz  : modified SPP peak-caller package (The original SPP-peak caller package was written by Peter Kharchenko[2])
(2) run_spp.R          : The script to compute the frag length, data quality characteristics based on cross-correlation analysis and/or peak calling
(3) run_spp_nodups.R   : (FOR DATASETS WHERE DUPLICATES ARE REMOVED i.e. MAX 1 READ STARTING AT ANY GENOMIC LOCATION) The script to compute the frag length, data quality characteristics based on cross-correlation analysis and/or peak calling
(4) README.txt         : This README

============================
INSTALLATION:
============================
(1) First make sure that you have installed R (version 2.10 or higher)

(2) Also, you must have the Boost C++ libraries installed. Most linux distributions have these preinstalled.
If not, you can easily get these from your standard package manager for your linux distribution.
e.g synaptic package manager (apt-get) for ubuntu or emerge for gentoo.

(3) Install the following R packages
	- caTools
	- snow (if you want parallel processing)
from within R
install.packages([packageName],dependencies=TRUE)

(4) You can then install the SPP package spp_1.10.X.tar.gz
<From your bash shell>
	R CMD INSTALL spp_1.10.X.tar.gz

<From within R>
	install.packages('spp_1.10.X.tar.gz',dependencies=TRUE)

(5) If your alignment files are BAM, you must have the samtools executable in your path so that the R script run_spp.R can call it using the system() command
You can get samtools from here http://samtools.sourceforge.net/
You can add the following line to your .bashrc file
	export PATH="<path_to_samtools_executable>:${PATH}"
	
(6) Run run_spp.R
	Rscript run_spp.R <options>

===========================
GENERAL USAGE
===========================
Usage: Rscript run_spp.R <options>

MANDATORY ARGUMENTS
-c=<ChIP_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)

MANDATORY ARGUMENTS FOR PEAK CALLING
-i=<Input_alignFile>, full path and name (or URL) of tagAlign/BAM file (can be gzipped) (FILE EXTENSION MUST BE tagAlign.gz, tagAlign, bam or bam.gz)

OPTIONAL ARGUMENTS
-s=<min>:<step>:<max> , strand shifts at which cross-correlation is evaluated, default=-100:5:600
-speak=<strPeak>, user-defined cross-correlation peak strandshift
-x=<min>:<max>, strand shifts to exclude (This is mainly to avoid region around phantom peak) default=10:(readlen+10)
-p=<nodes> , number of parallel processing nodes, default=0
-fdr=<falseDisoveryRate> , false discovery rate threshold for peak calling
-npeak=<numPeaks>, threshold on number of peaks to call
-tmpdir=<tempdir> , Temporary directory (if not specified R function tempdir() is used)
-filtchr=<chrnamePattern> , Pattern to use to remove tags that map to specific chromosomes e.g. _ will remove all tags that map to chromosomes with _ in their name

OUTPUT ARGUMENTS
-odir=<outputDirectory> name of output directory (If not set same as ChIP file directory is used)
-savn=<narrowpeakfilename> OR -savn NarrowPeak file name (fixed width peaks)
-savr=<regionpeakfilename> OR -savr RegionPeak file name (variable width peaks with regions of enrichment around peak summits)
-savd=<rdatafile> OR -savd, save Rdata file
-savp=<plotdatafile> OR -savp, save cross-correlation plot
-out=<resultfile>, append peakshift/phantomPeak results to a file
-rf, if plot or rdata or narrowPeak file exists replace it. If not used then the run is aborted if the plot or Rdata or narrowPeak file exists
-clean, if used it will remove the original chip and control files after reading them in. CAUTION: Use only if the script calling run_spp.R is creating temporary files

===========================
TYPICAL USAGE
===========================
(1) Determine strand cross-correlation peak / predominant fragment length OR print out quality measures
	
	Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>

-savp will create a pdf showing the cross-correlation plot
-out=<outFile> will create and/or append to a file named <outFile> several important characteristics of the dataset.
The file contains 11 tab delimited columns

COL1: Filename: tagAlign/BAM filename
COL2: numReads: effective sequencing depth i.e. total number of mapped reads in input file
COL3: estFragLen: comma separated strand cross-correlation peak(s) in decreasing order of correlation.
	  The top 3 local maxima locations that are within 90% of the maximum cross-correlation value are output.
	  In almost all cases, the top (first) value in the list represents the predominant fragment length.
	  If you want to keep only the top value simply run
	  sed -r 's/,[^\t]+//g' <outFile> > <newOutFile>
COL4: corr_estFragLen: comma separated strand cross-correlation value(s) in decreasing order (col2 follows the same order)
COL5: phantomPeak: Read length/phantom peak strand shift
COL6: corr_phantomPeak: Correlation value at phantom peak
COL7: argmin_corr: strand shift at which cross-correlation is lowest
COL8: min_corr: minimum value of cross-correlation
COL9: Normalized strand cross-correlation coefficient (NSC) = COL4 / COL8
COL10: Relative strand cross-correlation coefficient (RSC) = (COL4 - COL8) / (COL6 - COL8)
COL11: QualityTag: Quality tag based on thresholded RSC (codes: -2:veryLow,-1:Low,0:Medium,1:High,2:veryHigh)

You can run the program on multiple datasets in parallel and append all the quality information to the same <outFile> for a summary analysis.

NSC values range from a minimum of 1 to larger positive numbers. 1.1 is the critical threshold. 
Datasets with NSC values much less than 1.1 (< 1.05) tend to have low signal to noise or few peaks (this could be biological eg.a factor that truly binds only a few sites in a particular tissue type OR it could be due to poor quality)

RSC values range from 0 to larger positive values. 1 is the critical threshold.
RSC values significantly lower than 1 (< 0.8) tend to have low signal to noise. The low scores can be due to failed and poor quality ChIP, low read sequence quality and hence lots of mismappings, shallow sequencing depth (significantly below saturation) or a combination of these. Like the NSC, datasets with few binding sites (< 200) which is biologically justifiable also show low RSC scores.

Qtag is a thresholded version of RSC.

(2) Peak calling

Rscript run_spp.R -c=<ChIP_tagalign/BAM_file> -i=<control_tagalign/BAM_file> -fdr=<fdr> -odir=<peak_call_output_dir> -savr -savp -savd -rf
Rscript run_spp.R -c=<ChIP_tagalign/BAM_file> -i=<control_tagalign/BAM_file> -npeak=<npeaks> -odir=<peak_call_output_dir> -savr -savp -savd -rf

(3) For IDR analysis you want to call a large number of peaks (relaxed threshold) so that the IDR model has access to a sufficient noise component.

Rscript run_spp.R -c=<ChIP_tagalign/BAM_file> -i=<control_tagalign/BAM_file> -npeak=300000 -odir=<peak_call_output_dir> -savr -savp -rf -out=<resultFile>

===========================
NOTES:
===========================
- It is EXTREMELY important to filter out multi-mapping reads from the BAM/tagAlign files. Large number of multimapping reads can severly affect the phantom peak coefficient and peak calling results.

- If a dataset seems to have high PCR bottlenecking, then you might want to actually clamp the number of unique mappping reads per position to 1 or upto 5. If not the phantom peak coefficient can be artificially good.

- For the IDR rescue strategy, one needs to pool reads from replicates and then shuffle and subsample the mapped reads to create two balanced pseudoReplicates. This is much easier to implement on tagAlign/BED read-mapping files using the unix 'shuf' command. So it is recommended to use the tagAlign format.

- In most cases, you can simply use the maximum reported strand correlation peak as the predominant fragment length.
However, it is useful to manually take a look at the cross-correlation plot to make sure the selected max peak is not an artifact.

- Also, if there are problems with library size-selection, a dataset's cross-correlation profile can have multiple strong cross-correlation peaks. This is currently not autodetected.

===========================
INPUT FILE FORMATS:
===========================
(1) BAM format
This is a binary alignment format specified in http://samtools.sourceforge.net/SAM-1.3.pdf
You MUST have samtools installed to use run_spp.R with BAM files

(2) TagAlign files
This a text-based BED3+3 alignment format that is easier to manipulate. It contains 6 tab delimited columns.

chrom	 string	 Name of the chromosome
chromStart	 int	 The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
chromEnd	 int	 The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as
                     chromStart=0, chromEnd=100, and span the bases numbered 0-99.
sequence	 string	 Sequence of this read
score	 int	 Indicates uniqueness or quality (preferably 1000/alignmentCount).
strand	 char	 Orientation of this read (+ or -)

NOTE: You dont have to store the sequence of reads in the sequence field as the peak caller never really uses that field. You can just put the letter 'N' in that field. This saves space significantly.

For the IDR rescue strategy, one needs to use shuffled and subsampled version of the alignment files. This is much easier to implement on tagAlign text files using the unix 'shuf' command.
So it is recommended to preferably use the tagAlign format.

----------------------------------
CONVERTING BAM TO TAGALIGN FILES
----------------------------------
It is very quick to convert BAM files to gzipped tagAlign files using

samtools view -F 0x0204 -o - <bamFile> | awk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' | gzip -c > <gzip_TagAlignFileName>

===========================
OUTPUT FILE FORMATS:
===========================
(1) NarrowPeak/RegionPeak format

The output peak file is in BED6+4 format known as tagAlign. It consists of 10 tab-delimited columns

chrom	 string	 Name of the chromosome
chromStart	 int	 The starting position of the feature in the chromosome. The first base in a chromosome is numbered 0.
chromEnd	 int	 The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature.
                     For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
name	 string	 Name given to a region (preferably unique). Use '.' if no name is assigned.
score	 int	 Indicates how dark the peak will be displayed in the browser (1-1000). If '0', the DCC will assign this based on signal value. Ideally average signalValue per base spread between 100-1000.
strand	 char	 +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
signalValue	 float	 Measurement of overall (usually, average) enrichment for the region.
pValue	 float	 Measurement of statistical signficance (-log10). Use -1 if no pValue is assigned.
qValue	 float	 Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
peak	 int	 Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
