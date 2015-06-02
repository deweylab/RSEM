===========================
README for consistency analysis of peak calling on replicates
Qunhua Li and Anshul Kundaje (Oct,2010)
===========================
This set of programs are used for consistency analysis on peak calling results on multiple replicates of a dataset

================
DEPENDENCIES
================
unix, R version 2.9 or higher

================
FILES:
================
batch-consistency-analysis.r : for pairwise IDR analysis of replicates
batch-consistency-plot.r: for creating diagnostic and IDR plots
functions-all-clayton-12-13.r: helper function
genome_table.txt: This file MUST contain the size of each chromosome of the genome of the organism that the peak files are referring to

================
INPUT FILE FORMATS
================
(1) genome_table.txt
It contains two space delimited fields
Col1: chromosome name (These MUST match the chromosome names in the peak files)
Col2: chromosome size (in bp)

(1) Peak Files
Peak files MUST be in narrowPeak format (and unzipped ... the code currently doesnt handle gzipped peak files directly)

NarrowPeak files are in BED6+4 format. It consists of 10 tab-delimited columns

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

*NOTE*: the p-value and q-value columns MUST be in -log10() scale

The narrowPeak format has 3 columns that can be used to rank peaks 
(signal.value, p.value (-log_10) and q.value (-log_10)).
The peak summit column must have values relative to the start coordinate of the peaks.
You can use any of these columns but make sure that whichever measure you are using rank peaks is relatively continuous without too many ties.
e.g. For the SPP peak caller it is recommended to use signal.value column
e.g. PeakSeq peak caller has relatively continuous q.values without too many ties. So for PeakSeq it is better to use q.value
 
================
RUNNING INSTRUCTIONS
================
First make sure the genome_table.txt file contains the appropriate chromosome names and sizes. If not replace the contents of this file. Make sure the file continues to be named 'genome_table.txt'.
The file name is currently hardcoded. We will change this in the next release of the code.

(1) batch-consistency-analysis.r

This is used to run pairwise consistency analysis on a pair of replicate peak files

----------------
GENERAL USAGE:
----------------
Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] [peak.half.width] [outfile.prefix] [min.overlap.ratio] [is.broadpeak] [ranking.measure]

Typical usage for SPP peak caller peaks
Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] -1 [outfile.prefix] 0 F q.value

Typical usage for MACS peak caller peaks
Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] 200 [outfile.prefix] 0 F p.value

[peakfile1] and [peakfile2] are the peak calls for the pair of replicates in narrowPeak format. They must be uncompressed files.
e.g. /peaks/reps/chipSampleRep1_VS_controlSampleRep0.narrowPeak AND 
	 /peaks/reps/chipSampleRep2_VS_controlSampleRep0.narrowPeak

[peak.half.width]: Set this to -1 if you want to use the reported peak width in the peak files.
If you want to truncate peak widths to say 400 bp max then use a value of 200.

[outfile.prefix] is a prefix that will be used to name the output data for this pair of replicates.
The prefix must also include the PATH to the directory where you want to store the output data.
e.g. /consistency/reps/chipSampleRep1_VS_chipSampleRep2

[min.overlap.ratio]: fractional bp overlap (ranges from 0 to 1) between peaks in replicates to be considered as overlapping peaks. 
Set to 0 if you want to allow overlap to be defined as >= 1 bp overlap.
If set to say 0.5 this would mean that atleast 50% of the peak in one replicate should be covered by a peak in the other replicate to count as an overlap.

[is.broadpeak]: Is the peak file format narrowPeak or broadPeak. Set to F if it is narrowPeak/regionPeak or T if it is broadPeak.

[ranking.measure] is the ranking measure to use. It can take only one of the following values
signal.value , p.value or q.value

OUTPUT:
The results will be written to the directory contained in [outfile.prefix]
a. The output from EM fitting: suffixed by -em.sav
b. The output for plotting empirical curves: suffixed by -uri.sav
	Note: 1 and 2 are objects that can be loaded back to R for plotting or other purposes (e.g. retrieve data)
c. The parameters estimated from EM and the log of consistency analysis, suffixed by -Rout.txt
d. The number of peaks that pass specific IDR thresholds for the pairwise analysis: suffixed by npeaks-aboveIDR.txt
e. The full set of peaks that overlap between the replicates with local and global IDR scores: suffixed by overlapped-peaks.txt


(2) batch-consistency-plot.r

This is used to plot the IDR plots and diagnostic plots for a single or multiple pairs of replicates.

----------------
GENERAL USAGE:
----------------
Rscript batch-consistency-plot.r [npairs] [output.prefix] [input.file.prefix1] [input.file.prefix2] [input.file.prefix3] ....

[n.pairs] is the number of pairs of replicates that you want to plot on the same plot
e.g. 1 or 3 or ...

[output.prefix] is a prefix that will be used to name output data from this analysis. 
NOT TO BE CONFUSED with [outfile.prefix] in batch-consistency-analysis.r
The prefix must also include the PATH to the directory where you want to store the output data.
e.g. /consistency/plots/chipSampleAllReps

[input.file.prefix 1, 2, 3 ...] are the [outfile.prefix] values used to name the output from pairwise analysis on all replicates
e.g. /consistency/reps/chipSampleRep1_VS_chipSampleRep2
	 /consistency/reps/chipSampleRep1_VS_chipSampleRep3
	 /consistency/reps/chipSampleRep2_VS_chipSampleRep3

OUTPUT:
1. summary consistency plots in .ps format: suffixed by -plot.ps
These plots are very informative about the quality and similarity of the replicates.

===================================================
GETTING NUMBER OF PEAKS THAT PASS AN IDR THRESHOLD
===================================================
For each pairwise analysis, we have a *overlapped-peaks.txt file

The last column (Column 11) of the overlapped-peaks.txt file has the global IDR score for each pair of overlapping peaks
To get the number of peaks that pass an IDR threshold of T (e.g. 0.01) you simply find the number of lines that have a global IDR score <= T

awk '$11 <= 0.01 {print $0}' [overlappedPeaksFileName] | wc -l



