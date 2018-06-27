README for RSEM
===============

[Bo Li](http://bli25ucb.github.io/) \(bli at cs dot wisc dot edu\)

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Compilation & Installation](#compilation)
* [Usage](#usage)
    * [Build RSEM references using RefSeq, Ensembl, or GENCODE annotations](#built)
    * [Build RSEM references for untypical organisms](#untypical)
* [Example](#example-main)
* [Simulation](#simulation)
* [Generate Transcript-to-Gene-Map from Trinity Output](#gen_trinity)
* [Differential Expression Analysis](#de)
* [Prior-Enhanced RSEM (pRSEM)](#pRSEM)
* [Authors](#authors)
* [Acknowledgements](#acknowledgements)
* [License](#license)

* * *

## <a name="introduction"></a> Introduction

RSEM is a software package for estimating gene and isoform expression
levels from RNA-Seq data. The RSEM package provides an user-friendly
interface, supports threads for parallel computation of the EM
algorithm, single-end and paired-end read data, quality scores,
variable-length reads and RSPD estimation. In addition, it provides
posterior mean and 95% credibility interval estimates for expression
levels. For visualization, It can generate BAM and Wiggle files in
both transcript-coordinate and genomic-coordinate. Genomic-coordinate
files can be visualized by both UCSC Genome browser and Broad
Institute's Integrative Genomics Viewer (IGV). Transcript-coordinate
files can be visualized by IGV. RSEM also has its own scripts to
generate transcript read depth plots in pdf format. The unique feature
of RSEM is, the read depth plots can be stacked, with read depth
contributed to unique reads shown in black and contributed to
multi-reads shown in red. In addition, models learned from data can
also be visualized. Last but not least, RSEM contains a simulator.

## <a name="compilation"></a> Compilation & Installation

To compile RSEM, simply run
   
    make

For Cygwin users, run

    make cygwin=true

To compile EBSeq, which is included in the RSEM package, run

    make ebseq

To install RSEM, simply put the RSEM directory in your environment's PATH
variable. Alternatively, run

    make install

By default, RSEM executables are installed to `/usr/local/bin`. You
can change the installation location by setting `DESTDIR` and/or
`prefix` variables. The RSEM executables will be installed to
`${DESTDIR}${prefix}/bin`. The default values of `DESTDIR` and
`prefix` are `DESTDIR=` and `prefix=/usr/local`. For example,

    make install DESTDIR=/home/my_name prefix=/software

will install RSEM executables to `/home/my_name/software/bin`.

**Note** that `make install` does not install `EBSeq` related scripts,
such as `rsem-generate-ngvector`, `rsem-run-ebseq`, and
`rsem-control-fdr`. But `rsem-generate-data-matrix`, which generates
count matrix for differential expression analysis, is installed.

### Prerequisites

C++, Perl and R are required to be installed. 

To use the `--gff3` option of `rsem-prepare-reference`, Python is also
required to be installed.

To take advantage of RSEM's built-in support for the Bowtie/Bowtie
2/STAR alignment program, you must have
[Bowtie](http://bowtie-bio.sourceforge.net)/[Bowtie
2](http://bowtie-bio.sourceforge.net/bowtie2)/[STAR](https://github.com/alexdobin/STAR)
installed.

## <a name="usage"></a> Usage

### I. Preparing Reference Sequences

RSEM can extract reference transcripts from a genome if you provide it
with gene annotations in a GTF/GFF3 file.  Alternatively, you can provide
RSEM with transcript sequences directly.

Please note that GTF files generated from the UCSC Table Browser do not
contain isoform-gene relationship information.  However, if you use the
UCSC Genes annotation track, this information can be recovered by
downloading the knownIsoforms.txt file for the appropriate genome.
 
To prepare the reference sequences, you should run the
`rsem-prepare-reference` program.  Run 

    rsem-prepare-reference --help

to get usage information or visit the [rsem-prepare-reference
documentation page](rsem-prepare-reference.html).

#### <a name="built"></a> Build RSEM references using RefSeq, Ensembl, or GENCODE annotations

RefSeq and Ensembl are two frequently used annotations. For human and
mouse, GENCODE annotaions are also available. In this section, we show
how to build RSEM references using these annotations. Note that it is
important to pair the genome with the annotation file for each
annotation source. In addition, we recommend users to use the primary
assemblies of genomes. Without loss of generality, we use human genome as
an example and in addition build Bowtie indices. 

For **RefSeq**, the genome and annotation file in GFF3 format can be found
at RefSeq genomes FTP:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```

For example, the human genome and GFF3 file locate at the subdirectory
`vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.31_GRCh38.p5`. `GCF_000001405.31_GRCh38.p5`
is the latest annotation version when this section was written.

Download and decompress the genome and annotation files to your working directory:

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_genomic.fna.gz
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.31_GRCh38.p5/GCF_000001405.31_GRCh38.p5_genomic.gff.gz
```

`GCF_000001405.31_GRCh38.p5_genomic.fna` contains all top level
sequences, including patches and haplotypes. To obtain the primary
assembly, run the following RSEM python script:

```
rsem-refseq-extract-primary-assembly GCF_000001405.31_GRCh38.p5_genomic.fna GCF_000001405.31_GRCh38.p5_genomic.primary_assembly.fna
```

Then type the following command to build RSEM references:

```
rsem-prepare-reference --gff3 GCF_000001405.31_GRCh38.p5_genomic.gff \
		       --trusted-sources BestRefSeq,Curated\ Genomic \
		       --bowtie \
		       GCF_000001405.31_GRCh38.p5_genomic.primary_assembly.fna \
		       ref/human_refseq
```

In the above command, `--trusted-sources` tells RSEM to only extract
transcripts from RefSeq sources like `BestRefSeq` or `Curated Genomic`. By
default, RSEM trust all sources. There is also an
`--gff3-RNA-patterns` option and its default is `mRNA`. Setting
`--gff3-RNA-patterns mRNA,rRNA` will allow RSEM to extract all mRNAs
and rRNAs from the genome. Visit [here](rsem-prepare-reference.html)
for more details.

Because the gene and transcript IDs (e.g. gene1000, rna28655)
extracted from RefSeq GFF3 files are hard to understand, it is
recommended to turn on the `--append-names` option in
`rsem-calculate-expression` for better interpretation of
quantification results.

For **Ensembl**, the genome and annotation files can be found at
[Ensembl FTP](http://uswest.ensembl.org/info/data/ftp/index.html).

Download and decompress the human genome and GTF files:

```
ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.83.gtf.gz
```

Then use the following command to build RSEM references:

```
rsem-prepare-reference --gtf Homo_sapiens.GRCh38.83.gtf \
		       --bowtie \
		       Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		       ref/human_ensembl
```

If you want to use GFF3 file instead, which is unnecessary and not
recommended, you should add option `--gff3-RNA-patterns transcript`
because `mRNA` is replaced by `transcript` in Ensembl GFF3 files.

**GENCODE** only provides human and mouse annotations. The genome and
  annotation files can be found from [GENCODE
  website](http://www.gencodegenes.org/).

Download and decompress the human genome and GTF files:

```
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
```

Then type the following command:

```
rsem-prepare-reference --gtf gencode.v24.annotation.gtf \
		       --bowtie \
		       GRCh38.primary_assembly.genome.fa \
		       ref/human_gencode
```

Similar to Ensembl annotation, if you want to use GFF3 files (not
recommended), add option `--gff3-RNA-patterns transcript`.

#### <a name="untypical"></a> Build RSEM references for untypical organisms

For untypical organisms, such as viruses, you may only have a GFF3 file that containing only genes but not any transcripts. You need to turn on `--gff3-genes-as-transcripts` so that RSEM will make each gene as a unique transcript.

Here is an example command:

```
rsem-prepare-reference --gff3 virus.gff \
               --gff3-genes-as-transcripts \
               --bowtie \
               virus.genome.fa \
               ref/virus
```

### II. Calculating Expression Values

To calculate expression values, you should run the
`rsem-calculate-expression` program.  Run 

    rsem-calculate-expression --help

to get usage information or visit the [rsem-calculate-expression
documentation page](rsem-calculate-expression.html).

#### Calculating expression values from single-end data

For single-end models, users have the option of providing a fragment
length distribution via the `--fragment-length-mean` and
`--fragment-length-sd` options.  The specification of an accurate fragment
length distribution is important for the accuracy of expression level
estimates from single-end data.  If the fragment length mean and sd are
not provided, RSEM will not take a fragment length distribution into
consideration.

#### Using an alternative aligner

By default, RSEM automates the alignment of reads to reference
transcripts using the Bowtie aligner. Turn on `--bowtie2` for
`rsem-prepare-reference` and `rsem-calculate-expression` will allow
RSEM to use the Bowtie 2 alignment program instead. Please note that
indel alignments, local alignments and discordant alignments are
disallowed when RSEM uses Bowtie 2 since RSEM currently cannot handle
them. See the description of `--bowtie2` option in
`rsem-calculate-expression` for more details. Similarly, turn on
`--star` will allow RSEM to use the STAR aligner. To use an
alternative alignment program, align the input reads against the file
`reference_name.idx.fa` generated by `rsem-prepare-reference`, and
format the alignment output in SAM/BAM/CRAM format.  Then, instead of
providing reads to `rsem-calculate-expression`, specify the
`--alignments` option and provide the SAM/BAM/CRAM file as an
argument.

RSEM requires the alignments of a read to be adjacent. For paired-end
reads, RSEM also requires the two mates of any alignment be
adjacent. To check if your SAM/BAM/CRAM file satisfy the requirements,
run

    rsem-sam-validator <input.sam/input.bam/input.cram>

If your file does not satisfy the requirements, you can use
`convert-sam-for-rsem` to convert it into a BAM file which RSEM can
process. Run
 
    convert-sam-for-rsem --help

to get usage information or visit the [convert-sam-for-rsem
documentation
page](convert-sam-for-rsem.html).

Note that RSEM does ** not ** support gapped alignments. So make sure
that your aligner does not produce alignments with
intersions/deletions. In addition, you should make sure that you use
`reference_name.idx.fa`, which is generated by RSEM, to build your
aligner's indices.

### III. Visualization

RSEM includes a copy of SAMtools. When `--no-bam-output` is not
specified and `--sort-bam-by-coordinate` is specified, RSEM will
produce these three files:`sample_name.transcript.bam`, the unsorted
BAM file, `sample_name.transcript.sorted.bam` and
`sample_name.transcript.sorted.bam.bai` the sorted BAM file and
indices generated by the SAMtools included. All three files are in
transcript coordinates. When users in addition specify the
`--output-genome-bam` option, RSEM will produce three more files:
`sample_name.genome.bam`, the unsorted BAM file,
`sample_name.genome.sorted.bam` and
`sample_name.genome.sorted.bam.bai` the sorted BAM file and
indices. All these files are in genomic coordinates.

#### a) Converting transcript BAM file into genome BAM file

Normally, RSEM will do this for you via `--output-genome-bam` option
of `rsem-calculate-expression`. However, if you have run
`rsem-prepare-reference` and use `reference_name.idx.fa` to build
indices for your aligner, you can use `rsem-tbam2gbam` to convert your
transcript coordinate BAM alignments file into a genomic coordinate
BAM alignments file without the need to run the whole RSEM
pipeline.

Usage:

    rsem-tbam2gbam reference_name unsorted_transcript_bam_input genome_bam_output

reference_name	   		  : The name of reference built by `rsem-prepare-reference`				
unsorted_transcript_bam_input	  : This file should satisfy: 1) the alignments of a same read are grouped together, 2) for any paired-end alignment, the two mates should be adjacent to each other, 3) this file should not be sorted by samtools 
genome_bam_output		  : The output genomic coordinate BAM file's name

#### b) Generating a Wiggle file

A wiggle plot representing the expected number of reads overlapping
each position in the genome/transcript set can be generated from the
sorted genome/transcript BAM file output.  To generate the wiggle
plot, run the `rsem-bam2wig` program on the
`sample_name.genome.sorted.bam`/`sample_name.transcript.sorted.bam` file.

Usage:    

    rsem-bam2wig sorted_bam_input wig_output wiggle_name [--no-fractional-weight]

sorted_bam_input        : Input BAM format file, must be sorted  
wig_output              : Output wiggle file's name, e.g. output.wig  
wiggle_name             : The name of this wiggle plot  
--no-fractional-weight  : If this is set, RSEM will not look for "ZW" tag and each alignment appeared in the BAM file has weight 1. Set this if your BAM file is not generated by RSEM. Please note that this option must be at the end of the command line

#### c) Loading a BAM and/or Wiggle file into the UCSC Genome Browser or Integrative Genomics Viewer(IGV)

For UCSC genome browser, please refer to the [UCSC custom track help page](http://genome.ucsc.edu/goldenPath/help/customTrack.html).

For integrative genomics viewer, please refer to the [IGV home page](http://www.broadinstitute.org/software/igv/home). Note: Although IGV can generate read depth plot from the BAM file given, it cannot recognize "ZW" tag RSEM puts. Therefore IGV counts each alignment as weight 1 instead of the expected weight for the plot it generates. So we recommend to use the wiggle file generated by RSEM for read depth visualization.

Here are some guidance for visualizing transcript coordinate files using IGV:

1) Import the transcript sequences as a genome 

Select File -> Import Genome, then fill in ID, Name and Fasta file. Fasta file should be `reference_name.idx.fa`. After that, click Save button. Suppose ID is filled as `reference_name`, a file called `reference_name.genome` will be generated. Next time, we can use: File -> Load Genome, then select `reference_name.genome`.

2) Load visualization files

Select File -> Load from File, then choose one transcript coordinate visualization file generated by RSEM. IGV might require you to convert wiggle file to tdf file. You should use igvtools to perform this task. One way to perform the conversion is to use the following command:

    igvtools tile reference_name.transcript.wig reference_name.transcript.tdf reference_name.genome   
 
#### d) Generating Transcript Wiggle Plots

To generate transcript wiggle plots, you should run the
`rsem-plot-transcript-wiggles` program.  Run 

    rsem-plot-transcript-wiggles --help

to get usage information or visit the [rsem-plot-transcript-wiggles
documentation page](rsem-plot-transcript-wiggles.html).

#### e) Visualize the model learned by RSEM

RSEM provides an R script, `rsem-plot-model`, for visulazing the model learned.

Usage:
    
    rsem-plot-model sample_name output_plot_file

sample_name: the name of the sample analyzed    
output_plot_file: the file name for plots generated from the model. It is a pdf file    

The plots generated depends on read type and user configuration. It
may include fragment length distribution, mate length distribution,
read start position distribution (RSPD), quality score vs observed
quality given a reference base, position vs percentage of sequencing
error given a reference base and alignment statistics.

fragment length distribution and mate length distribution: x-axis is fragment/mate length, y axis is the probability of generating a fragment/mate with the associated length

RSPD: Read Start Position Distribution. x-axis is bin number, y-axis is the probability of each bin. RSPD can be used as an indicator of 3' bias

Quality score vs. observed quality given a reference base: x-axis is Phred quality scores associated with data, y-axis is the "observed quality", Phred quality scores learned by RSEM from the data. Q = -10log_10(P), where Q is Phred quality score and P is the probability of sequencing error for a particular base

Position vs. percentage sequencing error given a reference base: x-axis is position and y-axis is percentage sequencing error

Alignment statistics: It includes a histogram and a pie chart. For the histogram, x-axis shows the number of **isoform-level** alignments a read has and y-axis provides the number of reads with that many alignments. The inf in x-axis means number of reads filtered due to too many alignments. For the pie chart, four categories of reads --- unalignable, unique, **isoform-level**multi-mapping, filtered -- are plotted and their percentages are noted. In both the histogram and the piechart, numbers belong to unalignable, unique, multi-mapping, and filtered are colored as green, blue, gray and red. 
 
## <a name="example-main"></a> Example

Suppose we download the mouse genome from UCSC Genome Browser.  We do
not add poly(A) tails and use `/ref/mouse_0` as the reference name.
We have a FASTQ-formatted file, `mmliver.fq`, containing single-end
reads from one sample, which we call `mmliver_single_quals`.  We want
to estimate expression values by using the single-end model with a
fragment length distribution. We know that the fragment length
distribution is approximated by a normal distribution with a mean of
150 and a standard deviation of 35. We wish to generate 95%
credibility intervals in addition to maximum likelihood estimates.
RSEM will be allowed 1G of memory for the credibility interval
calculation.  We will visualize the probabilistic read mappings
generated by RSEM on UCSC genome browser. We will generate a list of
genes` transcript wiggle plots in `output.pdf`. The list is
`gene_ids.txt`. We will visualize the models learned in
`mmliver_single_quals.models.pdf`

The commands for this scenario are as follows:

    rsem-prepare-reference --gtf mm9.gtf --transcript-to-gene-map knownIsoforms.txt --bowtie --bowtie-path /sw/bowtie /data/mm9 /ref/mouse_0
    rsem-calculate-expression --bowtie-path /sw/bowtie --phred64-quals --fragment-length-mean 150.0 --fragment-length-sd 35.0 -p 8 --output-genome-bam --calc-ci --ci-memory 1024 /data/mmliver.fq /ref/mouse_0 mmliver_single_quals
    rsem-bam2wig mmliver_single_quals.sorted.bam mmliver_single_quals.sorted.wig mmliver_single_quals
    rsem-plot-transcript-wiggles --gene-list --show-unique mmliver_single_quals gene_ids.txt output.pdf 
    rsem-plot-model mmliver_single_quals mmliver_single_quals.models.pdf

## <a name="simulation"></a> Simulation

RSEM provides users the `rsem-simulate-reads` program to simulate RNA-Seq data based on parameters learned from real data sets. Run

    rsem-simulate-reads

to get usage information or read the following subsections.
 
### Usage: 

    rsem-simulate-reads reference_name estimated_model_file estimated_isoform_results theta0 N output_name [-q]

__reference_name:__ The name of RSEM references, which should be already generated by `rsem-prepare-reference`   	     

__estimated_model_file:__ This file describes how the RNA-Seq reads will be sequenced given the expression levels. It determines what kind of reads will be simulated (single-end/paired-end, w/o quality score) and includes parameters for fragment length distribution, read start position distribution, sequencing error models, etc. Normally, this file should be learned from real data using `rsem-calculate-expression`. The file can be found under the `sample_name.stat` folder with the name of `sample_name.model`. `model_file_description.txt` provides the format and meanings of this file.    

__estimated_isoform_results:__ This file contains expression levels for all isoforms recorded in the reference. It can be learned using `rsem-calculate-expression` from real data. The corresponding file users want to use is `sample_name.isoforms.results`. If simulating from user-designed expression profile is desired, start from a learned `sample_name.isoforms.results` file and only modify the `TPM` column. The simulator only reads the TPM column. But keeping the file format the same is required. If the RSEM references built are aware of allele-specific transcripts, `sample_name.alleles.results` should be used instead.   

__theta0:__ This parameter determines the fraction of reads that are coming from background "noise" (instead of from a transcript). It can also be estimated using `rsem-calculate-expression` from real data. Users can find it as the first value of the third line of the file `sample_name.stat/sample_name.theta`.   

__N:__ The total number of reads to be simulated. If `rsem-calculate-expression` is executed on a real data set, the total number of reads can be found as the 4th number of the first line of the file `sample_name.stat/sample_name.cnt`.   

__output_name:__ Prefix for all output files.   

__--seed seed:__ Set seed for the random number generator used in simulation. The seed should be a 32-bit unsigned integer.

__-q:__ Set it will stop outputting intermediate information.   

### Outputs:

output_name.sim.isoforms.results, output_name.sim.genes.results: Expression levels estimated by counting where each simulated read comes from.
output_name.sim.alleles.results: Allele-specific expression levels estimated by counting where each simulated read comes from.

output_name.fa if single-end without quality score;   
output_name.fq if single-end with quality score;   
output_name_1.fa & output_name_2.fa if paired-end without quality
score;   
output_name_1.fq & output_name_2.fq if paired-end with quality score.   

**Format of the header line**: Each simulated read's header line encodes where it comes from. The header line has the format:

    {>/@}_rid_dir_sid_pos[_insertL]

__{>/@}:__ Either '>' or '@' must appear. '>' appears if FASTA files are generated and '@' appears if FASTQ files are generated

__rid:__ Simulated read's index, numbered from 0   

__dir:__ The direction of the simulated read. 0 refers to forward strand ('+') and 1 refers to reverse strand ('-')   

__sid:__ Represent which transcript this read is simulated from. It ranges between 0 and M, where M is the total number of transcripts. If sid=0, the read is simulated from the background noise. Otherwise, the read is simulated from a transcript with index sid. Transcript sid's transcript name can be found in the `transcript_id` column of the `sample_name.isoforms.results` file (at line sid + 1, line 1 is for column names)   

__pos:__ The start position of the simulated read in strand dir of transcript sid. It is numbered from 0   

__insertL:__ Only appear for paired-end reads. It gives the insert length of the simulated read.   

### Example:

Suppose we want to simulate 50 millon single-end reads with quality scores and use the parameters learned from [Example](#example-main). In addition, we set theta0 as 0.2 and output_name as `simulated_reads`. The command is:

    rsem-simulate-reads /ref/mouse_0 mmliver_single_quals.stat/mmliver_single_quals.model mmliver_single_quals.isoforms.results 0.2 50000000 simulated_reads

## <a name="gen_trinity"></a> Generate Transcript-to-Gene-Map from Trinity Output

For Trinity users, RSEM provides a perl script to generate transcript-to-gene-map file from the fasta file produced by Trinity.

### Usage:

    extract-transcript-to-gene-map-from-trinity trinity_fasta_file map_file

trinity_fasta_file: the fasta file produced by trinity, which contains all transcripts assembled.    
map_file: transcript-to-gene-map file's name.    

## <a name="de"></a> Differential Expression Analysis

Popular differential expression (DE) analysis tools such as edgeR and
DESeq do not take variance due to read mapping uncertainty into
consideration. Because read mapping ambiguity is prevalent among
isoforms and de novo assembled transcripts, these tools are not ideal
for DE detection in such conditions.

EBSeq, an empirical Bayesian DE analysis tool developed in UW-Madison,
can take variance due to read mapping ambiguity into consideration by
grouping isoforms with parent gene's number of isoforms. In addition,
it is more robust to outliers. For more information about EBSeq
(including the paper describing their method), please visit [EBSeq's
website](http://www.biostat.wisc.edu/~ningleng/EBSeq_Package).


RSEM includes EBSeq in its folder named `EBSeq`. To use it, first type

    make ebseq

to compile the EBSeq related codes. 

EBSeq requires gene-isoform relationship for its isoform DE
detection. However, for de novo assembled transcriptome, it is hard to
obtain an accurate gene-isoform relationship. Instead, RSEM provides a
script `rsem-generate-ngvector`, which clusters transcripts based on
measures directly relating to read mappaing ambiguity. First, it
calculates the 'unmappability' of each transcript. The 'unmappability'
of a transcript is the ratio between the number of k mers with at
least one perfect match to other transcripts and the total number of k
mers of this transcript, where k is a parameter. Then, Ng vector is
generated by applying Kmeans algorithm to the 'unmappability' values
with number of clusters set as 3. This program will make sure the mean
'unmappability' scores for clusters are in ascending order. All
transcripts whose lengths are less than k are assigned to cluster
3. Run

    rsem-generate-ngvector --help

to get usage information or visit the [rsem-generate-ngvector
documentation
page](rsem-generate-ngvector.html).

If your reference is a de novo assembled transcript set, you should
run `rsem-generate-ngvector` first. Then load the resulting
`output_name.ngvec` into R. For example, you can use 

    NgVec <- scan(file="output_name.ngvec", what=0, sep="\n")

. After that, set "NgVector = NgVec" for your differential expression
test (either `EBTest` or `EBMultiTest`).


For users' convenience, RSEM also provides a script
`rsem-generate-data-matrix` to extract input matrix from expression
results:

    rsem-generate-data-matrix sampleA.[genes/isoforms].results sampleB.[genes/isoforms].results ... > output_name.counts.matrix

The results files are required to be either all gene level results or
all isoform level results. You can load the matrix into R by

    IsoMat <- data.matrix(read.table(file="output_name.counts.matrix"))

before running either `EBTest` or `EBMultiTest`.

Lastly, RSEM provides two scripts, `rsem-run-ebseq` and
`rsem-control-fdr`, to help users find differential expressed
genes/transcripts. First, `rsem-run-ebseq` calls EBSeq to calculate related statistics
for all genes/transcripts. Run 

    rsem-run-ebseq --help

to get usage information or visit the [rsem-run-ebseq documentation
page](rsem-run-ebseq.html). Second,
`rsem-control-fdr` takes `rsem-run-ebseq` 's result and reports called
differentially expressed genes/transcripts by controlling the false
discovery rate. Run

    rsem-control-fdr --help

to get usage information or visit the [rsem-control-fdr documentation
page](rsem-control-fdr.html). These
two scripts can perform DE analysis on either 2 conditions or multiple
conditions.

Please note that `rsem-run-ebseq` and `rsem-control-fdr` use EBSeq's
default parameters. For advanced use of EBSeq or information about how
EBSeq works, please refer to [EBSeq's
manual](http://www.bioconductor.org/packages/devel/bioc/vignettes/EBSeq/inst/doc/EBSeq_Vignette.pdf).

Questions related to EBSeq should
be sent to <a href="mailto:nleng@wisc.edu">Ning Leng</a>.

## <a name="pRSEM"></a> Prior-Enhanced RSEM (pRSEM)

### I. Overview

[Prior-enhanced RSEM (pRSEM)](https://deweylab.github.io/pRSEM/) uses complementary information (e.g. ChIP-seq data) to allocate RNA-seq multi-mapping fragments. We included pRSEM code in the subfolder `pRSEM/` as well as in RSEM's scripts `rsem-prepare-reference` and `rsem-calculate-expression`. 

### II. Demo

To get a quick idea on how to use pRSEM, you can try [this demo](https://github.com/pliu55/pRSEM_demo). It provides a single script, named `run_pRSEM_demo.sh`, which allows you to run all pRSEM's functions. It also contains detailed descriptions of pRSEM's workflow, input and output files.

### III. Installation

To compile pRSEM, type

    make pRSEM

Note that you need to first compile RSEM before compiling pRSEM. Currently, pRSEM has only been tested on Linux.


### IV. Example

To run pRSEM on the [RSEM example above](#example-main), you need to provide:
- __ChIP-seq sequencing file(s) in FASTQ format__ or __a ChIP-seq peak file in BED format__. They will be used by pRSEM to obtain complementatry information for allocating RNA-seq multi-mapping fragments.
- __a genome mappability file in bigWig format__ to let pRSEM build a training
  set of isoforms to learn prior. Mappability can be obtained from UCSC's 
  ENCODE composite track for [human hg19](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig) 
  and [mouse mm9](http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig). For other genomes, you 
  can generate the mappability file by following [this tutorial] (http://wiki.bits.vib.be/index.php/Create_a_mappability_track#Install_and_run_the_GEM_library_tools).

Assuming you would like to use RNA Pol II's ChIP-seq sequencing files `/data/mmliver_PolIIRep1.fq.gz` and `/data/mmliver_PolIIRep2.fq.gz`, with ChIP-seq control `/data/mmliver_ChIPseqCtrl.fq.gz`. Also, assuming the mappability file for mouse genome is `/data/mm9.bigWig` and you prefer to use STAR located at `/sw/STAR` to align RNA-seq fragments and use Bowtie to align ChIP-seq reads. Then, you can use the following commands to run pRSEM:

    rsem-prepare-reference --gtf mm9.gtf \
                           --star \
                           --star-path /sw/STAR \
                           -p 8 \
                           --prep-pRSEM \
                           --bowtie-path /sw/bowtie \
                           --mappability-bigwig-file /data/mm9.bigWig \
                           /data/mm9 \
                           /ref/mouse_0
  
    rsem-calculate-expression --star \
                              --star-path /sw/STAR \
                              --calc-pme \
                              --run-pRSEM \
                              --chipseq-target-read-files /data/mmliver_PolIIRep1.fq.gz,/data/mmliver_PolIIRep2.fq.gz \
                              --chipseq-control-read-files /data/mmliver_ChIPseqCtrl.fq.gz \
                              --bowtie-path /sw/bowtie \
                              -p 8 \
                              /data/mmliver.fq \
                              /ref/mouse_0 \
                              mmliver_single_quals


To find out more about pRSEM options and examples, you can use the commands:

    rsem-prepare-reference --help

and 

    rsem-calculate-expression --help


### V. System Requirements
- Linux
- Perl version >= 5.8.8
- Python version >= 2.7.3
- R version >= 3.3.1
- Bioconductor 3.3


### VI. Required External Packages
All the following packages will be automatically installed when compiling pRSEM.
- [data.table 1.9.6](https://cran.r-project.org/web/packages/data.table/index.html): an extension of R's data.frame, heavily used by pRSEM.
- [GenomicRanges 1.24.3](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html): efficient representing and manipulating genomic intervals, heavily used by pRSEM.
- [ShortRead 1.30.0](https://bioconductor.org/packages/release/bioc/html/ShortRead.html): guessing the encoding of ChIP-seq FASTQ file's quality score.
- [caTools 1.17.1](https://cran.r-project.org/web/packages/caTools/index.html): used for SPP Peak Caller.
- [SPP Peak Caller](https://code.google.com/archive/p/phantompeakqualtools/):
  ChIP-seq peak caller. Source code was slightly modified in terms of included headers in order to be compiled under R v3.3.1.
- [IDR](https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz?attredirects=0):
  calculating Irreproducible Discovery Rate to call peaks from multiple ChIP-seq replicates.


## <a name="authors"></a> Authors

[Bo Li](http://bli25ucb.github.io/) and [Colin Dewey](https://www.biostat.wisc.edu/~cdewey/) designed the RSEM algorithm. [Bo Li](http://bli25ucb.github.io/) implemented the RSEM software. [Peng Liu](https://www.biostat.wisc.edu/~cdewey/group.html) contributed the STAR aligner options and prior-enhanced RSEM (pRSEM).

## <a name="acknowledgements"></a> Acknowledgements

RSEM uses the [Boost C++](http://www.boost.org/) and
[SAMtools](http://www.htslib.org/) libraries. RSEM includes
[EBSeq](http://www.biostat.wisc.edu/~ningleng/EBSeq_Package/) for
differential expression analysis.

We thank earonesty, Dr. Samuel Arvidsson, John Marshall, and Michael
R. Crusoe for contributing patches.

We thank Han Lin, j.miller, Jo&euml;l Fillon, Dr. Samuel G. Younkin,
Malcolm Cook, Christina Wells, Uro&#353; &#352;ipeti&#263;,
outpaddling, rekado, and Josh Richer for suggesting possible fixes.

**Note** that `bam_sort.c` of SAMtools is slightly modified so that
  `samtools sort -n` will not move the two mates of paired-end
  alignments apart. In addition, we turn on the `--without-curses`
  option when configuring SAMtools and thus SAMtools' curses-based
  `tview` subcommand is not built.

## <a name="license"></a> License

RSEM is licensed under the [GNU General Public License
v3](http://www.gnu.org/licenses/gpl-3.0.html).
