#!/bin/bash 

#
#  pliu 20150507
#
#  prepare a small example for testing pRSEM
#
#  take mm10's two smallest regular chroms: chr18 and chr19
#

nlines=100000

mousedir='/tier2/deweylab/scratch/pliu/mouse/'
chipseqdir="$mousedir/chipseq/01_download/"
rnaseqdir="$mousedir/rnaseq/01_download/"
fgtf="$mousedir/rnaseq/02_filtered_gnc/female.gtf"

# for id in 'MelInputIggmusRep1' 'MelPol2IggmusRep1' 'MelPol2IggmusRep2'
# do 
#   file="${id}.fastq.gz"
#   zcat $chipseqdir/$file | head -$nlines | gzip - > $file
# done

# for id in 'MelRibozerogRep1Rd1' 'MelRibozerogRep1Rd2' 
# do
#   file="${id}.fastq.gz"
#   zcat $rnaseqdir/$file | head -$nlines | gzip - > $file
# done

cat $fgtf | awk -F"\t" '{if( ($1=="chr18") || ($1=="chr19")) print $0}' > 'chr1819.gtf'
