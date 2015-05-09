#!/bin/env perl

#
#  pliu 20150508
#
#  test run for pRSEM
#

use strict;
use warnings;

my $starbin      = '/tier2/deweylab/pliu/local/STAR-2.4.0h/STAR';
my $bowtiepath   = '/tier2/deweylab/pliu/local/bin';
my $datdir       = '/tier2/deweylab/scratch/pliu/dev';
my $fgtf         = "$datdir/chr1819.gtf";
my $genomedir    = "$datdir/mm10";
my $star_gnmdir  = "$datdir/star_genome";
my $star_alndir  = "$datdir/star_aln";
my $rsem_refdir  = "$datdir/rsem_ref";
my $rsem_exprdir = "$datdir/rsem_expr";
my $nthreads     = 16;
my $runid        = 'test';

my $chipseq_input_rep1  = "$datdir/MelInputIggmusRep1.fastq.gz";
my $chipseq_target_rep1 = "$datdir/MelPol2IggmusRep1.fastq.gz";
my $chipseq_target_rep2 = "$datdir/MelPol2IggmusRep2.fastq.gz";
my $rnaseq_rd1 = "$datdir/MelRibozerogRep1Rd1.fastq.gz";
my $rnaseq_rd2 = "$datdir/MelRibozerogRep1Rd2.fastq.gz";


sub main {
 # prepRSEMRef();

 # prepStarGenome();

 # runStarAlignment();

   runRSEM();
}


sub runRSEM {
  my $bin = '/ua/pliu/dev/RSEM/rsem-calculate-expression';

  if ( -e $rsem_exprdir ) {
    system("rm -fr $rsem_exprdir");
  }
  mkdir $rsem_exprdir;
  chdir $rsem_exprdir;

  my $nsteps = 1e+4;
  my $fbam = "$star_alndir/Aligned.toTranscriptome.out.bam";

  my $cmd = "$bin --bam ".
                ' --estimate-rspd '.
                ' --no-bam-output '.
                ' --seed 12345 '.
                " --num-threads $nthreads ".
                ' --paired-end '.
                ' --forward-prob 0  '.

                ## added options by myself for PME
                ## save .ofg, .omit file for Gibbs sampling
                ' --keep-intermediate-files '. 
                ' --calc-pme '.
                " --gibbs-number-of-samples $nsteps" .
                ' --quiet' .
                ######

                " $fbam $rsem_refdir/$runid $runid "; 

 #print "$cmd\n";
  system($cmd);
}


sub runStarAlignment {
  if ( not -e $star_alndir ) {
    mkdir $star_alndir;
  }
  chdir $star_alndir;

  my $cmd = "$starbin ".
              " --genomeDir $star_gnmdir " .
              " --readFilesIn $rnaseq_rd1 $rnaseq_rd2 " .
              ' --outSAMunmapped Within   ' .
              ' --outFilterType BySJout   ' .
              ' --outSAMattributes NH HI AS NM MD ' .
              ' --outFilterMultimapNmax 20        ' .
              ' --outFilterMismatchNmax 999       ' .
              ' --outFilterMismatchNoverLmax 0.04 ' .
              ' --alignIntronMin 20        ' .
              ' --alignIntronMax 1000000   ' .
              ' --alignMatesGapMax 1000000 ' .
              ' --alignSJoverhangMin 8     ' .
              ' --alignSJDBoverhangMin 1   ' .
              ' --sjdbScore 1              ' .
              ' --readFilesCommand zcat ' .
              " --runThreadN $nthreads " .
              ' --outSAMtype BAM Unsorted ' .
              ' --quantMode TranscriptomeSAM ' .
              ' --outSAMheaderHD \@HD VN:1.4 SO:coordinate ';

  system($cmd);
}


sub prepStarGenome {
  if ( not -e $star_gnmdir ) {
    mkdir $star_gnmdir;
  }

  my @fasta_files = (<$genomedir/*.fa>, <$genomedir/*.fasta>);

  my $cmd = "$starbin --runThreadN $nthreads " .
                    " --runMode genomeGenerate " .
                    " --genomeDir $star_gnmdir       " .
                    " --genomeFastaFiles @fasta_files " .
                    " --sjdbGTFfile $fgtf             " .
                    " --sjdbOverhang 100              " .
                    " --outFileNamePrefix $star_gnmdir ";
  print $cmd, "\n";
  system($cmd);
}


sub prepRSEMRef {
  my $bin = '/ua/pliu/dev/RSEM/rsem-prepare-reference';

  if ( -e $rsem_refdir ) {
    system("rm -fr $rsem_refdir");
  }
  
  mkdir $rsem_refdir;
  chdir $rsem_refdir;

  my $cmd = "$bin --gtf $fgtf " .
                " --bowtie-path $bowtiepath " .
                ' --index-genome ' .
                "$genomedir $rsem_refdir/$runid";

 #print "\n$cmd\n\n";
  system($cmd)
}

main();
