#!/usr/bin/env perl

package rsem_perl_utils;

use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(runCommand);
our @EXPORT_OK = qw(runCommand collectResults showVersionInfo getSAMTOOLS hasPolyA);

my $version = "RSEM v1.3.1"; # Update version info here
my $samtools = "samtools-1.3"; # If update to another version of SAMtools, need to change this

# command, {err_msg}
sub runCommand {
    print $_[0]."\n";
    my $status = system($_[0]);

    if ($? == -1) {
	my @arr = split(/[ \t]+/, $_[0]);
	print "$arr[0] : $!!\n";
	print "Please check if you have compiled the associated codes by typing related \"make\" commands and/or made related executables ready to use.\n";
	exit(-1);
    }

    if ($status != 0) {
        my $errmsg = "";
        if (scalar(@_) > 1) { $errmsg .= $_[1]."\n"; }
	$errmsg .= "\"$_[0]\" failed! Plase check if you provide correct parameters/options for the pipeline!\n";
	print $errmsg;
        exit(-1);
    }
    print "\n";
}

my @allele_title = ("allele_id", "transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "AlleleIsoPct", "AlleleGenePct", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "AlleleIsoPct_from_pme_TPM", "AlleleGenePct_from_pme_TPM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation");

my @transcript_title = ("transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "IsoPct", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "IsoPct_from_pme_TPM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation");

my @gene_title = ("gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM", "FPKM", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "TPM_coefficient_of_quartile_variation", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound", "FPKM_coefficient_of_quartile_variation");

# type, inpF, outF
sub collectResults {
    my $local_status;
    my ($inpF, $outF);
    my @results = ();
    my $line;

    $inpF = $_[1];
    $outF = $_[2];

    $local_status = open(INPUT, $inpF);
    if ($local_status == 0) { print "Fail to open file $inpF!\n"; exit(-1); }
    
    @results = ();
    
    while ($line = <INPUT>) {
	chomp($line);
	my @local_arr = split(/\t/, $line);
	push(@results, \@local_arr); 
    }

    close(INPUT);

    $local_status = open(OUTPUT, ">$outF");
    if ($local_status == 0) { print "Fail to create file $outF!\n"; exit(-1); }

    my $n = scalar(@results);
    my $m = scalar(@{$results[0]});

    $" = "\t";

    my @out_arr = ();
    for (my $i = 0; $i < $n; $i++) {
	if ($_[0] eq "allele") { push(@out_arr, $allele_title[$i]); }
	elsif ($_[0] eq "isoform") { push(@out_arr, $transcript_title[$i]); }
	elsif ($_[0] eq "gene") { push(@out_arr, $gene_title[$i]); }
	else { print "A bug on 'collectResults' is detected!\n"; exit(-1); }
    }
    print OUTPUT "@out_arr\n";

    for (my $i = 0; $i < $m; $i++) {
	@out_arr = ();
	for (my $j = 0; $j < $n; $j++) { push(@out_arr, $results[$j][$i]); }
	print OUTPUT "@out_arr\n"; 
    }

    close(OUTPUT);
}

sub showVersionInfo {
    print "Current version: $version\n";
    exit(0);
}

sub getSAMTOOLS {
    return $samtools;
}

sub hasPolyA {
    open(INPUT, $_[0]);
    my $line = <INPUT>; chomp($line);
    close(INPUT);
    my ($fullLen, $totLen) = split(/ /, $line);
    return $fullLen < $totLen;
}

1;
