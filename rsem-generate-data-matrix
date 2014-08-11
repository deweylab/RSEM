#!/usr/bin/env perl

use strict;

if (scalar(@ARGV) == 0) {
    print "Usage: rsem-generate-data-matrix sampleA.[alleles/genes/isoforms].results sampleB.[alleles/genes/isoforms].results ... > output_name.matrix\n";
    print "All result files should have the same file type. The 'expected_count' columns of every result file are extracted to form the data matrix.\n";
    exit(-1);
}

my $offsite = 4; # for new file formats
if ($ARGV[0] =~ /alleles.results$/) { $offsite = 5; }

my $line;
my $n = scalar(@ARGV);
my $M = -1;
my @matrix = ();

# 0, file_name; 1, reference of expected count array; 2, reference of transcript_id/gene_id array
sub loadData {
    open(INPUT, $_[0]);
    my $line = <INPUT>; # The first line contains only column names
    while ($line = <INPUT>) {
	chomp($line); 
	my @fields = split(/\t/, $line);
	push(@{$_[2]}, "\"$fields[0]\"");
	push(@{$_[1]}, $fields[$offsite]);
    }
    close(INPUT);

    if (scalar(@{$_[1]}) == 0) {
	print STDERR "Nothing is detected! $_[0] may not exist or is empty.\n";
	exit(-1);
    }
}

#0, M; 1, reference of @ids_arr; 2, reference of @ids
sub check {
    my $size = $_[0];
    for (my $i = 0; $i < $size; $i++) { 
	if ($_[1]->[$i] ne $_[2]->[$i]) {
	    return 0;
	}
    }
    return 1;
}

my @ids_arr = ();

for (my $i = 0; $i < $n; $i++) {
    my (@ids, @ecs) = ();
    &loadData($ARGV[$i], \@ecs, \@ids);

    if ($M < 0) { 
	$M = scalar(@ids); 
	@ids_arr = @ids;
    }
    elsif (!&check($M, \@ids_arr, \@ids)) { 
	print STDERR "Number of lines among samples are not equal!\n"; 
	exit(-1); 
    }

    my $colname;
    if (substr($ARGV[$i], 0, 2) eq "./") { $colname = substr($ARGV[$i], 2); }
    else { $colname = $ARGV[$i]; }
    $colname = "\"$colname\"";
    @ecs = ($colname, @ecs);
    push(@matrix, \@ecs);
}

@ids_arr = ("", @ids_arr);
@matrix = (\@ids_arr, @matrix);

for (my $i = 0; $i <= $M; $i++) {
    for (my $j = 0; $j < $n; $j++) { print "$matrix[$j][$i]\t"; }
    print "$matrix[$n][$i]\n";
}
