#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);

my ($min_reads) = (10);
GetOptions (
	"min-reads:s" => \$min_reads,
);

while (<>) {
	my $line = $_;
	my @columns = split("\t",$line);
	my $chr = $columns[0];
	my $start = $columns[1];
	my $end = $start + 1;
	my $num_reads = $columns[3];
	my $calls = $columns[4];
	my $id = "mpileup_number_" . $.; # you can use this in the fourth column of the output
	if($num_reads < $min_reads){ # not enough coverage to have good confidence in the call
		next;
	}
	my $num_ref = 0;
	while ($calls =~ /[,.]/g) { $num_ref++ }
	my $num_var = $num_reads - $num_ref;
	my $varAlleleFreq;
	if($num_reads == 0){
		$varAlleleFreq = "NA";
	}
	else {
		$varAlleleFreq = ($num_var/$num_reads)*100;
	}	
	print("$chr\t$start\t$end\t$num_reads\t$varAlleleFreq\n");
}

