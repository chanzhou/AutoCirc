#! /usr/bin/perl -w

#Function: 
#Input: unmapped_sample.sam, 
#output: the X-mers from both ends of each sequences in fastq format.
#Note: the original read should be kept as part of the first anchors identifier 
#to simplify downstream analysis.

use strict;
my $usage="$0 unmapped_sample.sam len_of_mers[default is 20]> of\n or \n$0 unmapped_sample.sam len_of_mers[default is 20] | gzip > of.gz\n";

if(@ARGV<1){die $usage;}

my $len_mer;

if(@ARGV<2){ $len_mer=20;}
else{$len_mer=$ARGV[1];}

open IF, $ARGV[0];
my $line;
while($line=<IF>){
	chomp($line);
	my $seq=0;
	my $quality=0;
	my $identifier=0;
	my $head_seq=0;
	my $tail_seq=0;
	my $head_quality=0;
	my $tail_quality=0;

	if($line=~/^\S+\t\d+\t\S+\t\d+/){
		#warn "line=$line\n";
		my @fields=split "\t", $line;
		$identifier=$fields[0];
		$seq=$fields[9];
		$quality=$fields[10];
		
		my $minus_len_mer=(-1)*$len_mer;

		#fetch the 1st 20mer and last 20mer
		$head_seq=substr $seq, 0, $len_mer;
		$head_quality=substr $quality, 0, $len_mer;

		$tail_seq=substr $seq, $minus_len_mer;
		$tail_quality=substr $quality, $minus_len_mer;
	
		my $new_identifier=join "_", ($identifier, $seq);
		print "\@$new_identifier\n";
		print "$head_seq\n";
		print "+\n";
		print "$head_quality\n";

		print "\@$identifier\n";
		print "$tail_seq\n";
		print "+\n";
		print "$tail_quality\n";
	}

}
close IF;
