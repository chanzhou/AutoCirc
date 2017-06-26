#! /usr/bin/perl

use strict;
#08/26/13

my $usage="$0 gtag_circRNA_in_bed exon_boundaries_circRNA_in_bed > of_in_bed\n";

if(@ARGV<2){die $usage;}
open GTAG, $ARGV[0];
my $line;
my %circR2score_infor;

while($line=<GTAG>){
	chomp($line);
	my @fields=split /\t/, $line;
	my $chr=$fields[0];
	my $start=$fields[1];
	my $end=$fields[2];
	my $name=$fields[3];
	my $score=$fields[4];
	my $site_strand=$fields[5];
	#my $read_strand=$fields[6];
	my $distance=$fields[7];
	my $read_tags=$fields[8];

	my $new_line=join "\t", ($chr, $start,$end,$name, $score, $site_strand, $distance, $read_tags);
	my $circR=join ":", ($chr,$start,$end,$site_strand);

	if ( ($chr ne "#chrom") and (not exists $circR2score_infor{$circR})){
		$circR2score_infor{$circR}[0]=$score;
		$circR2score_infor{$circR}[1]=$new_line;
		$circR2score_infor{$circR}[2]=$read_tags;
	}  

}
close GTAG;

open EXONRNA, $ARGV[1];
while($line=<EXONRNA>){
	chomp($line);
	my @fields=split /\t/,$line;
	my $chr=$fields[0];
	my $start=$fields[1];
	my $end=$fields[2];
	my $name=$fields[3];
	my $score=$fields[4];
	my $site_strand=$fields[5];
	my $distance=$fields[7];#revise
	my $read_tags=$fields[8];#revise
	
	my $new_line;
	if($read_tags){
		$new_line=join "\t", ($chr, $start,$end,$name, $score, $site_strand, $distance, $read_tags);
	}else{
		$new_line=join "\t", ($chr, $start,$end,$name, $score, $site_strand, $distance);
	
	}

	my $circR=join ":", ($chr,$start,$end,$site_strand);

	if ( ($chr ne "#chrom") and (not exists $circR2score_infor{$circR})){
		$circR2score_infor{$circR}[0]=$score;
		$circR2score_infor{$circR}[1]=$new_line;
		$circR2score_infor{$circR}[2]=$read_tags;
	}elsif( ($chr ne "#chrom") and (exists $circR2score_infor{$circR})){
		my $previous_score=$circR2score_infor{$circR}[0];
		if($score > $previous_score){
			$circR2score_infor{$circR}[0]=$score;
			if($read_tags){
				$circR2score_infor{$circR}[1]=$new_line;
				$circR2score_infor{$circR}[2]=$read_tags;
			}else{
				my $previous_tags=$circR2score_infor{$circR}[2];
				my $full_new_line=join "|", ($new_line,$previous_tags);
				$circR2score_infor{$circR}[1]=$full_new_line;
			}
		}#end-if-score

	
	}#end-if-elsif  
}

my $key;
foreach $key (sort keys %circR2score_infor){
	print "$circR2score_infor{$key}[1]\n";
}
