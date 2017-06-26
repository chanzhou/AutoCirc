#! /usr/bin/perl 
use strict;
#revised on 10/30/2015

#Function:
#infile: circRNAs overlap or inside annotated genes
#outfile: circRNA info and label-if-start/end meets  the exon boundaries of certain 
my $usage="$0 circ.refseq/2genes.bed > of\n";

if(@ARGV<1){die $usage;}
open IF, $ARGV[0];
my $line;
#chr15	25307114	25312573	hESC_m6A_circRNA_105	2	+	#genomic_dist	#read_tags	chr15	25307478	25307575	NR_003322	+	25307575	25307575	1	97	25307478,	SNORD116-7
my %correct_circ=();
my %problem_circ=();
while($line=<IF>){
	chomp($line);
	my @fields=split /\t/, $line;
	my $circ_start=$fields[1];
	my $circ_end=$fields[2];
	my $circ_strand=$fields[5];
	my $numexons=$fields[15];
	my @exonlens=split /\,/, $fields[16];
	my @exonstarts=split /\,/, $fields[17];
	
	my $start_flag=-1;
	my $end_flag=-1;

	my $junc=join ":", ($circ_start, $circ_end, $circ_strand);

	for(my $i=0;$i<$numexons;$i++){
		if($circ_start == $exonstarts[$i]){$start_flag=$i+1;}
		if($circ_end == ($exonstarts[$i]+$exonlens[$i])){$end_flag=$i+1;}
	}

	my $flag=$start_flag*$end_flag;
	if($start_flag >= 0 and $end_flag >=0 and (not exists $correct_circ{$junc})){
		$correct_circ{$junc}[0]=$line;
	}elsif($start_flag >= 0 and $end_flag >=0 and (exists $correct_circ{$junc})){
		push @{$correct_circ{$junc}}, $line;
	}elsif(($start_flag>0 or $end_flag>0) and ($flag<0) 
		and (not exists $correct_circ{$junc})
		and (not exists $problem_circ{$junc})
		){
			$problem_circ{$junc}[0]=$line;
			warn "$start_flag\t$end_flag\t$flag\n";
	}elsif(($start_flag>0 or $end_flag>0) and ($flag<0) 
		and (not exists $correct_circ{$junc})
		and (exists $problem_circ{$junc})
		){
			push @{$problem_circ{$junc}};$line;	
	}
}
close IF;

my $key;
foreach $key (sort keys %problem_circ){
	if(not exists $correct_circ{$key}){
		for (my $i=0;$i<@{$problem_circ{$key}};$i++){print "$problem_circ{$key}[$i]\n";}

	}
	
}
