#! /usr/bin/perl -w

use strict;

my $usage="$0 circ2refseq > of\n";

if(@ARGV<1){die $usage;}
open IF, $ARGV[0];
my $line;
my %circ2genes=();#hash; key: circ_id, value: start and end of linear transcript 
my %circ2info=(); #hash; key: circ_id, value: all information of circ
my %excluded_circ=(); #circ_id, value=0;
my %kept_circ=(); #circ_id;

while($line=<IF>){
	chomp($line);
#chr10	101659675	101659823	chr10junc23	1	-	148	HWI-D00449:117:HVGGKADXX:1:2109:13540:50069	chr10	101635333	101769676	NM_015221	-	101636907	101731881	17	1760,551,199,180,333,129,105,131,200,18,148,100,194,1992,123,155,82	101635333,101639567,101643767,101645443,101646056,101648581,101654702,101656023,101657842,101658499,101659675,101667751,101668709,101714970,101728871,101731736,101769594,	DNMBP
	my @fields=split /\t/,$line;
	my $chr=$fields[0];
	my $start=$fields[1];
	my $end=$fields[2];
	my $circ_id=$fields[3];
	my $score=$fields[4];
	my $strand=$fields[5];
	my $size=$fields[6];
	my $reads=$fields[7];
	my $tx_start=$fields[9];
	my $tx_end=$fields[10];
	my $tx_id=$fields[11];

	my $circ_info=join "\t", ($chr, $start, $end, $circ_id, $score, $strand, $size, $reads);
	if(not exists $circ2info{$circ_id}){
		$circ2info{$circ_id}=$circ_info;
		#warn "circ_info: \n$circ2info{$circ_id}\n\n";
	}

	if(not exists $circ2genes{$circ_id}){
		$circ2genes{$circ_id}[0]=$tx_start;
		$circ2genes{$circ_id}[1]=$tx_end;
		$circ2genes{$circ_id}[2]=$tx_id;
	}elsif( (not exists $excluded_circ{$circ_id}) and ($tx_start < $circ2genes{$circ_id}[1] and $circ2genes{$circ_id}[0] < $tx_end)){ #overlap
		my $max_tx_end  = $tx_end >= $circ2genes{$circ_id}[1] ? $tx_end : $circ2genes{$circ_id}[1] ;
		my $min_tx_start= $tx_start <=$circ2genes{$circ_id}[0] ? $tx_start : $circ2genes{$circ_id}[0];

		$circ2genes{$circ_id}[0]= $min_tx_start;
		$circ2genes{$circ_id}[1]= $max_tx_end;
		my $tx_ids=join ":", ($circ2genes{$circ_id}[2], $tx_id);
		$circ2genes{$circ_id}[2]=$tx_ids;

		#warn "tx_start=$tx_start\ttx_end=$tx_end\n";
		#warn "tx_gene: start -end : $circ2genes{$circ_id}[0]\t$circ2genes{$circ_id}[1]\n";
		#warn "min=$min_tx_start\tmax=$max_tx_end\n";;

	}else{ #not overlap, exlude this circ	
		if(not exists $excluded_circ{$circ_id}) { 
			$excluded_circ{$circ_id}=0; 
			#warn "$circ_id:\t$circ2genes{$circ_id}[2]\t$tx_id\n";
		}	
	}

	if((not exists $kept_circ{$circ_id}) and ($tx_start <= $start and $end <= $tx_end) ){
			$kept_circ{$circ_id}=$tx_id;
	}
}
close IF;

my $key;
foreach $key (sort keys %excluded_circ){
	if (not exists $kept_circ{$key}){
		print "$circ2info{$key}\n";
	}
}

