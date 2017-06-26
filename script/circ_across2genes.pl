#! /usr/bin/perl -w

use strict;

my $usage="$0 circ2refseq > of\n";

if(@ARGV<1){die $usage;}
open IF, $ARGV[0];
my $line;
my %circ2genes=();#hash-array: gene-name, No. of genes
my %circ2info=(); #hash-array: info

while($line=<IF>){
	chomp($line);
#chr10	101659675	101659823	chr10junc23	1	-	148	HWI-D00449:117:HVGGKADXX:1:2109:13540:50069	chr10	101635333	101769676	NM_015221	-	101636907	101731881	17	1760,551,199,180,333,129,105,131,200,18,148,100,194,1992,123,155,82	101635333,101639567,101643767,101645443,101646056,101648581,101654702,101656023,101657842,101658499,101659675,101667751,101668709,101714970,101728871,101731736,101769594,	DNMBP
	my @fields=split /\t/,$line;
	my $circ_id=$fields[3];
	my $gname=$fields[18];
	if(not exists $circ2genes{$circ_id}){
		$circ2genes{$circ_id}[0]=$gname;
		$circ2info{$circ_id}[0]=$line;
	}elsif( (exists $circ2genes{$circ_id}) and ($gname ne $circ2genes{$circ_id}[0])){
		push @{$circ2genes{$circ_id}}, $gname;
		push @{$circ2info{$circ_id}}, $line;
		#warn "$circ_id\t$circ2genes{$circ_id}[0]\t$gname\n";
	}elsif( (exists $circ2genes{$circ_id}) and ($gname eq $circ2genes{$circ_id}[0])){
		push @{$circ2info{$circ_id}}, $line;
	}
}
close IF;

my $key;
foreach $key (sort keys %circ2genes){
	my $num=@{$circ2genes{$key}};
	if($num>1){
		for (my $i=0;$i<@{$circ2info{$key}};$i++){print "$circ2info{$key}[$i]\n";}
	}
}
