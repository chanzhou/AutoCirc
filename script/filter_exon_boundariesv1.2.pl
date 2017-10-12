#! /usr/bin/perl -w

use strict;
#revised on 10/11/2017

#input: **.circ.refseq.[chr*].bed (e.g. Ot2119R1.circ.refseq.chr1.bed)

my $usage="$0 infile window_size > of\n";
if(@ARGV<2){die $usage;}
my $window=$ARGV[1];
open IF, $ARGV[0];
my $line;

my %circRNAs=();

while($line=<IF>){
#chr14	100705790	100705802	*junc434	1	-	15	11	SRR2086054.70425637	chr14	100705101	100745371	NM_003403	+	5	1159,163,61,159,1617	0,23539,35933,37725,38653	YY1
#old-chr10	103432670	103436195	chr10junc35	1	+	3	3524	HWI-D00449:117:HVGGKADXX:1:2216:6375:97086	NM_022039	-	9	747,142,141,66,95,133,186,96,606	103370420,103371384,103372093,103384501,103427642,103432671,103433244,103436097,103454137,	FBXW4
	chomp($line);
	my @fields=split /\t/, $line;
	my $chr=$fields[0];
	my $pseudo_start=$fields[1];
	my $pseudo_end=$fields[2];
	my $junc_name=$fields[3];
	my $coverage=$fields[4];
	my $read_strand=$fields[5];
	my $extension_overlaps=$fields[6];
	my $genomic_distance=$fields[7];
	my $readtags=$fields[8];
	my $gstart=$fields[10];
	my $gend=$fields[11];
	my $gid=$fields[12];
	my $gstrand=$fields[13];
	my $exon_num=$fields[14];
	my @exon_sizes=split /\,/, $fields[15];
	my @relative_exon_starts=split /\,/, $fields[16];
	my @exon_starts;
	for(my $i=0; $i<@exon_sizes;$i++){
		$exon_starts[$i]=$relative_exon_starts[$i]+$gstart;
	}

	my $exonA_index=0;
	my $exonB_index=0;
	my $start_exonA=0;
	my $end_exonB=0;
	
	#filter out by the genomic distances
	if ($genomic_distance >= 100 and $genomic_distance <=100000){
		for(my $i=0; $i<@exon_sizes; $i++){
			if( ($pseudo_start <= ($exon_starts[$i] + $window)) and 
			    (($exon_starts[$i] - $window)<=$pseudo_start)){
				$exonA_index=$i+1;
				$start_exonA=$exon_starts[$i];
			}

			if( ($pseudo_end <=($exon_starts[$i]+$exon_sizes[$i]+$window)) and
			    (($exon_starts[$i]+$exon_sizes[$i]-$window)<=$pseudo_end) ){
		    		$exonB_index=$i+1;
				$end_exonB=$exon_starts[$i]+$exon_sizes[$i];
	   		 }
		}#end-for

		my $num_exon_circRNA=abs($exonB_index-$exonA_index)+1;
		my $distance_to_last_exon=-1;
		if($gstrand eq "+"){
			$distance_to_last_exon=$exon_num-$exonB_index;
		}elsif($gstrand eq "-"){
			#the 1st exon block in the file is actually the last exon of gene since it is in the anti-sense strand
			$distance_to_last_exon=$exonA_index-1;
		}
	
		#store in a hash array (make sure to merge the circR which supported by reads from both strand
		#...to be revised...
		if($exonA_index>0 and $exonB_index>0){
			$genomic_distance=abs($end_exonB-$start_exonA);
			my $span_num_exon=$exonB_index-$exonA_index+1;
			my $coord_gid=join ":", ($start_exonA, $end_exonB, $gid);#key
			if(not exists $circRNAs{$coord_gid}){
				$circRNAs{$coord_gid}[0]=$chr;
				$circRNAs{$coord_gid}[1]=$junc_name;
				$circRNAs{$coord_gid}[2]=$coverage;
				$circRNAs{$coord_gid}[3]=$gstrand;
				$circRNAs{$coord_gid}[4]=$genomic_distance;
				$circRNAs{$coord_gid}[5]=$gid;
				$circRNAs{$coord_gid}[6]=$exon_num;
				$circRNAs{$coord_gid}[7]=$exonA_index;
				$circRNAs{$coord_gid}[8]=$exonB_index;
				$circRNAs{$coord_gid}[9]=$span_num_exon;
				$circRNAs{$coord_gid}[10]=$distance_to_last_exon;
				$circRNAs{$coord_gid}[11]=$read_strand;
				$circRNAs{$coord_gid}[12]=$readtags;
			}elsif(exists $circRNAs{$coord_gid}){
				my $previous_junc_name=$circRNAs{$coord_gid}[1];
				my $new_junc_name;
				if($junc_name le $previous_junc_name){
					$new_junc_name=join ":", ($junc_name, $previous_junc_name);
				}else{
					$new_junc_name=join ":", ($previous_junc_name, $junc_name);
				}
				my $new_coverage=$coverage+$circRNAs{$coord_gid}[2];
				my $pre_readtags=$circRNAs{$coord_gid}[12];
				my $new_readtags=join "|", ($readtags, $pre_readtags);
				
				$circRNAs{$coord_gid}[1]=$new_junc_name;
				$circRNAs{$coord_gid}[2]=$new_coverage;
				$circRNAs{$coord_gid}[12]=$new_readtags;
			}

			#print "$chr\t$start_exonA\t$end_exonB\t$junc_name\t$coverage\t$gstrand\t$genomic_distance\t$gid\t$exon_num\t$exonA_index\t$exonB_index\t$span_num_exon\t$distance_to_last_exon\n";
		}#end-print
	
	}#end-if genomic distance
}
close IF;

my $key;
foreach $key (sort keys %circRNAs){
	my @key_fields=split /:/, $key;
	my $start=$key_fields[0];
	my $end=$key_fields[1];
	my $gid=$key_fields[2];
	print "$circRNAs{$key}[0]\t$start\t$end\t$circRNAs{$key}[1]\t$circRNAs{$key}[2]\t$circRNAs{$key}[3]\t$circRNAs{$key}[11]\t$circRNAs{$key}[4]\t$circRNAs{$key}[12]\t$circRNAs{$key}[5]\t$circRNAs{$key}[6]\t$circRNAs{$key}[7]\t$circRNAs{$key}[8]\t$circRNAs{$key}[9]\t$circRNAs{$key}[10]\n";
}
