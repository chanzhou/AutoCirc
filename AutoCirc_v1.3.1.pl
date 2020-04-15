#! /usr/bin/perl -w

###############################################################################
#Copyright 2017-2022 Chan Zhou
###############################################################################

#Program: AutoCirc is a computational tool to Automatically detect Circular RNAs (circRNAs) from RNA-Seq data. 
#website: https://github.com/chanzhou/AutoCirc 
#Author: Chan Zhou
#Contact: zhouchan99@gmail.com
#
#Reference:
#Identification and characterization of m6A circular RNA epitranscriptomes
#Chan Zhou, Benoit Molinie, Kaveh Daneshvar, Joshua V. Pondick, Jinkai Wang, 
#Nicholas O. Van Wittenberghe, Yi Xing, Cosmas C. Giallourakis, and Alan C. Mullen
#(submitted)
#BioRxiv version: http://biorxiv.org/content/early/2017/03/10/115899
#Last updated: 04/09/2018
#Pls read README before using this pipeline.

our $VERSION="version 1.3.1";

use strict;
use POSIX;
use warnings; 
use Getopt::Long qw(GetOptions);

my $USAGE =<<EOF;
Usage:
	[perl] $0 [options]

Options:
	
	-g/--genome     Absolute path of the reference genome file in fasta format (e.g. /genome_reference/hg19.fa) 
	-I/--gindex     Absolute path of the genome index of bowtie2 (e.g. /bowtie2_index/hg19 )
	--bam		Bam file of unmapped reads (output by bowtie/bowtie2/tophat)
	-b/--bed	Reference gene annotation in standard bed format 
			default: without reference gene annotation
			The absolute path to the annotated gene file
	-s/--seed	The size of the seed (e.g. anchor size)
			default: 20 (recommend values: 15-25; this value should be smaller than the (length of reads-10)/2)
	--mis		Maximum allowed mismatches during the search process (0, 1, or 2) 
			default: 0
	--min		Minimum distance between the two splice sites
			default: 100
	--max 		Maximum distance between the two splice sites
			default: 100000
	-o 		Output folder name 
			default: "AutoCirc_output"
	
	-h/--help	This help
	-v/--version	Print version

Examples :
	perl $0  -g /human_genome/hg19.fa -I /bowtie2_index/hg19 --bam unmapped.bam -b /hg19/Annotation/Refseq.gene.strandarded.bed --mis 0 --min 100 --max 100000 -s 20 -o autocirc_output 
	or
	perl $0  -g /human_genome/hg19.fa -I /bowtie2_index/hg19 --bam unmapped.bam -b /hg19/Annotation/Refseq.gene.strandarded.bed
	or
	perl $0  -g /human_genome/hg19.fa --bam unmapped.bam -s 15 -o autocirc_output (if do not have the reference gene annotations)   

NOTE: 
	1. The reference genome file should be the same one used in the previous Botwie/Bowtie2/Tophat mapping.
	2. The reference gene annotation should be in standard bed format with full 12 fields 
	   ( https://genome.ucsc.edu/FAQ/FAQformat.html#format1 )
EOF

sub main::HELP_MESSAGE {
	print $USAGE, $/;
	exit 1;
}

sub main::VERSION {
	print $VERSION,$/;
	exit 1;
}

my $unmapped_bam;
my $unmapped_sam;
my $genome_index;
my $ref_genome;
my $mismatch=0;
my $min=100;
my $max=100000;
my $ref_gene_bed="NA";
my $out_folder="AutoCirc_Output";
my $anchor_size=20;
my $help;
my $version;

GetOptions(
	    'bam=s' => \$unmapped_bam,
	    'gindex|I=s'=>\$genome_index,
	    'genome|g=s'=>\$ref_genome,	
	    'mis=s'=>\$mismatch,	
	    'min=s'=>\$min,
	    'max=s'=>\$max,
	    'bed|b=s'=>\$ref_gene_bed,
	    'seed|s=s'=>\$anchor_size,
	    'outfolder|o=s'=>\$out_folder,
	    'help|h' => \$help,
	    'version|v'=>\$version,
) or die $USAGE;

if ($help) {
	main::HELP_MESSAGE;
}

if ($version) {
	main::VERSION;
}

if ( (not $unmapped_bam) or (not $genome_index) ) {
	
	main::HELP_MESSAGE;
}

if ($mismatch == 0 or $mismatch == 1 or $mismatch == 2){ }
else{
	main::HELP_MESSAGE;
}

if ( ($ref_gene_bed ne "NA") and (not $ref_gene_bed =~/\S+\.bed/) ){
	main::HELP_MESSAGE;
} 

if( (not $ref_genome =~/^\/\S+/) and (not $ref_genome =~/^\~\/\S+/)){
	warn "\nPls provide the absolute path of a reference genome file in fasta format\n\n";
	main::HELP_MESSAGE;
}
if( (not $genome_index =~/^\/\S+/) and (not $genome_index =~/^\~\/\S+/)){
	warn "\nPls provide the absolute path of the genome index of bowtie2\n\n";
	main::HELP_MESSAGE;
}

## step 1 take in bowtie/tophat unmapped bam/sam file
#         extract the anchors from two ends of unmapped-reads

if ($anchor_size>20 or $anchor_size<15){
	main::HELP_MESSAGE;
}

if($unmapped_bam =~/(\S+)\.bam/){
	$unmapped_sam=join "\.", ($1, "sam");
	warn "Now convert the BAM file to SAM file:\n";
	system "mkdir $out_folder";	
	warn "samtools view -h $unmapped_bam -o ./$out_folder/$unmapped_sam\n";
	GENE_SAM: my $generate_sam=`samtools view -h $unmapped_bam -o ./$out_folder/$unmapped_sam`;

	until($generate_sam ne 0){ 
		print "......... \n";
		sleep(60);
	}
}else{
	main::HELP_MESSAGE;
}

chdir($out_folder);
warn "cd $out_folder\n";


warn "Now generate the two end anchors of unmapped reads in fastq file\n";
warn "../script/unmapped2anchors.pl $unmapped_sam $anchor_size > unmapped_anchors.qfa\n";
my $generate_anchors=`../script/unmapped2anchors.pl $unmapped_sam $anchor_size > unmapped_anchors.qfa`;

until($generate_anchors ne 0){
	print ".........\n";
	sleep(20);
}

##step 2: bowtie2 map the anchores seq. into ref. genome

my $map_anchors;
$map_anchors=`bowtie2 -p4 --reorder --mm -D20 --score-min=C,-15,0 -q -x $genome_index -U unmapped_anchors.qfa 1>bowtie2.anchors.out.sam 2>bowtie2.anchors.log`;
warn "Now mapping the anchors into the reference genome\n";
warn "bowtie2 -p4 --reorder --mm -D20 --score-min=C,-15,0 -q -x $genome_index -U unmapped_anchors.qfa 1> bowtie2.anchors.out.sam 2> bowtie2.anchors.log\n";

until($map_anchors ne 0){
	print ".........\n";
	sleep(60);
}

##step 3  -- detect all circ meeting the GT/AG rule

warn "Now detecting the potential circRNAs following GT/AG rules in the back splice sites\n";
warn "../script/autocirc_gtag.v1.1 bowtie2.anchors.out.sam  $ref_genome $mismatch $min $max $anchor_size circ.gtag.bed linear.gtag.bed backsplice.raw.bed 2> clog.AutoCirc_gtag\n";
my $generate_gtagcirc=`../script/autocirc_gtag.v1.1 bowtie2.anchors.out.sam  $ref_genome $mismatch $min $max $anchor_size circ.gtag.bed linear.gtag.bed backsplice.raw.bed 2> clog.AutoCirc_gtag`;

until($generate_gtagcirc ne 0){
	print ".........\n";
	sleep(60);
}

##step 4 -- detect circ meet the annotated exon boundaries
if($ref_gene_bed=~/\S+\.bed/){
	warn "Now detecting the potential circRNAs meeting the exon boundaries of annotated genes\n";
	warn "intersectBed -a backsplice.raw.bed -b $ref_gene_bed -f 1 -wa -wb | cut -f 1-13,15,19-22 > backwards.annotatedexons.bed\n";
	
	my $generate_circ_annotatedexons=`intersectBed -a backsplice.raw.bed -b $ref_gene_bed -f 1 -wa -wb | cut -f 1-13,15,19-22 > backwards.annotatedexons.bed`;
	until($generate_circ_annotatedexons ne 0){
		print ".........\n";
		sleep(10);
	}

	warn "../script/filter_exon_boundariesv1.2.pl backwards.annotatedexons.bed 10 > circ.exonboundaries.unstrand.bed\n";
	my $filter_exonboundaries=`../script/filter_exon_boundariesv1.2.pl backwards.annotatedexons.bed 10 > circ.exonboundaries.unstrand.bed`;
	until($filter_exonboundaries ne 0){
		print ".........\n";
		sleep(10);
	}

	##step 5 -- combine the circ identifed based on GT/AG rule and these found by annotated exon boundaries
	warn "Now combining the circRNAs detected by GT/AG rules or exon boundaries of annotated genes\n";
	warn "../script/combine.pl circ.gtag.bed circ.exonboundaries.unstrand.bed > circ.combined.bed\n";
	my $combine=`../script/combine.pl circ.gtag.bed circ.exonboundaries.unstrand.bed > circ.combined.bed`;
	until($combine ne 0){
		print ".........\n";
		sleep(5);
	}
 

	##step 6: remove pseudo-circ

	warn "Now excluding the pseudo circRNAs which across two duplicated genes or multiple genes\n";
	my $rem_1=`intersectBed -a circ.combined.bed -b $ref_gene_bed -s -wa -wb | cut -f 1-12,14-16,18-21 > circ.ref_gene.bed`;
	until($rem_1 ne 0){
		print ".........\n";
		sleep(3);
	}

	my $rem_2=`../script/circ_across_mgenes.pl circ.ref_gene.bed > circ.pseudo.bed 2>err`;
	until($rem_2 ne 0){
		print ".........\n";
		sleep(3);
	}

	my $rem_3=`../script/my_join.pl -a circ.combined.bed -b circ.pseudo.bed -F 4 -f 4 | awk '\$9==\"\"' | cut -f 1-8 > circ.final.bed`;
	until($rem_3 ne 0){
		print ".........\n";
		sleep(3);
	}
	
	system "mkdir tmp";
	system "mv backwards.annotatedexons.bed ./tmp";
	system "mv circ.pseudo.bed ./tmp";
	system "mv circ.combined.bed ./tmp";
	system "mv circ.exonboundaries.unstrand.bed ./tmp";
	system "rm *sam -rf";

}elsif($ref_gene_bed eq "NA"){ #do not input bed file 
	system "cp circ.gtag.bed circ.final.bed";
	warn "Without reference gene annotation, AutoCirc ab initio detect circRNAs based on GT/GT rules\n\n";
}
warn "Finished\n"
