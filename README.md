# AutoCirc


## What is AutoCirc
AutoCirc is a computational program to automatically identify back-splice junctions of potential circular RNAs (circRNAs) from RNA-seq data. Its core algorithm is implemented by C++, and it is wrapped in Perl, so it runs faster than most other circRNA detection tools under the same computer environment. AutoCirc is designed to detect circRNAs for all species as long as their genomes are available. It does not require gene annotations, but gene annotations will improve detection of circRNAs. AutoCirc runs on Linux/Unix OS.

Version: 1.3.1

Last Modified: 04/08/2020

Authors: Chan Zhou (chan.zhou@umassmed.edu)

Maintainer: Chan Zhou 

## A Workflow of AutoCirc 
![pipeline](https://github.com/chanzhou/AutoCirc/blob/master/AutoCirc_Workflow.jpg)

AutoCirc detects circRNAs from the BAM file of unmapped reads (output of Bowtie2 / Bowtie / TopHat) by the following steps: 
1)	Extract the 20 nucleotide anchors from the left and right ends of unmapped reads and map them to the reference genome with Bowtie2. 
2)	Scan the mapped anchors from both ends of unmapped reads to identify forward splice junctions (splicing from the 3’ end of an upstream sequence to the 5’ end of a downstream sequence) and back splice junctions (splicing from the 3’ end of a downstream sequence to the 5’ end of an upstream sequence) with unambiguous GT-AG breakpoint detection. We required a genomic distance between the two splice sites of no more than 100 kb and at least 100 bp (as default values). 
3)	Identify the back splice junctions that meet the annotated exon boundaries (if known gene annotation is provided).
4)	Merge all the back splice junctions either found based on the GT-AG rule of splice donor and acceptor or based on the exon boundaries of annotated genes. 
5)	Remove all forward splice junctions and filter out the pseudo-back splice junctions, which span the 3’ end and 5’ end of the duplicated genes.

Finally, AutoCirc will output the back splice junctions that contain the canonical GT-AG splice donor and acceptor sequences in the flanking intron boundaries or that occur at known exon boundaries (if inputting reference annotation, e.g. RefSeq used here). These back splice junctions are considered to represent circRNAs.

## Prerequisites

To use AutoCirc, you will need the following programs in your PATH:
* Bowtie2 (version 2.2.5 or higher)
* bedtools
* samtools

You will also need Perl version 5.16.3 or higher.


## RNA-seq data
PolyA-depleted, ribosome-depleted, or RNase R-treated RNA-seq data in FASTQ format.
The lengths of reads should be at least 50 nucleotides (nt). 
To identify reliable back splice junctions of potential circRNAs, we recommend using RNA-seq with read lengths of at least 80 nt . 

## Data Preparation
To obtain the unmapped reads, the user should map all reads to a reference genome either by *Bowtie2 or TopHat*. 
#### *For paired end RNA-seq data, each end of RNA-seq data will be processed separately as a single end RNA-seq data*.

For example:
* Bowtie2:

*Step 1: Mapping reads to reference genome*

```bash
bowtie2 -p4 --very-sensitive --mm -D20 --score-min=C,-15,0 -q -x <bowtie2_genome_index> -U <reads.fastq> 1>bowtie2out.sam 2>bowtie2.log 
```

*Step 2: extract the unmapped reads in BAM format*

```bash
samtools view -hbuS bowtie2out.sam | samtools view -hf 4 - | samtools view -Sb - | samtools sort - unmapped
```

### OR

* TopHat:
```bash
tophat -p 4  <options> -o tophat_out <bowtie2_genome_index> <reads.fastq>
```

You will use the *unmapped.bam* file obtained by Bowtie2 or TopHat (in the output folder of TopHat).

## Installation 

### Download AutoCirc
```bash
git clone https://github.com/chanzhou/AutoCirc.git 
cd AutoCirc
chmod 755 AutoCirc_v1.3.pl
chmod 755 script/*
```
Then run the pipeline inside the AutoCirc folder, which includes all the required scripts. 

### Usage
```bash
./AutoCirc_v1.3.pl  [options]

Options:	
	-g/--genome     Absolute path of the reference genome file (Fasta format) 
	-I/--gindex     Absolute path of the genome index of bowtie2 (e.g. /bowtie2_index/hg19)
	--bam		Bam file of unmapped reads (output by bowtie or bowtie2 or tophat)
	-b/--bed	Reference gene annotation (standard bed format) 
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
	
	-h/--help	This usage
	-v/--version	Print version
```

### Examples :
```bash
./AutoCirc_v1.3.pl  -g /human_genome/hg19.fa -I /bowtie2_index/hg19 --bam unmapped.bam -b /hg19/Annotation/refGene.bed --mis 0 --min 100 --max 100000 -s 20 -o autocirc_output 
```
or

```bash
./AutoCirc_v1.3.pl  -g /human_genome/hg19.fa -I /bowtie2_index/hg19 --bam unmapped.bam -b /hg19/Annotation/refGene.bed
```
or

```bash
./AutoCirc_v1.3.pl  -g /human_genome/hg19.fa --bam unmapped.bam -s 15 -o autocirc_output   
```

## NOTE: 
1. The reference genome file should be the same one used in the previous Botwie/Bowtie2/Tophat mapping.
2. The reference gene annotation should be in standard bed format with full 12 fields ( https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and can be obtained from the UCSC genome browser.
3. The minimum and maximum distance between the two splice sites are suggested to be set as the numbers close to the minimum and the maximum size of annotated genes in the studied genome.  The current default option“--min 100 --max 100000” is set for human genome. If you will study other species, pls change these settings accordingly. 


## Output
AutoCirc will output both the final predicted back splice junction file (circ.final.bed) and also some intermediate tmp files and log files. All are included in the specified output folder (default is “AutoCirc_output”). The output includes all posible candidate predicted back splice junctions, so please exclude the lowly expressed candidates with low expression levels by your own criteria. Usually, we only keep the circRNAs with at least two read counts after merging the results from replicates.

The Final output file “circ.final.bed” follows the standard BED format as follows:

| Field(column)| Description                                   |
| :-----------:| :---------------------------------------------|
| Chr	       | Chromosome                                    |
| Start	       | Start of junction                             |
| End	       | End of junction                               |
| Junc_name    | The name of back splice junction              |
| Coverage     | Number of reads mapping to the junction       |
| Site_strand  | Strand of junction                            |
| Distance     | Distance between the two ends of junction     |
| Read_tags    | All the tags of reads mapping to the junction |


## Citation

* Chan Zhou*, Benoit Molinie*, Kaveh Daneshvar, Joshua V. Pondick, Jinkai Wang, Nicholas O. Van Wittenberghe, Yi Xing, Cosmas C. Giallourakis+, and Alan C. Mullen+. Genome-Wide Maps of m6A circRNAs Identify Widespread and Cell-Type-Specific Methylation Patterns that Are Distinct from mRNAs. Cell Reports, 20(9), 2262–2276. doi:10.1016/j.celrep.2017.08.027 (published on Aug 29th, 2017) 

* bioRxiv (March, 2017) version:  http://biorxiv.org/content/early/2017/03/10/115899

## License
Copyright (C) 2017 The Mullen Lab. See the LICENSE file for license rights and limitations.


