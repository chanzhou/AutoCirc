#!/usr/bin/perl
#######################################################################
## @version 0.31
##
## Change log 
##     Add the regexp support for the split separator
##
## @version 0.3
##
## change log
##     Add the connector for multi-keyword connection
##
## 
## @version 0.2
##
## change log
##     change parameter parser to Getopt::Long modular
##     Add the debug mode
##     Add one-to-more hit
##     Get ride of the bug for no hit
##                           
##                                  Sat May 14 22:45:29 CST 2005
#######################################################################

use warnings;
use strict;
use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case bundling);

use vars qw($field_first $field_second $seperator_first $seperator_second $seperator_output $gap_num);
use vars qw( $file_first $file_second $seperator_connect $debug);

GetOptions(
	   "F=s"   => \$field_first,
	   "f=s"   => \$field_second,
           "S=s"   => \$seperator_first,
	   "s=s"   => \$seperator_second,
           "o=s"   => \$seperator_output,
	   "c=s"   => \$seperator_connect,
           "N=i"   => \$gap_num,
           "d"     => \$debug,
           "a=s"   => \$file_first,
	   "b=s"   => \$file_second
	  )   || Usage();

usage() unless (defined($file_first) && -f $file_first) && (defined($file_second) && -f $file_second);

my %index=();

$field_first="1" unless defined($field_first);
$field_second="1" unless defined($field_second);

$seperator_first="\t" unless defined($seperator_first);
$seperator_first=~s/^\"|\"$//g;

$seperator_second="\t" unless defined($seperator_second);
$seperator_second=~s/^\"|\"$//g;

$seperator_output="\t" unless defined($seperator_output);
$seperator_output=~s/^\"|\"$//g;

$seperator_connect="\t" unless defined($seperator_connect);
$seperator_connect=~s/^\"|\"$//g;

$gap_num=1 unless $gap_num;

open(FILE1,$file_first) or die "cannot open the file $file_first";
open(FILE2,$file_second) or die "cannot open the file $file_second";

my @keya=split(/\D+/,$field_first);
my @keyb=split(/\D+/,$field_second);

while (my $line=<FILE2>) {
  chomp($line);
  my @items=split(/$seperator_second/,$line);

  my $key_file2="";

  for (my $i=0;$i<@keyb;$i++) {
    $items[$keyb[$i]-1] = "" unless defined($items[$keyb[$i]-1]);
    $key_file2.=$items[$keyb[$i]-1];
    $key_file2.="$seperator_connect" if $i<@keyb-1;
  }

  print ">>> key2: [$key_file2] sep2: [$seperator_second]\n" if $debug;

  push(@{$index{$key_file2}},$line) if $key_file2 ne "";
}

while (my $line=<FILE1>) {
  chomp($line);
  my @items=split(/$seperator_first/,$line);

  my $key_file1="";

  for (my $i=0;$i<@keya;$i++) {
    $items[$keya[$i]-1]="" unless defined($items[$keya[$i]-1]);
    $key_file1.=$items[$keya[$i]-1];
    $key_file1.="$seperator_connect" if $i<@keya-1;
  }

  print "<<< key1: [$key_file1] sep1: [$seperator_first]\n" if $debug;

  if ($key_file1 ne "" && exists($index{$key_file1})) {

    for (my $i=0;$i<@{$index{$key_file1}};$i++) {
      print "$line$seperator_output";
      print $index{$key_file1}->[$i]."\n";
    }
  } else {
    print "$line$seperator_output";
    print "$seperator_output*-" x ($gap_num) if $gap_num>1;
    print "\n";
  }

}
close(FILE1);
close(FILE2);

exit(0);
################################################
sub usage
  {
    print "This program join two file together used the common key\n";
    print "Usage:\n\n";
    print "$0 [-d] [-F 1,2] [-f 1,2] [-S \"\\t\"] [-s \"\\t\"] [-c \"\\t\"] [-N num] [-o \"\\t\"] <-a file1> <-b file2>\n";
    print "  d: debug mode\n";
    print "  F: the field group of key in file1, separate by ','\n";
    print "  f: the field group of key in file2, separate by ','\n";
    print "  c: the connector of keylist, used when connect with many keywords [default :\\t]\n";
    print "  S: the seperator of file1. WARN: regexp support warning \\, quote with ', [default : \\t]\n";
    print "  s: the seperator of file2. WARN: regexp support warning \\, quote with '. [default : \\t]\n";
    print "  o: the connector for two file [default : \\t]\n";
    print "  N: the number of \\t if lost data in file2.default:1\n";
    exit(0);
  }
